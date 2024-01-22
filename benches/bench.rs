use std::iter::repeat_with;

use criterion::{criterion_group, criterion_main, Criterion};
use galois_fields::Z64;
use rand::Rng;
use rand_xoshiro::rand_core::SeedableRng;
use rare::{
    dense_poly::{DensePoly, DensePoly2},
    rat::Rat,
    rec_linear::LinearRec,
    rec_newton::{NewtonPoly, NewtonPoly2, NewtonRec},
    rec_thiele::{ThieleRat, ThieleRec},
    rec_rat_mod::RatRecMod,
    sparse_poly::SparsePoly,
    traits::{Eval, Rec, TryEval},
    _test_util::{gen_dense_poly1, gen_dense_poly2, gen_dense_rat1, gen_dense_rat2, gen_sparse_rat}, rec_linear_multivar::rec_coeff
};

fn rec_poly1<const P: u64>(
    polys: &[DensePoly<Z64<P>>],
) -> Vec<NewtonPoly<Z64<P>>> {
    let mut res = Vec::with_capacity(polys.len());
    let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
    let rec = NewtonRec::new(1);

    for poly in polys {
        let p = (|x: Z64<P>| poly.eval(&x))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        res.push(p);
    }
    res
}

fn rec_poly1_expanded<const P: u64>(
    polys: &[DensePoly<Z64<P>>],
) -> Vec<DensePoly<Z64<P>>> {
    let mut res = Vec::with_capacity(polys.len());
    let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
    let rec = NewtonRec::new(1);

    for poly in polys {
        let p = (|x: Z64<P>| poly.eval(&x))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        res.push(p.into());
    }
    res
}

fn rec_poly2<const P: u64>(
    polys: &[DensePoly2<Z64<P>>],
) -> Vec<NewtonPoly2<Z64<P>>> {
    let mut res = Vec::with_capacity(polys.len());
    let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
    let rec = NewtonRec::new(1);

    for poly in polys {
        let p = (|x| poly.eval(&x)).rec_with_ran(rec, &mut rng).unwrap();
        res.push(p);
    }
    res
}

fn rec_poly2_expanded<const P: u64>(
    polys: &[DensePoly2<Z64<P>>],
) -> Vec<DensePoly2<Z64<P>>> {
    let mut res = Vec::with_capacity(polys.len());
    let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
    let rec = NewtonRec::new(1);

    for poly in polys {
        let p = (|x| poly.eval(&x)).rec_with_ran(rec, &mut rng).unwrap();
        res.push(p.into());
    }
    res
}

fn rec_rat1<const P: u64>(
    rats: &[Rat<DensePoly<Z64<P>>>],
) -> Vec<ThieleRat<Z64<P>>> {
    let mut res = Vec::with_capacity(rats.len());
    let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
    let rec = ThieleRec::new(1);

    for poly in rats {
        let p = (|x: Z64<P>| poly.try_eval(&x))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        res.push(p);
    }
    res
}

fn rec_rat1_expanded<const P: u64>(
    rats: &[Rat<DensePoly<Z64<P>>>],
) -> Vec<Rat<DensePoly<Z64<P>>>> {
    let mut res = Vec::with_capacity(rats.len());
    let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
    let rec = ThieleRec::new(1);

    for rat in rats {
        let p = (|x: Z64<P>| rat.try_eval(&x))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        res.push(p.into());
    }
    res
}

fn rec_rat1_linear<const P: u64>(
    rats: &[Rat<DensePoly<Z64<P>>>],
) -> Vec<Rat<DensePoly<Z64<P>>>> {
    let mut res = Vec::with_capacity(rats.len());
    let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

    for rat in rats {
        let rec = LinearRec::new(
            rat.num().len(),
            rat.den().len().try_into().unwrap(),
        );
        let p = (|x: Z64<P>| rat.try_eval(&x))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        res.push(p);
    }
    res
}

fn rec_rat2<const P: u64>(
    rats: &[Rat<DensePoly2<Z64<P>>>],
) -> Vec<Rat<SparsePoly<Z64<P>, 2>>> {
    let mut res = Vec::with_capacity(rats.len());
    let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
    let rec = RatRecMod::new(1);

    for rat in rats {
        let p = (|x: [Z64<P>; 2]| rat.try_eval(&x))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        res.push(p);
    }
    res
}

fn rec_rat2_linear<const P: u64>(
    rats: &[Rat<SparsePoly<Z64<P>, 2>>],
) -> Vec<Rat<SparsePoly<Z64<P>, 2>>> {
    let mut res = Vec::with_capacity(rats.len());
    let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

    for rat in rats {
        let nneeded = rat.num().len() + rat.den().len() - 1;
        let pts = Vec::from_iter(
            repeat_with(|| {
                let pt = [(); 2].map(|_| rng.gen());
                if let Some(val) = rat.try_eval(&pt) {
                    Some((pt, val))
                } else {
                    None
                }
            }).filter_map(|pt| pt)
            .take(nneeded)
        );
        res.push(rec_coeff(&rat, &pts).unwrap());
    }
    res
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

    const N: u32 = 8;
    const NPOLYS: usize = 32;
    const P: u64 = 1152921504606846883;
    let mut polys: [DensePoly<Z64<P>>; NPOLYS] = Default::default();
    for poly in &mut polys {
        *poly = gen_dense_poly1(&[N], &mut rng)
    }
    c.bench_function("poly1", |b| b.iter(|| rec_poly1(&polys)));

    c.bench_function("poly1 expanded", |b| {
        b.iter(|| rec_poly1_expanded(&polys))
    });

    let mut polys: [DensePoly2<Z64<P>>; NPOLYS] = Default::default();
    for poly in &mut polys {
        *poly = gen_dense_poly2(&[N / 2, N / 2], &mut rng)
    }

    c.bench_function("poly2", |b| b.iter(|| rec_poly2(&polys)));

    c.bench_function("poly2 expanded", |b| {
        b.iter(|| rec_poly2_expanded(&polys))
    });
    std::mem::drop(polys);

    let mut rats: [Rat<DensePoly<Z64<P>>>; NPOLYS] = Default::default();
    for rat in &mut rats {
        *rat = gen_dense_rat1(&[N / 2], &mut rng)
    }

    c.bench_function("rat1 thiele", |b| b.iter(|| rec_rat1(&rats)));

    c.bench_function("rat1 thiele expanded", |b| {
        b.iter(|| rec_rat1_expanded(&rats))
    });

    c.bench_function("rat1 linear", |b| b.iter(|| rec_rat1_linear(&rats)));

    let mut rats: [Rat<DensePoly2<Z64<P>>>; NPOLYS] = Default::default();
    for rat in &mut rats {
        *rat = gen_dense_rat2(&[N / 2, N / 2], &mut rng)
    }

    c.bench_function("rat2", |b| b.iter(|| rec_rat2(&rats)));

    let mut rats: [Rat<SparsePoly<Z64<P>, 2>>; NPOLYS] = Default::default();
    for rat in &mut rats {
        *rat = gen_sparse_rat(N / 2, N / 2, &mut rng)
    }

    c.bench_function(
        "rat2 (known degrees)",
        |b| b.iter(|| rec_rat2_linear(&rats))
    );
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
