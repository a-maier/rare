use std::{
    num::NonZeroUsize,
    ops::{AddAssign, ControlFlow},
};

use galois_fields::Z64;
use log::{debug, trace};
use paste::paste;
use rand::Rng;

use crate::{
    dense_poly::{DensePoly, DensePoly1},
    rat::Rat,
    rec_linear::LinearRec,
    rec_thiele::ThieleRec,
    sparse_poly::{SparseMono, SparsePoly},
    traits::{Eval, One, Rec, Shift, TryEval, WithVars, Zero},
};

/// Rational function reconstruction
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct RatRecMod {
    extra_pts: usize,
}

trait RecWithRanAndShift<const N: usize> {
    fn rec_with_ran_and_shift<F, const P: u64>(
        &self,
        f: F,
        rng: impl Rng,
        shift: [Z64<P>; N],
    ) -> Option<Rat<SparsePoly<Z64<P>, N>>>
    where
        F: FnMut([Z64<P>; N]) -> Option<Z64<P>>;
}

impl RatRecMod {
    pub fn new(extra_pts: usize) -> Self {
        Self { extra_pts }
    }
}

macro_rules! impl_rec_with_ran_and_shift {
    ( $($x:literal, $y:literal), * ) => {
        $(
            paste! {
                impl RecWithRanAndShift<$x> for RatRecMod {
                    fn rec_with_ran_and_shift<F, const P: u64>(
                        &self,
                        mut f: F,
                        mut rng: impl Rng,
                        shift: [Z64<P>; $x],
                    ) -> Option<Rat<SparsePoly<Z64<P>, $x>>>
                    where F: FnMut([Z64<P>; $x]) -> Option<Z64<P>>
                    {
                        debug!("Shift: {shift:?}");

                        // initial [x0, x1, ..., t]
                        let base_coord: [Z64<P>; $x] = rng.gen();
                        debug!("Starting multivariate rational reconstruction with base point {base_coord:?}");

                        // varying t while keeping x0, ... fixed
                        let coord = transformed_coord_iter(base_coord, shift);

                        // reconstruct in t
                        let thiele_rec = ThieleRec::new(self.extra_pts);
                        let rat = thiele_rec.rec_from_seq(
                            coord.filter_map(|(t, z)| f(z).map(|v| (t, v)))
                        )?;
                        let rat: Rat<_> = rat.into();
                        debug!("From Thiele reconstruction: {}", rat.with_vars(&["t"]));
                        debug_assert!(!rat.den().coeff(0).is_zero());

                        // first point to be used for starting the reconstruction
                        let pt = RatPt::<P, $y>::new([Z64::zero(); $y], rat);

                        let mut rec = [<RecHelper $x>]::new(pt, self.extra_pts, base_coord, shift);
                        use NumOrDen::*;
                        let num_rec = rec.rec(|x| f(x), Num)?;
                        debug!("Reconstructed numerator: {num_rec}");
                        let den_rec = rec.rec(|x| f(x), Den)?;
                        debug!("Reconstructed denominator: {den_rec}");
                        let [mut num_rec, mut den_rec] = rec.into_rec();

                        // fix normalisation
                        let norm = den_rec.term(0).coeff.inv();
                        num_rec *= &norm;
                        den_rec *= &norm;
                        Some(Rat::from_num_den_unchecked(num_rec, den_rec))
                    }
                }
            }
        )*
    };
}

// transform from [x0, x1, ..., t] to
// [z0, z1, z2, ...] = [t*x0 + s0, t*x1 + s1, t*x2 + s2, ..., t + sM]
// with a `shift` = [s0, s1, s2, ..., sM]
fn x_to_z<const N: usize, const P: u64>(
    mut coord: [Z64<P>; N],
    shift: &[Z64<P>; N],
) -> [Z64<P>; N] {
    let (t, xs) = coord.split_last_mut().unwrap();
    for x in xs {
        *x *= *t
    }
    for (z, s) in coord.iter_mut().zip(shift.iter()) {
        *z += s;
    }
    coord
}

// iterate over coordinates [x0, x1, ..., t] with varying values of t
// transforming to z
// with a shift s = [s0, s1, ...]
// returns an iterator over tuples (t, z)
fn transformed_coord_iter<const N: usize, const P: u64>(
    mut xcoord: [Z64<P>; N],
    shift: [Z64<P>; N],
) -> impl Iterator<Item = (Z64<P>, [Z64<P>; N])> {
    std::iter::repeat_with(move || {
        xcoord[N - 1] += Z64::one();
        (xcoord[N - 1], x_to_z(xcoord, &shift))
    })
    .take(P as usize)
}

// find a shift such that the shifted function fs(z) = f(z + s)
// has a constant leading term in the denominator
// In other words, find s such that f(s) is well-defined
fn find_shift<F, const P: u64, const N: usize>(
    mut f: F,
    mut rng: impl Rng,
) -> [Z64<P>; N]
where
    F: FnMut([Z64<P>; N]) -> Option<Z64<P>>,
{
    let mut shift = [Z64::zero(); N];
    for i in (0..N).cycle() {
        if f(shift).is_some() {
            return shift;
        }
        shift[i] = rng.gen()
    }
    unreachable!()
}

macro_rules! impl_rec_with_ran {
    ( $($x:literal), *) => {
        $(

            impl<F, const P: u64> Rec<RatRecMod, [Z64<P>; $x]> for F
            where F: FnMut([Z64<P>; $x]) -> Option<Z64<P>> {
                type Output = Option<Rat<SparsePoly<Z64<P>, $x>>>;

                fn rec_with_ran(
                    &mut self,
                    rec: RatRecMod,
                    mut rng: impl Rng
                ) -> Self::Output {
                    let shift = find_shift(|z| (self)(z), &mut rng);
                    rec.rec_with_ran_and_shift(
                        |z| (self)(z),
                        rng,
                        shift
                    )
                }
            }

            impl<F, const P: u64> Rec<RatRecMod, [[Z64<P>; $x]; 1]> for F
            where F: TryEval<[Z64<P>; $x], Output = Z64<P>> {
                type Output = Option<Rat<SparsePoly<Z64<P>, $x>>>;

                fn rec_with_ran(
                    &mut self,
                    rec: RatRecMod,
                    rng: impl Rng
                ) -> Self::Output {
                    (|x: [Z64<P>; $x]| self.try_eval(&x)).rec_with_ran(rec, rng)
                }
            }
        )*
    };
}

impl<F, const P: u64> Rec<RatRecMod, [Z64<P>; 1]> for F
where F: FnMut([Z64<P>; 1]) -> Option<Z64<P>> {
    type Output = Option<Rat<SparsePoly<Z64<P>, 1>>>;

    fn rec_with_ran(
        &mut self,
        rec: RatRecMod,
        rng: impl Rng
    ) -> Self::Output {
        let rec = ThieleRec::new(rec.extra_pts);
        rec.rec_univariate_with_ran(|x: Z64<P>| (self)([x]), rng)
            .map(|r| r.into())
    }
}

impl<F, const P: u64> Rec<RatRecMod, [[Z64<P>; 1]; 1]> for F
where F: TryEval<[Z64<P>; 1], Output = Z64<P>> {
    type Output = Option<Rat<SparsePoly<Z64<P>, 1>>>;

    fn rec_with_ran(
        &mut self,
        rec: RatRecMod,
        rng: impl Rng
    ) -> Self::Output {
        (|x: [Z64<P>; 1]| self.try_eval(&x)).rec_with_ran(rec, rng)
    }
}


#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
enum NumOrDen {
    Num,
    Den,
}

macro_rules! impl_rec_helper {
    ( $($x:literal, $y:literal), * ) => {
        $(
            paste! {
                use crate::dense_poly::[<DensePoly $x>];
                use crate::rec_newton::[<NewtonPolyRec $y>];

                #[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
                struct [<RecHelper $x>]<const P: u64> {
                    known_pts: Vec<RatPt<P, $y>>,
                    extra_pts: usize,
                    base_coord: [Z64<P>; $x],
                    shift: [Z64<P>; $x],
                    num: SparsePoly<Z64<P>, $x>,
                    den: SparsePoly<Z64<P>, $x>,
                    num_coeff_known: usize,
                    den_coeff_known: usize,
                }

                impl<const P: u64> [<RecHelper $x>]<P> {
                    fn new(
                        first_pt: RatPt<P, $y>,
                        extra_pts: usize,
                        base_coord: [Z64<P>; $x],
                        shift: [Z64<P>; $x],
                    ) -> Self {
                        Self {
                            known_pts: vec![first_pt],
                            extra_pts,
                            base_coord,
                            shift,
                            ..Default::default()
                        }
                    }

                    fn rec<F>(
                        &mut self,
                        mut f: F,
                        wot: NumOrDen
                    ) -> Option<&SparsePoly<Z64<P>, $x>>
                    where F: FnMut([Z64<P>; $x]) -> Option<Z64<P>>
                    {
                        debug!("Reconstruct {wot:?} coefficients");
                        use NumOrDen::*;
                        let base_x_coord: [Z64<P>; $y] = self.base_coord[..$y].try_into().unwrap();

                        let nnum = self.known_pts[0].num_coeff.len();
                        let nden = self.known_pts[0].den_coeff.len();
                        // on to the meat:
                        // iteratively reconstruct the highest-degree polynomial in the numerator / denominator
                        // storing also the results for all other coefficients for later use
                        let max_idx = match wot {
                            Num => nnum,
                            Den => usize::from(nden)
                        };
                        for coeff_idx in (0..max_idx).rev() {
                            debug!("Reconstruct coefficient of t^{coeff_idx}");
                            let mut poly = [<NewtonPolyRec $y>]::<P>::new(self.extra_pts);
                            // what do we know at this point?
                            // TODO: maybe use a BTreeMap so we don't have to sort?
                            self.known_pts.sort_unstable();
                            let known = self.known_pts.iter();
                            let known: Vec<_> = match wot {
                                Num => known.map(|pt| (pt.offset, pt.num_coeff[coeff_idx])).collect(),
                                Den => known.map(|pt| (pt.offset, pt.den_coeff[coeff_idx])).collect(),
                            };
                            let mut known = known.iter().peekable();

                            let first_pt = known.next().unwrap();
                            let mut next_offset = first_pt.0;
                            let mut next_val = first_pt.1;
                            loop {
                                let coord = add(base_x_coord, next_offset);
                                match poly.add_pt(&coord, next_val) {
                                    ControlFlow::Continue(n) => {
                                        next_offset[n] += Z64::one();
                                        next_offset[(n + 1)..].fill(Z64::zero());
                                        // skip known points with a lesser offset
                                        // they are not needed for this reconstruction
                                        // TODO: some way to do this with `skip_while`?
                                        while let Some((offset, _)) = known.peek() {
                                            if offset < &next_offset {
                                                known.next();
                                            } else {
                                                break;
                                            }
                                        }
                                        // TODO: accept points with larger offset in the same coordinate
                                        match known.peek() {
                                            Some((offset, val)) if offset == &next_offset => {
                                                trace!("Using known value at offset {next_offset:?}: {val}");
                                                next_val = *val;
                                            },
                                            _ => {
                                                trace!("Need new value at offset {next_offset:?}");
                                                let mut coord = self.base_coord;
                                                for (x, s) in coord[..$y].iter_mut().zip(next_offset.iter()) {
                                                    *x += s;
                                                }

                                                let eval_pts = transformed_coord_iter(coord, self.shift);
                                                let linear_rec = LinearRec::new(
                                                    nnum - self.num_coeff_known,
                                                    NonZeroUsize::new(nden - self.den_coeff_known).unwrap()
                                                );
                                                trace!("Known numerator: {}", self.num);
                                                let rat = linear_rec.rec_from_seq_with_subtr(
                                                    eval_pts.filter_map(|(t, z)| {
                                                        f(z).map(
                                                            |f_z| {
                                                                let sub = self.num.eval(&z) - f_z * self.den.eval(&z);
                                                                trace!("Subtraction for {z:?}: {sub}");
                                                                (t, f_z, sub)
                                                            }
                                                        )
                                                    }),
                                                )?;
                                                trace!("From rational reconstruction in t: {}", rat.with_vars(&["t"]));

                                                next_val = match wot {
                                                    Num => *rat.num().coeff(coeff_idx),
                                                    Den => *rat.den().coeff(coeff_idx),
                                                };
                                                let pt = RatPt::<P, $y>::new(next_offset, rat);

                                                self.known_pts.push(pt);
                                            },
                                        }

                                    },
                                    ControlFlow::Break(()) => {
                                        let rec_poly: [<DensePoly $y>]<Z64<P>> = poly.into_poly().into();
                                        // turn it back from a polynomial in [x0, x1, ..., t] to [z0, z1, ...]
                                        // first expand it out
                                        let rec_poly: SparsePoly<Z64<P>, $y> = rec_poly.into();
                                        // rescale, adding back the power in t
                                        let terms = rec_poly.into_terms().into_iter()
                                            .map(|c| {
                                                let mut powers = [0; $x];
                                                powers[0..$y].copy_from_slice(&c.powers);
                                                powers[$y] = coeff_idx as u32 - c.powers.iter().sum::<u32>();
                                                SparseMono::new(c.coeff, powers)
                                            });
                                        let rec_poly = SparsePoly::from_terms(terms.collect());
                                        debug!("Reconstructed numerator contribution {rec_poly}");
                                        if !self.shift.is_zero() {
                                            self.add_shift_subtr_from(wot, rec_poly.clone());
                                        }
                                        match wot {
                                            Num => {
                                                self.num += rec_poly;
                                                self.num_coeff_known += 1
                                            },
                                            Den => {
                                                self.den += rec_poly;
                                                self.den_coeff_known += 1
                                            },
                                        };
                                        break;
                                    },
                                }
                            }
                        }
                        match wot {
                            NumOrDen::Num => {
                                Some(&self.num)
                            },
                            NumOrDen::Den => {
                                Some(&self.den)
                            },
                        }
                    }

                    fn add_shift_subtr_from(
                        &mut self,
                        wot: NumOrDen,
                        poly: SparsePoly<Z64<P>, $x>
                    ) {
                        // calculate P(z + s) - P(Z)
                        let collected = [<DensePoly $x>]::from(poly);
                        let shifted = collected.clone().shift(self.shift);
                        let subtr = shifted - &collected;
                        trace!("shift subtraction in z: {subtr}");
                        // turn it into a polynomial in [t, x0, x1, ...]
                        let subtr: SparsePoly<Z64<P>, $x> = subtr.into();
                        let terms = subtr.into_terms().into_iter()
                            .map(|c| {
                                let mut powers = [0; $x];
                                powers[0] = c.powers.iter().sum::<u32>();
                                powers[1..].copy_from_slice(&c.powers[0..$y]);
                                SparseMono::new(c.coeff, powers)
                            });
                        let subtr = SparsePoly::from_terms(terms.collect());
                        let subtr = [<DensePoly $x>]::from(subtr);
                        trace!("shift subtraction in [t, x1, ...]: {subtr}");
                        match wot {
                            NumOrDen::Num => self.subtr_from_num(subtr),
                            NumOrDen::Den => self.subtr_from_den(subtr),
                        }
                    }

                    pub(crate) fn into_rec(self) -> [SparsePoly<Z64<P>, $x>; 2] {
                        [self.num, self.den]
                    }

                    fn subtr_from_num(&mut self, poly: [<DensePoly $x>]<Z64<P>>) {
                        let base_x_coord: [Z64<P>; $y] = self.base_coord[..$y].try_into().unwrap();
                        for pt in &mut self.known_pts {
                            let coord = add(base_x_coord, pt.offset);
                            for coeff_idx in 0..poly.len() {
                                pt.num_coeff[coeff_idx] -= poly.coeff(coeff_idx).eval(&coord);
                            }
                        }
                    }

                    fn subtr_from_den(&mut self, poly: [<DensePoly $x>]<Z64<P>>) {
                        // TODO: code duplication
                        let base_x_coord: [Z64<P>; $y] = self.base_coord[..$y].try_into().unwrap();
                        for pt in &mut self.known_pts {
                            let coord = add(base_x_coord, pt.offset);
                            for coeff_idx in 0..poly.len() {
                                pt.den_coeff[coeff_idx] -= poly.coeff(coeff_idx).eval(&coord);
                            }
                        }
                    }

                }
            }
        )*
    };
}

fn add<T: AddAssign, const N: usize>(mut a: [T; N], b: [T; N]) -> [T; N] {
    for (a, b) in a.iter_mut().zip(b.into_iter()) {
        *a += b;
    }
    a
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct RatPt<const P: u64, const N: usize> {
    offset: [Z64<P>; N],
    num_coeff: Vec<Z64<P>>,
    den_coeff: Vec<Z64<P>>,
}

impl<const P: u64, const N: usize> RatPt<P, N> {
    fn new(pt: [Z64<P>; N], rat: Rat<DensePoly<Z64<P>>>) -> Self {
        let (num, den) = rat.into_num_den();
        Self {
            offset: pt,
            num_coeff: num.into_coeff(),
            den_coeff: den.into_coeff(),
        }
    }
}

impl_rec_with_ran!(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
impl_rec_with_ran_and_shift!(
    16, 15, 15, 14, 14, 13, 13, 12, 12, 11, 11, 10, 10, 9, 9, 8, 8, 7, 7, 6, 6,
    5, 5, 4, 4, 3, 3, 2, 2, 1
);
impl_rec_helper!(
    16, 15, 15, 14, 14, 13, 13, 12, 12, 11, 11, 10, 10, 9, 9, 8, 8, 7, 7, 6, 6,
    5, 5, 4, 4, 3, 3, 2, 2, 1
);

#[cfg(test)]
mod tests {
    use std::iter::repeat_with;

    use rand::Rng;
    use rand_xoshiro::rand_core::SeedableRng;

    use crate::{
        sparse_poly::SparseMono,
        traits::{One, TryEval, Zero},
    };

    use super::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn rec_rat2_small_noshift() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_POW: u32 = 2;
        const MAX_COEFF: u64 = 29;
        const P: u64 = 1152921504606846883;
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = RatRecMod::new(1);

        for _ in 0..NTESTS {
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let coeff = repeat_with(|| {
                SparseMono::new(
                    Z64::<P>::new_unchecked(rng.gen_range(0..MAX_COEFF)),
                    [(); 2].map(|_| rng.gen_range(0..=MAX_POW)),
                )
            })
            .take(nterms)
            .collect();
            let num = SparsePoly::from_terms(coeff);

            let den = if num.is_zero() {
                One::one()
            } else {
                let mut den: SparsePoly<_, 2> = Zero::zero();
                while den.is_zero() {
                    let coeff = repeat_with(|| {
                        SparseMono::new(
                            Z64::new_unchecked(rng.gen_range(0..MAX_COEFF)),
                            [(); 2].map(|_| rng.gen_range(0..=MAX_POW)),
                        )
                    })
                    .take(nterms)
                    .collect();
                    den = SparsePoly::from_terms(coeff);
                }
                // ensure we don't need a shift
                let mut den = den.into_terms();
                den[0].powers = Zero::zero();
                den[0].coeff = One::one();
                SparsePoly::from_terms(den)
            };

            let rat = Rat::from_num_den_unchecked(num, den);
            eprintln!("trying to reconstruct {rat}");

            let res =
                (|z| rat.try_eval(&z)).rec_with_ran(rec, &mut rng).unwrap();
            eprintln!("reconstructed {res}");
            assert_eq!(rat, res);
        }
    }

    #[test]
    fn rec_rat2_simple() {
        log_init();

        const P: u64 = 1152921504606846883;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = RatRecMod::new(1);

        let x = SparseMono::new(Z64::one(), [1, 0]);
        let y = SparseMono::new(Z64::one(), [0, 1]);

        // 1 / y
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![y]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let res = (|z: [Z64<P>; 2]| rat.try_eval(&z))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        eprintln!("reconstructed {res}");
        assert_eq!(res, rat);

        // 1 / y^2
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![y * y]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let res = (|z: [Z64<P>; 2]| rat.try_eval(&z))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        eprintln!("reconstructed {res}");
        assert_eq!(res, rat);

        // 1 / x
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![x]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let res = (|z: [Z64<P>; 2]| rat.try_eval(&z))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        eprintln!("reconstructed {res}");
        assert_eq!(res, rat);
    }

    #[test]
    fn rec_rat2_small_shift() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_POW: u32 = 2;
        const MAX_COEFF: u64 = 29;
        const P: u64 = 1152921504606846883;
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = RatRecMod::new(1);

        for _ in 0..NTESTS {
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let coeff = repeat_with(|| {
                SparseMono::new(
                    Z64::<P>::new_unchecked(rng.gen_range(0..MAX_COEFF)),
                    [(); 2].map(|_| rng.gen_range(0..=MAX_POW)),
                )
            })
            .take(nterms)
            .collect();
            let mut num = SparsePoly::from_terms(coeff);

            let mut den = if num.is_zero() {
                One::one()
            } else {
                let mut den: SparsePoly<_, 2> = Zero::zero();
                while den.is_zero() {
                    let coeff = repeat_with(|| {
                        SparseMono::new(
                            Z64::new_unchecked(rng.gen_range(0..MAX_COEFF)),
                            [(); 2].map(|_| rng.gen_range(0..=MAX_POW)),
                        )
                    })
                    .take(nterms)
                    .collect();
                    den = SparsePoly::from_terms(coeff);
                }
                den
            };

            // normalise
            let norm = den.term(0).coeff.inv();
            num *= norm;
            den *= norm;

            let rat = Rat::from_num_den_unchecked(num, den);
            eprintln!("trying to reconstruct {rat}");

            let res =
                (|z| rat.try_eval(&z)).rec_with_ran(rec, &mut rng).unwrap();
            eprintln!("reconstructed {res}");
            assert!(sample_eq(&rat, &res, &mut rng));
        }
    }

    #[test]
    fn rec_rat3() {
        log_init();

        const NTESTS: u32 = 40;
        const MAX_POW: u32 = 2;
        const MAX_COEFF: u64 = 29;
        const P: u64 = 1152921504606846883;
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = RatRecMod::new(1);

        for _ in 0..NTESTS {
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let coeff = repeat_with(|| {
                SparseMono::new(
                    Z64::<P>::new_unchecked(rng.gen_range(0..MAX_COEFF)),
                    [(); 3].map(|_| rng.gen_range(0..=MAX_POW)),
                )
            })
            .take(nterms)
            .collect();
            let mut num = SparsePoly::from_terms(coeff);

            let mut den = if num.is_zero() {
                One::one()
            } else {
                let mut den: SparsePoly<_, 3> = Zero::zero();
                while den.is_zero() {
                    let coeff = repeat_with(|| {
                        SparseMono::new(
                            Z64::new_unchecked(rng.gen_range(0..MAX_COEFF)),
                            [(); 3].map(|_| rng.gen_range(0..=MAX_POW)),
                        )
                    })
                    .take(nterms)
                    .collect();
                    den = SparsePoly::from_terms(coeff);
                }
                den
            };

            // normalise
            let norm = den.term(0).coeff.inv();
            num *= norm;
            den *= norm;

            let rat = Rat::from_num_den_unchecked(num, den);
            eprintln!("trying to reconstruct {rat}");

            let res =
                (|z| rat.try_eval(&z)).rec_with_ran(rec, &mut rng).unwrap();
            eprintln!("reconstructed {res}");
            assert!(sample_eq3(&rat, &res, &mut rng));
        }
    }

    fn sample_eq<const P: u64>(
        orig: &Rat<SparsePoly<Z64<P>, 2>>,
        rec: &Rat<SparsePoly<Z64<P>, 2>>,
        mut rng: impl Rng,
    ) -> bool {
        const TESTS: usize = 10;
        for _ in 0..TESTS {
            let pt: [Z64<P>; 2] = rng.gen();
            if orig.try_eval(&pt) != rec.try_eval(&pt) {
                return false;
            }
        }
        true
    }

    fn sample_eq3<const P: u64>(
        orig: &Rat<SparsePoly<Z64<P>, 3>>,
        rec: &Rat<SparsePoly<Z64<P>, 3>>,
        mut rng: impl Rng,
    ) -> bool {
        const TESTS: usize = 10;
        for _ in 0..TESTS {
            let pt: [Z64<P>; 3] = rng.gen();
            if orig.try_eval(&pt) != rec.try_eval(&pt) {
                return false;
            }
        }
        true
    }
}
