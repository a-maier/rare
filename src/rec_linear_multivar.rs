use galois_fields::Z64;

use crate::{
    matrix::Matrix,
    rat::Rat,
    rec_linear::gauss_solve,
    sparse_poly::{SparseMono, SparsePoly},
    traits::{Eval, One},
};

/// Reconstruct coefficients of a multivariate rational function over
/// a finite field from known points
// TODO: make a trait, impl for Rat<DensePoly>
pub fn rec_coeff<T, const P: u64, const N: usize>(
    rat: &Rat<SparsePoly<T, N>>,
    pts: &[([Z64<P>; N], Z64<P>)],
) -> Option<Rat<SparsePoly<Z64<P>, N>>> {
    assert!(rat.den().len() > 0);
    // Number of unknown coefficients - one of them is set to 1 to fix the normalisation
    let num_unknown = rat.num().len() + rat.den().len() - 1;
    // Need at least one equation for each coefficient in rat
    if pts.len() < num_unknown {
        return None;
    }

    // construct equations num(x) - q(x) * den(x) from known points x, q(x)
    // by convention the first coefficient in den is normalised to 1
    let mut lhs = Vec::with_capacity(num_unknown * num_unknown);
    let mut rhs = Vec::with_capacity(num_unknown);
    for (x, q_x) in pts {
        for term in rat.num().terms() {
            lhs.push(eval_pow(term, &x))
        }
        let (normalised, unknown) = rat.den().terms().split_first().unwrap();
        for term in unknown {
            lhs.push(-*q_x * eval_pow(term, &x))
        }
        rhs.push(q_x * eval_pow(normalised, &x));
    }
    let lhs = Matrix::from_vec(num_unknown, lhs);
    let Some(res) = gauss_solve(lhs, rhs) else {
        return None;
    };

    let (num_coeff, den_coeff) = res.split_at(rat.num().len());
    let num = num_coeff
        .iter()
        .copied()
        .zip(rat.num().terms().iter().map(|t| t.powers))
        .map(|(coeff, powers)| SparseMono { powers, coeff })
        .collect();
    let num = SparsePoly::from_raw_terms(num);
    let den = std::iter::once(Z64::one())
        .chain(den_coeff.iter().copied())
        .zip(rat.den().terms().iter().map(|t| t.powers))
        .map(|(coeff, powers)| SparseMono { powers, coeff })
        .collect();
    let den = SparsePoly::from_raw_terms(den);
    Some(Rat::from_num_den_unchecked(num, den))
}

fn eval_pow<T, const P: u64, const N: usize>(
    term: &SparseMono<T, N>,
    x: &[Z64<P>; N],
) -> Z64<P> {
    SparseMono {
        coeff: Z64::<P>::one(),
        powers: term.powers,
    }
    .eval(x)
}

#[cfg(test)]
mod tests {
    use rand::{SeedableRng, Rng};

    use crate::{_test_util::gen_sparse_rat, traits::TryEval};

    use super::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    fn gen_pts<const P: u64, const N: usize>(
        rat: &impl TryEval<[Z64<P>; N], Output = Z64<P>>,
        num: usize,
        mut rng: impl Rng,
    ) -> Vec<([Z64<P>; N], Z64<P>)> {
        let mut pts = Vec::with_capacity(num);
        while pts.len() < num {
            let pt = [(); N].map(|_| rng.gen::<Z64<P>>());
            if let Some(val) = rat.try_eval(&pt) {
                pts.push((pt, val))
            }
        }
        pts
    }

    #[test]
    fn rec_rat2_small() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_POW: u32 = 2;
        const P: u64 = 1152921504606846883;
        const N: usize = 2;
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        for _ in 0..NTESTS {
            let rat = gen_sparse_rat::<P, N>(MAX_POW, MAX_POW, &mut rng);
            eprintln!("trying to reconstruct {rat}");
            let nneeded = rat.num().len() + rat.den().len() - 1;
            let pts = gen_pts(&rat, nneeded, &mut rng);
            let rec = rec_coeff(&rat, &pts).unwrap();
            assert_eq!(rat, rec);
        }
    }

    // TODO: duplication with rec_rat_mod
    #[test]
    fn rec_rat2_simple() {
        log_init();

        const P: u64 = 1152921504606846883;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

        let x = SparseMono::new(Z64::one(), [1, 0]);
        let y = SparseMono::new(Z64::one(), [0, 1]);

        // 1 / y
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![y]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let nneeded = rat.num().len() + rat.den().len() - 1;
        let pts = gen_pts(&rat, nneeded, &mut rng);
        let res = rec_coeff(&rat, &pts).unwrap();
        eprintln!("reconstructed {res}");
        assert_eq!(res, rat);

        // 1 / y^2
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![y * y]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let nneeded = rat.num().len() + rat.den().len() - 1;
        let pts = gen_pts(&rat, nneeded, &mut rng);
        let res = rec_coeff(&rat, &pts).unwrap();
        eprintln!("reconstructed {res}");
        assert_eq!(res, rat);

        // 1 / x
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![x]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let nneeded = rat.num().len() + rat.den().len() - 1;
        let pts = gen_pts(&rat, nneeded, &mut rng);
        let res = rec_coeff(&rat, &pts).unwrap();
        eprintln!("reconstructed {res}");
        assert_eq!(res, rat);
    }


}
