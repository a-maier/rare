use ffnt::Z64;

use crate::{
    matrix::Matrix,
    rat::Rat,
    sparse_poly::{SparseMono, SparsePoly},
    traits::{Eval, One, Zero},
};

/// Reconstruction using a linear system of equations built from the given points
pub(crate) trait RecLinear<Pts> {
    type Output;

    /// Attempt the reconstruction
    fn rec_linear(&self, pts: &[Pts]) -> Option<Self::Output>;
}

impl<T, const P: u64, const N: usize> RecLinear<([Z64<P>; N], Z64<P>)> for Rat<SparsePoly<T, N>> {
    type Output = Rat<SparsePoly<Z64<P>, N>>;

    fn rec_linear(
        &self,
        pts: &[([Z64<P>; N], Z64<P>)]
    ) -> Option<Rat<SparsePoly<Z64<P>, N>>> {
        assert!(!self.den().is_empty());
        let num_coeff = self.num().len() + self.den().len();
        // Need at least one equation for each *unknown* coefficient in rat
        // one coefficient is set to one by convention
        if pts.len() < num_coeff - 1 {
            return None;
        }

        // construct equations num(x) - q(x) * den(x) from known points x, q(x)
        // by convention the coefficient of the highest-order monomial in
        // the numerator is set to one
        let mut eqs = Vec::with_capacity(num_coeff * pts.len());
        for (x, q_x) in pts {
            for term in self.den().terms() {
                eqs.push(-*q_x * eval_pow(term, x));
            }
            for term in self.num().terms() {
                eqs.push(eval_pow(term, x));
            }
        }
        let mut eqs = Matrix::from_vec(pts.len(), eqs);
        eqs.row_reduce();
        let mut den_coeffs = Vec::with_capacity(self.den().len());
        let mut num_coeffs = Vec::with_capacity(self.num().len());
        for (i, row) in eqs.rows().enumerate() {
            if !row[i].is_one() || row[(i+1)..(eqs.ncols() - 1)].iter().any(|z| !z.is_zero()){
                return None;
            }
            let coeff = -*row.last().unwrap();
            if i < usize::from(self.den().len()) {
                den_coeffs.push(coeff);
            } else {
                num_coeffs.push(coeff);
            }
        }
        num_coeffs.push(Z64::one());
        let num = num_coeffs
            .iter()
            .copied()
            .zip(self.num().terms().iter().map(|t| t.powers))
            .map(|(coeff, powers)| SparseMono { powers, coeff })
            .collect();
        let num = SparsePoly::from_raw_terms(num);
        let den = den_coeffs
            .iter()
            .copied()
            .zip(self.den().terms().iter().map(|t| t.powers))
            .map(|(coeff, powers)| SparseMono { powers, coeff })
            .collect();
        let den = SparsePoly::from_raw_terms(den);
        Some(Rat::from_num_den_unchecked(num, den))
    }
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

    use crate::{_test_util::{gen_sparse_rat, sample_eq}, traits::TryEval};

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
            let rec = rat.rec_linear(&pts).unwrap();
            eprintln!("reconstructed {rec}");
            assert!(sample_eq(&rat, &rec, &mut rng));
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
        let res = rat.rec_linear(&pts).unwrap();
        eprintln!("reconstructed {res}");
        assert!(sample_eq(&rat, &res, &mut rng));

        // 1 / y^2
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![y * y]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let nneeded = rat.num().len() + rat.den().len() - 1;
        let pts = gen_pts(&rat, nneeded, &mut rng);
        let res = rat.rec_linear(&pts).unwrap();
        eprintln!("reconstructed {res}");
        assert!(sample_eq(&rat, &res, &mut rng));

        // 1 / x
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![x]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let nneeded = rat.num().len() + rat.den().len() - 1;
        let pts = gen_pts(&rat, nneeded, &mut rng);
        let res = rat.rec_linear(&pts).unwrap();
        eprintln!("reconstructed {res}");
        assert!(sample_eq(&rat, &res, &mut rng));
    }


}
