use std::num::NonZeroUsize;

use ffnt::Z64;
use log::debug;

use crate::{
    matrix::Matrix,
    rand::pt_iter,
    rat::Rat,
    sparse_poly::{SparseMono, SparsePoly},
    traits::{Eval, One, Zero, Rec}, dense_poly::DensePoly,
};

/// Reconstruction using a linear system of equations built from the given points
pub(crate) trait RecLinear<Pts> {
    type Output;

    /// Attempt the reconstruction
    fn rec_linear(&self, pts: Pts) -> Option<Self::Output>;
}

/// Reconstruct a rational function by solving linear systems of equations
///
/// The degrees of both numerator and denominator have to be known.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct LinearRec {
    num_len: usize,
    den_len: NonZeroUsize,
}

impl LinearRec {
    pub fn new(num_len: usize, den_len: NonZeroUsize) -> Self {
        Self { num_len, den_len }
    }
}

impl<Pts, const P: u64> RecLinear<Pts> for LinearRec
where
    Pts: IntoIterator<Item = (Z64<P>, Z64<P>)>,
{
    type Output = Rat<DensePoly<Z64<P>>>;

    fn rec_linear(
        &self,
        pts: Pts
    ) -> Option<Rat<DensePoly<Z64<P>>>> {
        let pts = pts.into_iter();
        debug!(
            "1d rational function reconstruction with known degrees {}/{}",
            self.num_len as i64 - 1,
            usize::from(self.den_len) - 1
        );
        if self.num_len == 0 {
            return Some(Zero::zero());
        }
        let num_coeffs = self.num_len + usize::from(self.den_len);
        let neqs = std::cmp::max(pts.size_hint().0, num_coeffs - 1);
        debug!("Reconstructing with {neqs} equations");
        let mut eqs = Vec::with_capacity(neqs * num_coeffs);
        for (x, q_x) in pts.take(neqs) {
            let mut den_coeff = -q_x;
            for _ in 0..usize::from(self.den_len) {
                eqs.push(den_coeff);
                den_coeff *= x;
            }
            let mut x_to_i = One::one();
            for _ in 0..self.num_len {
                eqs.push(x_to_i);
                x_to_i *= x;
            }
        }
        let (den_coeffs, num_coeffs) = solve_eqs(
            neqs,
            eqs,
            usize::from(self.den_len)
        )?;
        let num = DensePoly::from_coeff(num_coeffs);
        let den = DensePoly::from_coeff(den_coeffs);
        Some(Rat::from_num_den_unchecked(num, den))
    }
}

impl<Pts, T, const P: u64, const N: usize> RecLinear<Pts> for Rat<SparsePoly<T, N>>
where
    Pts: IntoIterator<Item = ([Z64<P>; N], Z64<P>)>,
{
    type Output = Rat<SparsePoly<Z64<P>, N>>;

    fn rec_linear(
        &self,
        pts: Pts
    ) -> Option<Rat<SparsePoly<Z64<P>, N>>> {
        let pts = pts.into_iter();
        assert!(!self.den().is_empty());
        if self.num().len() == 0 {
            return Some(Zero::zero());
        }
        let num_coeffs = self.num().len() + self.den().len();
        let neqs = std::cmp::max(pts.size_hint().0, num_coeffs - 1);
        debug!("Reconstructing with {neqs} equations");

        // construct equations num(x) - q(x) * den(x) from known points x, q(x)
        // by convention the coefficient of the highest-order monomial in
        // the numerator is set to one
        let mut eqs = Vec::with_capacity(neqs * num_coeffs);
        for (x, q_x) in pts.take(neqs) {
            for term in self.den().terms() {
                eqs.push(-q_x * eval_pow(term, &x));
            }
            for term in self.num().terms() {
                eqs.push(eval_pow(term, &x));
            }
        }
        let (den_coeffs, num_coeffs) = solve_eqs(neqs, eqs, self.den().len())?;
        let num = num_coeffs
            .into_iter()
            .zip(self.num().terms().iter().map(|t| t.powers))
            .map(|(coeff, powers)| SparseMono { powers, coeff })
            .collect();
        let num = SparsePoly::from_raw_terms(num);
        let den = den_coeffs
            .into_iter()
            .zip(self.den().terms().iter().map(|t| t.powers))
            .map(|(coeff, powers)| SparseMono { powers, coeff })
            .collect();
        let den = SparsePoly::from_raw_terms(den);
        Some(Rat::from_num_den_unchecked(num, den))
    }
}

// Solve `neqs` linear equations Σ c_i x_i = 0
//
// Returns `None` iff there is no non-trivial solution. Otherwise
// returns one non-trivial solution with the first `nden` x_i in one
// vector and the remainder in a second vector. The solution is chosen
// such that as many x_i with the highest i as possible are set to zero
// and one set to unity. TODO: implement that
fn solve_eqs<const P: u64>(
    neqs: usize,
    eqs: Vec<Z64<P>>,
    nden: usize,
) -> Option<(Vec<Z64<P>>, Vec<Z64<P>>)> {
    debug_assert!(nden <= neqs);
    if eqs.len() < neqs * (neqs + 1) {
        return None;
    }
    let mut eqs = Matrix::from_vec(neqs, eqs);
    eqs.row_reduce();
    eqs.trim_end();
    debug!("After row reduction:\n{eqs}");
    let nnum =  eqs.ncols() - nden;
    let mut den_coeffs = Vec::with_capacity(nden);
    let mut num_coeffs = Vec::with_capacity(nnum);
    for (i, row) in eqs.rows().enumerate() {
        // TODO: fix
        if !row[i].is_one() || row[(i+1)..(eqs.ncols() - 1)].iter().any(|z| !z.is_zero()){
            return None;
        }
        let coeff = -*row.last().unwrap();
        if i < nden {
            den_coeffs.push(coeff);
        } else {
            num_coeffs.push(coeff);
        }
    }
    num_coeffs.push(Z64::one());
    Some((den_coeffs, num_coeffs))
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

impl<F, const P: u64> Rec<LinearRec, Z64<P>> for F
where
    F: FnMut(Z64<P>) -> Option<Z64<P>>,
{
    type Output = Option<Rat<DensePoly<Z64<P>>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: LinearRec,
        rng: impl ::rand::Rng,
    ) -> Self::Output {
        reconstructor.rec_linear(
            pt_iter(rng).filter_map(|pt| (self)(pt).map(|fy| (pt, fy))),
        )
    }
}

impl<F, const P: u64> Rec<LinearRec, [Z64<P>; 1]> for F
where
    F: FnMut([Z64<P>; 1]) -> Option<Z64<P>>,
{
    type Output = Option<Rat<DensePoly<Z64<P>>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: LinearRec,
        rng: impl ::rand::Rng,
    ) -> Self::Output {
        (|pt| (self)([pt])).rec_with_ran(reconstructor, rng)
    }
}

#[cfg(test)]
mod tests {
    use rand::{SeedableRng, Rng};

    use crate::{_test_util::{gen_sparse_rat, sample_eq, gen_dense_rat1}, traits::TryEval};

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
            let rec = rat.rec_linear(pts).unwrap();
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
        let res = rat.rec_linear(pts).unwrap();
        eprintln!("reconstructed {res}");
        assert!(sample_eq(&rat, &res, &mut rng));

        // 1 / y^2
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![y * y]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let nneeded = rat.num().len() + rat.den().len() - 1;
        let pts = gen_pts(&rat, nneeded, &mut rng);
        let res = rat.rec_linear(pts).unwrap();
        eprintln!("reconstructed {res}");
        assert!(sample_eq(&rat, &res, &mut rng));

        // 1 / x
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![x]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let nneeded = rat.num().len() + rat.den().len() - 1;
        let pts = gen_pts(&rat, nneeded, &mut rng);
        let res = rat.rec_linear(pts).unwrap();
        eprintln!("reconstructed {res}");
        assert!(sample_eq(&rat, &res, &mut rng));
    }

    #[test]
    fn rec_rat_small() {
        log_init();

        const NTESTS: u32 = 10;
        const MAX_POW: u32 = 1;
        const P: u64 = 29;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

        for _ in 0..NTESTS {
            let rat = gen_dense_rat1(&[MAX_POW], &mut rng);
            eprintln!("trying to reconstruct {rat}");
            let rec = LinearRec::new(
                rat.num().len(),
                rat.den().len().try_into().unwrap(),
            );
            let reconstructed = (|x: Z64<P>| rat.try_eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            eprintln!("{reconstructed}");
            let reconstructed: Rat<DensePoly<Z64<P>>> = reconstructed.into();
            eprintln!("{reconstructed}");
            assert!(sample_eq(&rat, &reconstructed, &mut rng))
        }
    }

    #[test]
    fn rec_rat_large() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_POW: u32 = 8;
        const P: u64 = 1152921504606846883;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

        for _ in 0..NTESTS {
            let rat = gen_dense_rat1(&[MAX_POW], &mut rng);
            let rec = LinearRec::new(
                rat.num().len(),
                rat.den().len().try_into().unwrap(),
            );
            let reconstructed = (|x: Z64<P>| rat.try_eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            let reconstructed: Rat<DensePoly<Z64<P>>> = reconstructed.into();
            assert!(sample_eq(&rat, &reconstructed, &mut rng))
        }
    }

    #[test]
    fn rec_unknown_degree() {
        log_init();
        const P: u64 = 1152921504606846883;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

        let rec = LinearRec::new(3, 3.try_into().unwrap());

        let mut rat = |x: Z64<P>| Some(x - Z64::one());
        let reconstructed = rat.rec_with_ran(rec, &mut rng).unwrap();

    }

}
