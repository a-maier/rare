use std::{
    num::NonZeroUsize,
    ops::{Mul, MulAssign, SubAssign},
};

use ffnt::Z64;
use log::{debug, trace};
use num_traits::Inv;

use crate::{
    dense_poly::DensePoly,
    matrix::Matrix,
    rand::pt_iter,
    rat::Rat,
    traits::{One, Rec, Zero},
};

/// Reconstruct a rational function by solving linear systems of equations
///
/// The degrees of both numerator and denominator have to be known.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct LinearRec {
    num_len: usize,
    den_len: NonZeroUsize,
}

impl Default for LinearRec {
    fn default() -> Self {
        Self {
            num_len: Default::default(),
            den_len: One::one(),
        }
    }
}

impl LinearRec {
    pub fn new(num_len: usize, den_len: NonZeroUsize) -> Self {
        Self { num_len, den_len }
    }

    pub fn rec_from_seq<I, const P: u64>(
        &self,
        pts: I,
    ) -> Option<Rat<DensePoly<Z64<P>>>>
    where
        I: IntoIterator<Item = (Z64<P>, Z64<P>)>,
    {
        debug!(
            "1d rational function reconstruction with known degrees {}/{}",
            self.num_len as i64 - 1,
            usize::from(self.den_len) - 1
        );
        if self.num_len == 0 {
            return Some(Zero::zero());
        }
        let ncoeffs = self.num_len + usize::from(self.den_len) - 1;
        let mut eqs = Vec::with_capacity(ncoeffs * ncoeffs);
        let mut rhs = Vec::with_capacity(ncoeffs);
        for (x, q_x) in pts.into_iter().take(ncoeffs) {
            rhs.push(q_x);
            let mut x_to_i = One::one();
            for _ in 0..self.num_len {
                eqs.push(x_to_i);
                x_to_i *= x;
            }
            let mut den_coeff = -q_x * x;
            for _ in 1..usize::from(self.den_len) {
                eqs.push(den_coeff);
                den_coeff *= x;
            }
        }
        if rhs.len() < ncoeffs {
            return None;
        }
        let eqs = Matrix::from_vec(ncoeffs, eqs);
        trace!("Solving system\n{eqs} == {rhs:?}");
        let mut coeffs = gauss_solve(eqs, rhs)?;
        let mut den_coeffs = Vec::with_capacity(usize::from(self.den_len));
        den_coeffs.push(One::one());
        den_coeffs.extend_from_slice(&coeffs[self.num_len..]);
        let den = DensePoly::from_coeff(den_coeffs);
        coeffs.truncate(self.num_len);
        let num = DensePoly::from_coeff(coeffs);
        Some(Rat::from_num_den_unchecked(num, den))
    }
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
        reconstructor.rec_from_seq(
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

// TODO: put rhs inside the matrix?
pub(crate) fn gauss_solve<T>(mut eqs: Matrix<T>, mut rhs: Vec<T>) -> Option<Vec<T>>
where
    T: Copy
        + Zero
        + One
        + Inv<Output = T>
        + MulAssign
        + Mul<Output = T>
        + SubAssign,
{
    debug_assert_eq!(eqs.nrows(), rhs.len());
    // forward substitution
    for nrow in 0..eqs.nrows() {
        if eqs[(nrow, nrow)].is_zero() {
            let pivot = ((nrow + 1)..eqs.nrows())
                .find(|&n| !eqs[(n, nrow)].is_zero())?;
            // TODO: don't we also need to swap the rhs?
            eqs.swap_rows(nrow, pivot);
        }
        debug_assert!(!eqs[(nrow, nrow)].is_zero());
        let pivot = std::mem::replace(&mut eqs[(nrow, nrow)], One::one());
        let inv_pivot = pivot.inv();
        let pivot_row = eqs.row_mut(nrow);
        for e in &mut pivot_row[nrow + 1..] {
            *e *= inv_pivot;
        }
        rhs[nrow] *= inv_pivot;
        for sub_row in (nrow + 1)..eqs.nrows() {
            let fact =
                std::mem::replace(&mut eqs[(sub_row, nrow)], Zero::zero());
            for col in (nrow + 1)..eqs.ncols() {
                let sub = fact * eqs[(nrow, col)];
                eqs[(sub_row, col)] -= sub;
            }
            let sub = fact * rhs[nrow];
            rhs[sub_row] -= sub;
        }
    }

    // backward substitution
    for nrow in (0..eqs.nrows()).rev().skip(1) {
        for ncol in (nrow + 1)..eqs.nrows() {
            let sub = eqs[(nrow, ncol)] * rhs[ncol];
            rhs[nrow] -= sub;
        }
    }
    Some(rhs)
}

#[cfg(test)]
mod tests {
    use super::*;

    use ffnt::Z64;
    use rand_xoshiro::rand_core::SeedableRng;

    use crate::{traits::{One, TryEval}, _test_util::{gen_dense_rat1, sample_eq}};

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn solve_linear() {
        const P: u64 = 7;
        let one: Z64<P> = One::one();
        let sys = Matrix::from_vec(2, vec![one, -one, one, one]);
        let rhs = vec![one, one];
        let x = gauss_solve(sys, rhs).unwrap();
        assert_eq!(x, [one, Zero::zero()]);
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
}
