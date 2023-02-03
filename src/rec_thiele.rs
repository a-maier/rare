use std::fmt::{Display, self};

use log::{debug, trace};
use galois_fields::{Z64, TryDiv};
use rand::{Rng, thread_rng};

use crate::{traits::{Zero, One, WithVars, Rec, TryEval}, rand::UniqueRandIter, dense_rat::DenseRat, dense_poly::DensePoly};

/// Univariate rational function reconstruction using Thiele interpolation
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ThieleRec {
    extra_pts: usize
}

impl Default for ThieleRec {
    fn default() -> Self {
        Self { extra_pts: 1 }
    }
}

fn a_next<const P: u64>(
    rat: &ThieleRat<Z64<P>>,
    y_last: Z64<P>,
    y: Z64<P>,
    f_y: Z64<P>
) -> Option<Z64<P>> {
    let mut a = f_y;
    // TODO: check if calculating `a` with less divisions is faster
    for (ai, yi) in &rat.coeffs {
        a = (y - yi).try_div(a - ai)?;
    }
    (y - y_last).try_div(a - rat.a_last)
}

impl ThieleRec {
    pub fn new(extra_pts: usize) -> Self {
        Self { extra_pts }
    }

    pub fn rec_from_seq<I, const P: u64>(
        &self,
        pts: I
    ) -> Option<ThieleRat<Z64<P>>>
    where I: IntoIterator<Item = (Z64<P>, Z64<P>)>
    {

        debug!("1d rational function reconstruction");
        let mut pts = pts.into_iter();
        let Some((y0, a0)) = pts.next() else {
            return None
        };
        trace!("Adding q({y0}) = {a0}");
        let mut rat = ThieleRat::from(a0);
        trace!("q(x1) = {rat}");
        let mut y_last = y0;
        let mut n_zeroes = 0;
        for (y, f_y) in pts {
            if let Some(a) = a_next(&rat, y_last, y, f_y) {
                n_zeroes = 0;
                rat.coeffs.push((rat.a_last, y_last));
                rat.a_last = a;
                trace!("q(x1) = {rat}");
                y_last = y;
            } else {
                n_zeroes += 1;
                trace!("zero ({n_zeroes})");
                if n_zeroes > self.extra_pts {
                    return Some(rat);
                }
            };
        }
        debug!("Reconstruction failed");
        None
    }

    pub fn rec_univariate_with_ran<F, const P: u64>(
        &self,
        poly: F,
        rng: impl Rng,
    ) -> Option<ThieleRat<Z64<P>>>
    where F: FnMut(Z64<P>) -> Option<Z64<P>>
    {
        self.rec1_with_ran(poly, rng)
    }

    pub fn rec1_with_ran<F, const P: u64>(
        &self,
        mut poly: F,
        rng: impl Rng,
    ) -> Option<ThieleRat<Z64<P>>>
    where F: FnMut(Z64<P>) -> Option<Z64<P>>
    {
        self.rec_from_seq(
            UniqueRandIter::new(rng).filter_map(
                |pt| poly(pt).map(|fy| (pt, fy))
            )
        )
    }

    pub fn rec_univariate<F, const P: u64>(
        &self,
        poly: F
    ) -> Option<ThieleRat<Z64<P>>>
    where F: FnMut(Z64<P>) -> Option<Z64<P>>
    {
        self.rec_univariate_with_ran(poly, thread_rng())
    }

}

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ThieleRat<T, U = T> {
    a_last: T,
    coeffs: Vec<(T, U)>,
}

impl<const P: u64> TryEval<Z64<P>> for ThieleRat<Z64<P>> {
    type Output = Z64<P>;

    fn try_eval(&self, x: &Z64<P>) -> Option<Z64<P>> {
        let mut res = self.a_last;
        for (a, y) in self.coeffs.iter().rev() {
            res = a + (x - y).try_div(res)?
        }
        Some(res)
    }
}

impl<T: Zero, U> Zero for ThieleRat<T, U> {
    fn zero() -> Self {
        Self{
            a_last: Zero::zero(),
            coeffs: Vec::new(),
        }
    }

    fn is_zero(&self) -> bool {
        self.coeffs.is_empty() && self.a_last.is_zero()
    }
}

impl<const P: u64> From<&ThieleRat<Z64<P>>> for DenseRat<Z64<P>> {
    fn from(p: &ThieleRat<Z64<P>>) -> Self {
        if p.is_zero() {
            return Zero::zero();
        }
        let mut res_num = DensePoly::from_coeff_unchecked(vec![p.a_last]);
        let mut res_den = One::one();
        for (a, y) in p.coeffs.iter().rev() {
            let xmy = DensePoly::from_coeff_unchecked(vec![-*y, One::one()]);
            (res_num, res_den) = (&res_num * a + xmy * res_den, res_num);
        }
        let norm = res_den.coeffs().last().unwrap().inv();
        res_num *= &norm;
        res_den *= &norm;
        DenseRat::from_num_den_unchecked(res_num, res_den)
    }
}

impl<const P: u64> From<ThieleRat<Z64<P>>> for DenseRat<Z64<P>> {
    fn from(p: ThieleRat<Z64<P>>) -> Self {
        Self::from(&p)
    }
}

impl<'a, 'b, T, U, S: Display> WithVars<'a, &'b [S; 1]> for ThieleRat<T, U>
where
    T: Display + One + Zero + 'a,
    U: Display + One + Zero + 'a,
{
    type Output = FmtThieleRat<'a, 'b, T, U, S>;

    fn with_vars(&'a self, vars: &'b[S; 1]) -> Self::Output {
        FmtThieleRat::new(self, vars)
    }
}

impl<const P: u64> Display for ThieleRat<Z64<P>> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.with_vars(&["x"]).fmt(f)
    }
}

impl<const P: u64> From<Z64<P>> for ThieleRat<Z64<P>> {
    fn from(source: Z64<P>) -> Self {
        Self {
            a_last: source,
            coeffs: Vec::new()
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct FmtThieleRat<'a, 'b, T: Display + One + Zero, U: Display + One + Zero, V: Display> {
    rat: &'a ThieleRat<T, U>,
    var: &'b [V],
}

impl<'a, 'b, T: Display + One + Zero, U: Display + One + Zero, V: Display> FmtThieleRat<'a, 'b, T, U, V> {
    fn new(rat: &'a ThieleRat<T, U>, var: &'b [V]) -> Self {
        Self { rat, var }
    }
}

impl<'a, 'b, V: Display, const P: u64> Display for FmtThieleRat<'a, 'b, Z64<P>, Z64<P>, V> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let var = &self.var[0];
        for (a, y) in &self.rat.coeffs {
            write!(f, "{a} + ({var} - {y}) / (")?;
        }
        write!(f, "{}", self.rat.a_last)?;
        for _ in 0..self.rat.coeffs.len() {
            write!(f, ")")?;
        }
        Ok(())
    }
}

impl<F, const P: u64> Rec<ThieleRec, Z64<P>> for F
where F: FnMut(Z64<P>) -> Option<Z64<P>> {
    type Output = Option<ThieleRat<Z64<P>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: ThieleRec,
        rng: impl ::rand::Rng
    ) -> Self::Output {
        reconstructor.rec_from_seq(
            UniqueRandIter::new(rng).filter_map(
                |pt| (self)(pt).map(|fy| (pt, fy))
            )
        )
    }
}

impl<F, const P: u64> Rec<ThieleRec, [Z64<P>; 1]> for F
where F: FnMut([Z64<P>; 1]) -> Option<Z64<P>> {
    type Output = Option<ThieleRat<Z64<P>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: ThieleRec,
        rng: impl ::rand::Rng
    ) -> Self::Output {
        (|pt| (self)([pt])).rec_with_ran(reconstructor, rng)
    }
}

#[cfg(test)]
mod tests {
    use std::iter::repeat_with;

    use rand::Rng;
    use rand_xoshiro::rand_core::SeedableRng;

    use crate::{dense_poly::DensePoly, traits::TryEval};

    use super::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn rec_rat_small() {
        log_init();

        const NTESTS: u32 = 10;
        const MAX_POW: u32 = 1;
        const P: u64 = 29;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = ThieleRec::new(1);

        for _ in 0..NTESTS {
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let coeff = repeat_with(|| rng.gen::<Z64<P>>()).take(nterms).collect();
            let num = DensePoly::from_coeff(coeff);

            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let mut coeff = Vec::from_iter(
                repeat_with(|| rng.gen::<Z64<P>>()).take(nterms - 1)
            );
            coeff.push(One::one());
            let den = DensePoly::from_coeff(coeff);

            let rat = DenseRat::from_num_den_unchecked(num, den);
            eprintln!("trying to reconstruct {rat}");
            let reconstructed =
                (|x: Z64<P>| rat.try_eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            eprintln!("{reconstructed}");
            let reconstructed: DenseRat<Z64<P>> = reconstructed.into();
            eprintln!("{reconstructed}");
            assert_eq!(rat, reconstructed)
        }
    }

    #[test]
    fn rec_rat_large() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_POW: u32 = 8;
        const P: u64 = 1152921504606846883;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = ThieleRec::new(1);

        for _ in 0..NTESTS {
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let coeff = repeat_with(|| rng.gen::<Z64<P>>()).take(nterms).collect();
            let num = DensePoly::from_coeff(coeff);

            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let mut coeff = Vec::from_iter(
                repeat_with(|| rng.gen::<Z64<P>>()).take(nterms - 1)
            );
            coeff.push(One::one());
            let den = DensePoly::from_coeff(coeff);

            let rat = DenseRat::from_num_den_unchecked(num, den);
            let reconstructed =
                (|x: Z64<P>| rat.try_eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            let reconstructed: DenseRat<Z64<P>> = reconstructed.into();
            assert_eq!(rat, reconstructed)
        }
    }
}
