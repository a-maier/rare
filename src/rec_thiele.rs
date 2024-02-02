use std::{fmt::{self, Display}, ops::ControlFlow};

use ffnt::{TryDiv, Z64};
use log::{debug, trace};
use rand::{thread_rng, Rng};

use crate::{
    dense_poly::DensePoly,
    rand::pt_iter,
    rat::Rat,
    traits::{One, Rec, TryEval, WithVars, Zero}, sparse_poly::SparsePoly,
};

/// Univariate rational function reconstruction using Thiele interpolation
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ThieleRec<const P: u64> {
    extra_pts: usize,
    rat: ThieleRat<Z64<P>>,
    y_last: Z64<P>,
    rec_started: bool,
    n_zeroes: usize,
}

impl<const P: u64> Default for ThieleRec<P> {
    fn default() -> Self {
        Self::new(1)
    }
}


impl<const P: u64> ThieleRec<P> {
    pub fn new(extra_pts: usize) -> Self {
        Self {
            extra_pts,
            rat: Default::default(),
            y_last: Z64::zero(),
            rec_started: false,
            n_zeroes: 0,
        }
    }

    pub fn add_pt(
        &mut self,
        y: Z64<P>,
        q_y: Z64<P>
    ) -> ControlFlow<()> {
        trace!("Adding q({y}) = {q_y}");
        if !self.rec_started {
            self.rec_started = true;
            self.rat = ThieleRat::from(q_y);
            self.y_last = y;
            return ControlFlow::Continue(());
        }
        if let Some(a) = self.a_next(y, q_y) {
            self.n_zeroes = 0;
            self.rat.coeffs.push((self.rat.a_last, self.y_last));
            self.rat.a_last = a;
            trace!("q(x1) = {}", self.rat);
            self.y_last = y;
        } else {
            self.n_zeroes += 1;
            trace!("zero ({})", self.n_zeroes);
            if self.n_zeroes > self.extra_pts {
                debug!("Reconstructed {}", self.rat);
                return ControlFlow::Break(())
            }
        };
        ControlFlow::Continue(())
    }

    fn a_next(
        &self,
        y: Z64<P>,
        f_y: Z64<P>,
    ) -> Option<Z64<P>> {
        let mut a = f_y;
        // TODO: check if calculating `a` with less divisions is faster
        for (ai, yi) in &self.rat.coeffs {
            a = (y - yi).try_div(a - ai)?;
        }
        (y - self.y_last).try_div(a - self.rat.a_last)
    }

    pub fn into_rat(self) -> ThieleRat<Z64<P>> {
        self.rat
    }

    pub fn rec_from_seq<I>(
        mut self,
        pts: I,
    ) -> Option<ThieleRat<Z64<P>>>
    where
        I: IntoIterator<Item = (Z64<P>, Z64<P>)>,
    {
        debug!("1d rational function reconstruction");
        for (y, q_y) in pts {
            if self.add_pt(y, q_y) == ControlFlow::Break(()) {
                return Some(self.into_rat())
            }
        }
        debug!("Reconstruction failed");
        None
    }

    pub fn rec_univariate_with_ran<F>(
        self,
        poly: F,
        rng: impl Rng,
    ) -> Option<ThieleRat<Z64<P>>>
    where
        F: FnMut(Z64<P>) -> Option<Z64<P>>,
    {
        self.rec1_with_ran(poly, rng)
    }

    pub fn rec1_with_ran<F>(
        self,
        mut poly: F,
        rng: impl Rng,
    ) -> Option<ThieleRat<Z64<P>>>
    where
        F: FnMut(Z64<P>) -> Option<Z64<P>>,
    {
        self.rec_from_seq(
            pt_iter(rng).filter_map(|pt| poly(pt).map(|fy| (pt, fy))),
        )
    }

    pub fn rec_univariate<F>(
        self,
        poly: F,
    ) -> Option<ThieleRat<Z64<P>>>
    where
        F: FnMut(Z64<P>) -> Option<Z64<P>>,
    {
        self.rec_univariate_with_ran(poly, thread_rng())
    }

    pub fn rec_started(&self) -> bool {
        self.rec_started
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
        Self {
            a_last: Zero::zero(),
            coeffs: Vec::new(),
        }
    }

    fn is_zero(&self) -> bool {
        self.coeffs.is_empty() && self.a_last.is_zero()
    }
}

impl<const P: u64> From<&ThieleRat<Z64<P>>> for Rat<DensePoly<Z64<P>>> {
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
        let norm = res_num
            .coeffs()
            .iter()
            .rfind(|c| !c.is_zero())
            .unwrap()
            .inv();
        res_num *= &norm;
        res_den *= &norm;
        Rat::from_num_den_unchecked(res_num, res_den)
    }
}

impl<const P: u64> From<ThieleRat<Z64<P>>> for Rat<DensePoly<Z64<P>>> {
    fn from(p: ThieleRat<Z64<P>>) -> Self {
        Self::from(&p)
    }
}

impl<const P: u64> From<&ThieleRat<Z64<P>>> for Rat<SparsePoly<Z64<P>, 1>> {
    fn from(p: &ThieleRat<Z64<P>>) -> Self {
        let (num, den) = Rat::<DensePoly<Z64<P>>>::from(p).into_num_den();
        Rat::from_num_den_unchecked(num.into(), den.into())
    }
}

impl<const P: u64> From<ThieleRat<Z64<P>>> for Rat<SparsePoly<Z64<P>, 1>> {
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

    fn with_vars(&'a self, vars: &'b [S; 1]) -> Self::Output {
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
            coeffs: Vec::new(),
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct FmtThieleRat<
    'a,
    'b,
    T: Display + One + Zero,
    U: Display + One + Zero,
    V: Display,
> {
    rat: &'a ThieleRat<T, U>,
    var: &'b [V],
}

impl<'a, 'b, T: Display + One + Zero, U: Display + One + Zero, V: Display>
    FmtThieleRat<'a, 'b, T, U, V>
{
    fn new(rat: &'a ThieleRat<T, U>, var: &'b [V]) -> Self {
        Self { rat, var }
    }
}

impl<'a, 'b, V: Display, const P: u64> Display
    for FmtThieleRat<'a, 'b, Z64<P>, Z64<P>, V>
{
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

impl<F, const P: u64> Rec<ThieleRec<P>, Z64<P>> for F
where
    F: FnMut(Z64<P>) -> Option<Z64<P>>,
{
    type Output = Option<ThieleRat<Z64<P>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: ThieleRec<P>,
        rng: impl ::rand::Rng,
    ) -> Self::Output {
        reconstructor.rec_from_seq(
            pt_iter(rng).filter_map(|pt| (self)(pt).map(|fy| (pt, fy))),
        )
    }
}

impl<F, const P: u64> Rec<ThieleRec<P>, [Z64<P>; 1]> for F
where
    F: FnMut([Z64<P>; 1]) -> Option<Z64<P>>,
{
    type Output = Option<ThieleRat<Z64<P>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: ThieleRec<P>,
        rng: impl ::rand::Rng,
    ) -> Self::Output {
        (|pt| (self)([pt])).rec_with_ran(reconstructor, rng)
    }
}

#[cfg(test)]
mod tests {
    use rand_xoshiro::rand_core::SeedableRng;

    use crate::{dense_poly::DensePoly, traits::TryEval, _test_util::{gen_dense_rat1, sample_eq}};

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

        for _ in 0..NTESTS {
            let rec = ThieleRec::new(1);
            let rat = gen_dense_rat1(&[MAX_POW], &mut rng);
            eprintln!("trying to reconstruct {rat}");
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
            let rec = ThieleRec::new(1);
            let rat = gen_dense_rat1(&[MAX_POW], &mut rng);
            let reconstructed = (|x: Z64<P>| rat.try_eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            let reconstructed: Rat<DensePoly<Z64<P>>> = reconstructed.into();
            assert!(sample_eq(&rat, &reconstructed, &mut rng))
        }
    }
}
