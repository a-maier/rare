use std::{
    fmt::{self, Display},
    ops::ControlFlow,
};

use galois_fields::Z64;
use log::{debug, trace};
use paste::paste;
use rand::{thread_rng, Rng};

use crate::{
    dense_poly::DensePoly,
    traits::{Eval, TryEval},
};
use crate::{
    rand::pt_iter,
    traits::{One, Rec, WithVars, Zero},
    util::{ALL_VARS, slice_start},
};

/// Univariate polynomial reconstruction using Newton interpolation
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct NewtonPolyRec<const P: u64> {
    poly: NewtonPoly<Z64<P>>,
    y_last: Option<Z64<P>>,
    extra_pts: usize,
}

impl<const P: u64> Default for NewtonPolyRec<P> {
    fn default() -> Self {
        Self {
            poly: Zero::zero(),
            y_last: None,
            extra_pts: 1,
        }
    }
}

impl<const P: u64> NewtonPolyRec<P> {
    pub fn new(extra_pts: usize) -> Self {
        Self {
            extra_pts,
            ..Default::default()
        }
    }

    pub fn add_pt(
        &mut self,
        y: &[Z64<P>; 1],
        f_y: Z64<P>,
    ) -> ControlFlow<(), usize> {
        trace!("adding p({}) = {f_y}", y[0]);
        if let Some(y_last) = self.y_last {
            let Some(a) = self.next_a(y[0], f_y) else {
                return ControlFlow::Continue(0)
            };
            self.poly.coeffs.push((self.poly.a_last, y_last));
            self.poly.a_last = a;
            trace!("p(x0) = {}", self.poly);
        } else {
            self.poly.a_last = f_y;
        }
        self.y_last = Some(y[0]);
        self.next_step()
    }

    fn next_a(&self, y: Z64<P>, f_y: Z64<P>) -> Option<Z64<P>> {
        let prefact = self
            .poly
            .coeffs
            .iter()
            .map(|(_ai, yi)| y - yi)
            .fold(y - self.y_last.unwrap(), |acc, x| acc * x);
        let prefact = prefact.try_inv()?;
        let poly_at_y = self.poly.eval(&y);
        Some(prefact * (f_y - poly_at_y))
    }

    fn next_step(&self) -> ControlFlow<(), usize> {
        if self.trailing_zeros() > self.extra_pts {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(0)
        }
    }

    fn trailing_zeros(&self) -> usize {
        if !self.poly.a_last.is_zero() {
            return 0;
        }
        1 + self
            .poly
            .coeffs
            .iter()
            .rev()
            .take_while(|(a, _)| a.is_zero())
            .count()
    }

    pub fn into_poly(mut self) -> NewtonPoly<Z64<P>> {
        let trailing_zeros = self.trailing_zeros();
        if trailing_zeros > 0 {
            self.poly
                .coeffs
                .truncate(1 + self.poly.coeffs.len() - trailing_zeros);
            if let Some((a, _y)) = self.poly.coeffs.pop() {
                self.poly.a_last = a;
            }
        }
        debug!("Reconstructed {}", self.poly);
        self.poly
    }
}

pub type NewtonPolyRec1<const P: u64> = NewtonPolyRec<P>;

macro_rules! impl_newton_rec_recursive {
    ( $($x:literal, $y:literal), * ) => {
        $(
            paste! {
                /// Multivariate polynomial reconstruction using Newton interpolation
                #[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
                pub struct [<NewtonPolyRec $x>]<const P: u64> {
                    poly: [<NewtonPoly $x>]<Z64<P>>,
                    next_a: [<NewtonPolyRec $y>]<P>,
                    last_y: Option<Z64<P>>,
                    cur_y: Option<Z64<P>>,
                    inv_y_prod: Z64<P>,
                    extra_pts: usize
                }

                impl<const P: u64> Default for [<NewtonPolyRec $x>]<P> {
                    fn default() -> Self {
                        Self {
                            poly: Zero::zero(),
                            cur_y: None,
                            last_y: None,
                            next_a: Default::default(),
                            inv_y_prod: One::one(),
                            extra_pts: 1
                        }
                    }
                }

                impl<const P: u64> [<NewtonPolyRec $x>]<P> {
                    pub fn new(extra_pts: usize) -> Self {
                        Self {
                            extra_pts,
                            next_a: [<NewtonPolyRec $y>]::new(extra_pts),
                            inv_y_prod: One::one(),
                            ..Default::default()
                        }
                    }

                    pub fn add_pt(
                        &mut self,
                        y: &[Z64<P>; $x],
                        f_y: Z64<P>
                    ) -> ControlFlow<(), usize> {
                        trace!("adding p({y:?}) = {f_y}");
                        if self.cur_y.is_none() {
                            // new value of first (outermost) variable
                            // cache new inverse product of y differences
                            let y_prod = self.poly.coeffs.iter()
                                .map(|(_ai, yi)| y[0] - yi)
                                .reduce(|acc, x| acc * x)
                                .unwrap_or(One::one());
                            let Some(inv_y_prod) = y_prod.try_inv() else {
                                return ControlFlow::Continue(0)
                            };
                            self.inv_y_prod = inv_y_prod;
                            self.cur_y = Some(y[0]);
                        }
                        if self.cur_y == Some(y[0]) {
                            // work on reconstructing the next coefficient
                            let next_a_val = self.inv_y_prod * (f_y - self.poly.eval(y));
                            let y_rest: &[Z64<P>; $y] = &y[1..].try_into().unwrap();
                            match self.next_a.add_pt(y_rest, next_a_val) {
                                ControlFlow::Continue(n) => ControlFlow::Continue(n + 1),
                                ControlFlow::Break(_) => {
                                    // add reconstructed coefficient to polynomial
                                    let a = std::mem::replace(
                                        &mut self.next_a,
                                        [<NewtonPolyRec $y>]::new(self.extra_pts)
                                    ).into_poly();
                                    let a_prev = std::mem::replace(&mut self.poly.a_last, a);
                                    if let Some(last_y) = self.last_y {
                                        self.poly.coeffs.push((a_prev, last_y));
                                    }
                                    trace!("p = {}", self.poly);
                                    self.last_y = Some(y[0]);
                                    self.cur_y = None;
                                    self.next_step()
                                },
                            }
                        } else {
                            // TODO:
                            // Got a new value of the outermost variable before
                            // finishing previous reconstructions
                            // Should probably signal an error here
                            debug!("Bad y value: {}", y[0]);
                            ControlFlow::Continue($x)
                        }
                    }

                    fn next_step(&self) -> ControlFlow<(), usize> {
                        if self.trailing_zeros() > self.extra_pts {
                            ControlFlow::Break(())
                        } else {
                            ControlFlow::Continue(0)
                        }
                    }

                    fn trailing_zeros(&self) -> usize {
                        if !self.poly.a_last.is_zero() {
                            return 0;
                        }
                        1 + self.poly.coeffs.iter().rev()
                            .take_while(|(a, _)| a.is_zero())
                            .count()
                    }

                    pub fn into_poly(mut self) -> [<NewtonPoly $x>]<Z64<P>> {
                        let trailing_zeros = self.trailing_zeros();
                        if trailing_zeros > 0 {
                            self.poly.coeffs.truncate(1 + self.poly.coeffs.len() - trailing_zeros);
                            if let Some((a, _y)) = self.poly.coeffs.pop() {
                                self.poly.a_last = a;
                            }
                        }
                        debug!("Reconstructed {}", self.poly);
                        self.poly
                    }
                }
            }

        )*
    };
}

/// Univariate polynomial reconstruction using Newton interpolation
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct NewtonRec {
    extra_pts: usize,
}

impl Default for NewtonRec {
    fn default() -> Self {
        Self { extra_pts: 1 }
    }
}

impl NewtonRec {
    pub fn new(extra_pts: usize) -> Self {
        Self { extra_pts }
    }

    pub fn rec_from_seq<I, const P: u64>(
        &self,
        pts: I,
    ) -> Option<NewtonPoly<Z64<P>>>
    where
        I: IntoIterator<Item = (Z64<P>, Z64<P>)>,
    {
        debug!("1d polynomial reconstruction");
        let mut rec = NewtonPolyRec::new(self.extra_pts);
        for (y, f_y) in pts.into_iter() {
            if rec.add_pt(&[y], f_y).is_break() {
                return Some(rec.into_poly());
            }
        }
        debug!("Reconstruction failed");
        None
    }

    pub fn rec_univariate_with_ran<F, const P: u64>(
        &self,
        poly: F,
        rng: impl Rng,
    ) -> Option<NewtonPoly<Z64<P>>>
    where
        F: FnMut(Z64<P>) -> Z64<P>,
    {
        self.rec1_with_ran(poly, rng)
    }

    pub fn rec1_with_ran<F, const P: u64>(
        &self,
        mut poly: F,
        rng: impl Rng,
    ) -> Option<NewtonPoly<Z64<P>>>
    where
        F: FnMut(Z64<P>) -> Z64<P>,
    {
        self.rec_from_seq(pt_iter(rng).map(|pt| (pt, poly(pt))))
    }

    pub fn rec_univariate<F, const P: u64>(
        &self,
        poly: F,
    ) -> Option<NewtonPoly<Z64<P>>>
    where
        F: FnMut(Z64<P>) -> Z64<P>,
    {
        self.rec_univariate_with_ran(poly, thread_rng())
    }
}

impl<F, const P: u64> Rec<NewtonRec, (Z64<P>,)> for F
where
    F: FnMut(Z64<P>) -> Z64<P>,
{
    type Output = Option<NewtonPoly1<Z64<P>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: NewtonRec,
        rng: impl Rng,
    ) -> Self::Output {
        reconstructor.rec_from_seq(pt_iter(rng).map(|pt| (pt, (self)(pt))))
    }
}

impl<F, const P: u64> Rec<NewtonRec, [Z64<P>; 1]> for F
where
    F: FnMut([Z64<P>; 1]) -> Z64<P>,
{
    type Output = Option<NewtonPoly1<Z64<P>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: NewtonRec,
        rng: impl Rng,
    ) -> Self::Output {
        reconstructor.rec_from_seq(pt_iter(rng).map(|pt| (pt, (self)([pt]))))
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct NewtonPoly<T, U = T> {
    a_last: T,
    coeffs: Vec<(T, U)>,
}

impl<T: Zero, U> Zero for NewtonPoly<T, U> {
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

impl<T: One, U> One for NewtonPoly<T, U> {
    fn one() -> Self {
        Self {
            a_last: One::one(),
            coeffs: Vec::new(),
        }
    }

    fn is_one(&self) -> bool {
        self.coeffs.is_empty() && self.a_last.is_one()
    }
}

impl<const P: u64> From<Z64<P>> for NewtonPoly<Z64<P>> {
    fn from(z: Z64<P>) -> Self {
        Self {
            a_last: z,
            coeffs: Vec::new(),
        }
    }
}

impl<const P: u64> From<&NewtonPoly<Z64<P>>> for DensePoly<Z64<P>> {
    fn from(p: &NewtonPoly<Z64<P>>) -> Self {
        if p.is_zero() {
            return Zero::zero();
        }
        let mut res = DensePoly::from_coeff_unchecked(vec![p.a_last]);
        for (a, y) in p.coeffs.iter().rev() {
            let xmy = Self::from_coeff_unchecked(vec![-*y, One::one()]);
            res = xmy * res + a;
        }
        res
    }
}

impl<const P: u64> TryEval<Z64<P>> for NewtonPoly<Z64<P>> {
    type Output = Z64<P>;

    fn try_eval(&self, x: &Z64<P>) -> Option<Z64<P>> {
        let mut res = self.a_last;
        for (a, y) in self.coeffs.iter().rev() {
            res = a + (x - y) * res;
        }
        Some(res)
    }
}
impl<const P: u64> Eval<Z64<P>> for NewtonPoly<Z64<P>> {}

impl<const P: u64> TryEval<[Z64<P>; 1]> for NewtonPoly<Z64<P>> {
    type Output = Z64<P>;

    fn try_eval(&self, x: &[Z64<P>; 1]) -> Option<Z64<P>> {
        self.try_eval(&x[0])
    }
}
impl<const P: u64> Eval<[Z64<P>; 1]> for NewtonPoly<Z64<P>> {}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct FmtNewtonPoly<
    'a,
    'b,
    T: Display + One + Zero,
    U: Display + One + Zero,
    V: Display,
> {
    poly: &'a NewtonPoly<T, U>,
    var: &'b [V],
}

impl<'a, 'b, T: Display + One + Zero, U: Display + One + Zero, V: Display>
    FmtNewtonPoly<'a, 'b, T, U, V>
{
    fn new(poly: &'a NewtonPoly<T, U>, var: &'b [V]) -> Self {
        Self { poly, var }
    }
}

impl<'a, 'b, V: Display, const P: u64> Display
    for FmtNewtonPoly<'a, 'b, Z64<P>, Z64<P>, V>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let var = &self.var[0];
        for (a, y) in &self.poly.coeffs {
            write!(f, "{a} + ({var} - {y}) * (")?;
        }
        write!(f, "{}", self.poly.a_last)?;
        for _ in 0..self.poly.coeffs.len() {
            write!(f, ")")?;
        }
        Ok(())
    }
}

pub type NewtonPoly1<T> = NewtonPoly<T>;

impl<'a, 'b, T, U, S: Display> WithVars<'a, &'b [S; 1]> for NewtonPoly<T, U>
where
    T: Display + One + Zero + 'a,
    U: Display + One + Zero + 'a,
{
    type Output = FmtNewtonPoly<'a, 'b, T, U, S>;

    fn with_vars(&'a self, vars: &'b [S; 1]) -> Self::Output {
        FmtNewtonPoly::new(self, vars)
    }
}

macro_rules! impl_newton_poly {
    ( $($x:literal), *) => {
        $(
            paste! {
                use crate::dense_poly::[<DensePoly $x>];

                impl<const P: u64> Display for [<NewtonPoly $x>]<Z64<P>> {
                    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                        let vars: &[_; $x] = slice_start(&ALL_VARS);
                        self.with_vars(vars).fmt(f)
                    }
                }

                impl<'a, 'b, V: Display, const P: u64> Display for FmtNewtonPoly<'a, 'b, [<NewtonPoly $x>]<Z64<P>>, Z64<P>, V> {
                    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                        let var = &self.var[0];
                        for (a, y) in &self.poly.coeffs {
                            let coeff = FmtNewtonPoly::new(a, &self.var[1..]);
                            write!(f, "{coeff} + ({var} - {y}) * (")?;
                        }
                        let coeff = FmtNewtonPoly::new(&self.poly.a_last, &self.var[1..]);
                        write!(f, "{coeff}")?;
                        for _ in 0..self.poly.coeffs.len() {
                            write!(f, ")")?;
                        }
                        Ok(())
                    }
                }

                impl<const P: u64> From<[<NewtonPoly $x>]<Z64<P>>> for [<DensePoly $x>]<Z64<P>> {
                    fn from(p: [<NewtonPoly $x>]<Z64<P>>) -> Self {
                        Self::from(&p)
                    }
                }
            }
        )*
    };
}

macro_rules! impl_newton_poly_recursive {
    ( $($x:literal, $y:literal), * ) => {
        $(
            paste! {
                pub type [<NewtonPoly $x>]<T> = NewtonPoly<[<NewtonPoly $y >]<T>, T>;

                impl<'a, 'b, S: Display, const P: u64> WithVars<'a, &'b [S; $x]> for [<NewtonPoly $x>]<Z64<P>>
                {
                    type Output = FmtNewtonPoly<'a, 'b, [<NewtonPoly $y >]<Z64<P>>, Z64<P>, S>;

                    fn with_vars(&'a self, vars: &'b[S; $x]) -> Self::Output {
                        FmtNewtonPoly::new(self, vars)
                    }
                }

                impl<const P: u64> From<&[<NewtonPoly $x>]<Z64<P>>> for [<DensePoly $x>]<Z64<P>> {
                    fn from(p: &[<NewtonPoly $x>]<Z64<P>>) -> Self {
                        debug!("convert to poly: {p}");

                        let mut res = Self::zero();
                        let mut prod = Self::one();
                        for (a, y) in &p.coeffs {
                            res += &prod * &DensePoly::from(a);
                            let my = [<DensePoly $y>]::from(-*y);
                            prod *= &Self::from_coeff_unchecked(vec![my, One::one()]);
                        }
                        res + &prod * &DensePoly::from(&p.a_last)
                    }
                }

                impl<const P: u64> TryEval<[Z64<P>; $x]> for [<NewtonPoly $x>]<Z64<P>> {
                    type Output = Z64<P>;

                    fn try_eval(&self, x: &[Z64<P>; $x]) -> Option<Z64<P>> {
                        // TODO: better split_first operating on arrays
                        let (x0, rest) = x.split_first().unwrap();
                        let mut xs = [Z64::zero(); $y];
                        xs.copy_from_slice(&rest);

                        let mut res = self.a_last.eval(&xs);
                        for (a, y) in self.coeffs.iter().rev() {
                            res = a.eval(&xs) + (x0 - y) * res;
                        }
                        Some(res)
                    }
                }
                impl<const P: u64> Eval<[Z64<P>; $x]> for [<NewtonPoly $x>]<Z64<P>> { }

                impl<F, const P: u64> Rec<NewtonRec, [Z64<P>; $x]> for F
                where F: FnMut([Z64<P>; $x]) -> Z64<P> {
                    type Output = Option<[<NewtonPoly $x>]<Z64<P>>>;

                    fn rec_with_ran(
                        &mut self,
                        rec: NewtonRec,
                        mut rng: impl Rng
                    ) -> Self::Output {
                        debug!("{}d reconstruction", $x);
                        let mut rec = [<NewtonPolyRec $x>]::new(rec.extra_pts);
                        let mut y: [Z64<P>; $x] = rng.gen();
                        while let ControlFlow::Continue(n) = rec.add_pt(&y, (self)(y)) {
                            if n > y.len() {
                                return None
                            }
                            y[n] += Z64::one();
                        }
                        Some(rec.into_poly())
                    }
                }

            }
        )*
    };
}

impl_newton_poly!(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
impl_newton_poly_recursive!(
    16, 15, 15, 14, 14, 13, 13, 12, 12, 11, 11, 10, 10, 9, 9, 8, 8, 7, 7, 6, 6,
    5, 5, 4, 4, 3, 3, 2, 2, 1
);
impl_newton_rec_recursive!(
    16, 15, 15, 14, 14, 13, 13, 12, 12, 11, 11, 10, 10, 9, 9, 8, 8, 7, 7, 6, 6,
    5, 5, 4, 4, 3, 3, 2, 2, 1
);

#[cfg(test)]
mod tests {
    use std::iter::repeat_with;

    use rand_xoshiro::rand_core::SeedableRng;

    use super::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn rec_poly_small() {
        log_init();

        const NTESTS: u32 = 10;
        const MAX_POW: u32 = 1;
        const P: u64 = 5;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = NewtonRec::new(1);
        for _ in 0..NTESTS {
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let coeff =
                repeat_with(|| rng.gen::<Z64<P>>()).take(nterms).collect();
            let poly = DensePoly::from_coeff(coeff);
            eprintln!("testing {poly}");
            let reconstructed = (|x: Z64<P>| poly.eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            eprintln!("{reconstructed}");
            let reconstructed: DensePoly<Z64<P>> = reconstructed.into();
            eprintln!("{reconstructed}");
            assert_eq!(poly, reconstructed)
        }
    }

    #[test]
    fn rec_poly_large() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_POW: u32 = 8;
        const P: u64 = 1152921504606846883;
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = NewtonRec::new(1);
        for _ in 0..NTESTS {
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let coeff =
                repeat_with(|| rng.gen::<Z64<P>>()).take(nterms).collect();
            let poly = DensePoly::from_coeff(coeff);
            let reconstructed = (|x: Z64<P>| poly.eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            let reconstructed: DensePoly<Z64<P>> = reconstructed.into();
            assert_eq!(poly, reconstructed)
        }
    }

    #[test]
    fn rec2_poly_small() {
        log_init();

        const NTESTS: u32 = 10;
        const MAX_POW: u32 = 1;
        const P: u64 = 5;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = NewtonRec::new(1);
        for _ in 0..NTESTS {
            let mut coeffs = Vec::new();
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            for _ in 0..nterms {
                let max_pow = rng.gen_range(0..=MAX_POW);
                let nterms = 2usize.pow(max_pow);
                let coeff =
                    repeat_with(|| rng.gen::<Z64<P>>()).take(nterms).collect();
                let poly = DensePoly::from_coeff(coeff);
                coeffs.push(poly);
            }
            let poly = DensePoly::from_coeff(coeffs);
            eprintln!("original: {poly}");
            let reconstructed =
                (|x| poly.eval(&x)).rec_with_ran(rec, &mut rng).unwrap();
            eprintln!("{reconstructed}");
            let reconstructed: DensePoly2<Z64<P>> = reconstructed.into();
            eprintln!("{reconstructed}");
            assert_eq!(poly, reconstructed)
        }
    }

    #[test]
    fn rec2_poly_large() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_POW: u32 = 5;
        const P: u64 = 1152921504606846883;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = NewtonRec::new(1);
        for _ in 0..NTESTS {
            let mut coeffs = Vec::new();
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            for _ in 0..nterms {
                let max_pow = rng.gen_range(0..=MAX_POW);
                let nterms = 2usize.pow(max_pow);
                let coeff =
                    repeat_with(|| rng.gen::<Z64<P>>()).take(nterms).collect();
                let poly = DensePoly::from_coeff(coeff);
                coeffs.push(poly);
            }
            let poly = DensePoly::from_coeff(coeffs);
            let reconstructed =
                (|x| poly.eval(&x)).rec_with_ran(rec, &mut rng).unwrap();
            let reconstructed: DensePoly2<Z64<P>> = reconstructed.into();
            assert_eq!(poly, reconstructed)
        }
    }

    #[test]
    fn rec3_poly_small() {
        log_init();

        const P: u64 = 5;
        let x1: DensePoly3<Z64<P>> =
            DensePoly3::from_coeff(vec![DensePoly2::zero(), DensePoly2::one()]);
        let x2: DensePoly2<Z64<P>> =
            DensePoly::from_coeff(vec![DensePoly::zero(), DensePoly::one()]);
        let x2: DensePoly3<Z64<P>> = DensePoly::from(x2);
        let x3 = DensePoly::from_coeff(vec![Z64::zero(), Z64::one()]);
        let x3: DensePoly2<_> = DensePoly::from(x3);
        let x3: DensePoly3<Z64<P>> = DensePoly::from(x3);
        let poly = x1 + x2 + x3;

        eprintln!("original: {poly}");
        let rec = NewtonRec::new(1);
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let reconstructed =
            (|x| poly.eval(&x)).rec_with_ran(rec, &mut rng).unwrap();
        eprintln!("{reconstructed}");
        let reconstructed: DensePoly3<Z64<P>> = reconstructed.into();
        eprintln!("{reconstructed}");
        assert_eq!(poly, reconstructed)
    }
}
