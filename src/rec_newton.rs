use std::fmt::{Display, self};

use galois_fields::Z64;
use log::{debug, trace};
use rand::{Rng, thread_rng};
use paste::paste;

use crate::{traits::{Zero, One, WithVars, Rec}, rand::pt_iter};
use crate::{dense_poly::DensePoly, traits::{Eval, TryEval}};

/// Univariate polynomial reconstruction using Newton interpolation
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct NewtonRec {
    extra_pts: usize
}

impl Default for NewtonRec {
    fn default() -> Self {
        Self { extra_pts: 1 }
    }
}

fn next_a<const P: u64>(
    poly: &NewtonPoly<Z64<P>>,
    y: Z64<P>,
    f_y: Z64<P>,
    y_last: Z64<P>
) -> Option<Z64<P>> {
    let prefact = poly.coeffs.iter()
        .map(|(_ai, yi)| y - yi)
        .fold(y - y_last, |acc, x| acc * x);
    let prefact = prefact.try_inv()?;
    let poly_at_y = poly.eval(&y);
    Some(prefact * (f_y - poly_at_y))
}

impl NewtonRec {
    pub fn new(extra_pts: usize) -> Self {
        Self { extra_pts }
    }

    pub fn rec_from_seq<I, const P: u64>(
        &self,
        pts: I
    ) -> Option<NewtonPoly<Z64<P>>>
    where I: IntoIterator<Item = (Z64<P>, Z64<P>)>
    {
        debug!("1d polynomial reconstruction");
        let mut pts = pts.into_iter();
        let Some((y0, a0)) = pts.next() else {
            return None
        };
        trace!("Adding p({y0}) = {a0}");
        let mut poly = NewtonPoly::from(a0);
        trace!("p(x1) = {poly}");
        let mut y_last = y0;
        for (y, f_y) in pts {
            trace!("Adding p({y}) = {f_y}");
            let Some(a) = next_a(&poly, y, f_y, y_last) else {
                continue
            };
            poly.coeffs.push((poly.a_last, y_last));
            poly.a_last = a;
            trace!("p(x1) = {poly}");
            if poly.a_last.is_zero() {
                let trailing_zeroes = poly.coeffs.iter().rev()
                    .take_while(|(a, _)| a.is_zero())
                    .count();
                if trailing_zeroes >= self.extra_pts {
                    poly.coeffs.truncate(poly.coeffs.len() - trailing_zeroes);
                    if let Some((a, _y)) = poly.coeffs.pop() {
                        poly.a_last = a;
                    }
                    debug!("Reconstructed {poly}");
                    return Some(poly);
                }
            }
            y_last = y;
        }
        debug!("Reconstruction failed");
        None
    }

    pub fn rec_univariate_with_ran<F, const P: u64>(
        &self,
        poly: F,
        rng: impl Rng,
    ) -> Option<NewtonPoly<Z64<P>>>
    where F: FnMut(Z64<P>) -> Z64<P>
    {
        self.rec1_with_ran(poly, rng)
    }

    pub fn rec1_with_ran<F, const P: u64>(
        &self,
        mut poly: F,
        rng: impl Rng,
    ) -> Option<NewtonPoly<Z64<P>>>
    where F: FnMut(Z64<P>) -> Z64<P>
    {
        self.rec_from_seq(
            pt_iter(rng).map(|pt| (pt, poly(pt)))
        )
    }

    pub fn rec_univariate<F, const P: u64>(
        &self,
        poly: F
    ) -> Option<NewtonPoly<Z64<P>>>
    where F: FnMut(Z64<P>) -> Z64<P>
    {
        self.rec_univariate_with_ran(poly, thread_rng())
    }
}

impl<F, const P: u64> Rec<NewtonRec, (Z64<P>,)> for F
where F: FnMut(Z64<P>) -> Z64<P> {
    type Output = Option<NewtonPoly1<Z64<P>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: NewtonRec,
        rng: impl Rng
    ) -> Self::Output {
        reconstructor.rec_from_seq(
            pt_iter(rng).map(|pt| (pt, (self)(pt)))
        )
    }
}

impl<F, const P: u64> Rec<NewtonRec, [Z64<P>; 1]> for F
where F: FnMut([Z64<P>; 1]) -> Z64<P> {
    type Output = Option<NewtonPoly1<Z64<P>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: NewtonRec,
        rng: impl Rng
    ) -> Self::Output {
        reconstructor.rec_from_seq(
            pt_iter(rng).map(|pt| (pt, (self)([pt])))
        )
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
            coeffs: Vec::new()
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
            coeffs: Vec::new()
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
        let mut res = DensePoly::from_coeff_unchecked(
            vec![p.a_last]
        );
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
impl<const P: u64> Eval<Z64<P>> for NewtonPoly<Z64<P>> { }

impl<const P: u64> TryEval<[Z64<P>; 1]> for NewtonPoly<Z64<P>> {
    type Output = Z64<P>;

    fn try_eval(&self, x: &[Z64<P>; 1]) -> Option<Z64<P>> {
        self.try_eval(&x[0])
    }
}
impl<const P: u64> Eval<[Z64<P>; 1]> for NewtonPoly<Z64<P>> { }

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct FmtNewtonPoly<'a, 'b, T: Display + One + Zero, U: Display + One + Zero, V: Display> {
    poly: &'a NewtonPoly<T, U>,
    var: &'b [V],
}

impl<'a, 'b, T: Display + One + Zero, U: Display + One + Zero, V: Display> FmtNewtonPoly<'a, 'b, T, U, V> {
    fn new(poly: &'a NewtonPoly<T, U>, var: &'b [V]) -> Self {
        Self { poly, var }
    }
}

impl<'a, 'b, V: Display, const P: u64> Display for FmtNewtonPoly<'a, 'b, Z64<P>, Z64<P>, V> {
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

    fn with_vars(&'a self, vars: &'b[S; 1]) -> Self::Output {
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
                        const VARS: [&str; $x] = {
                            let mut vars = [""; $x];
                            let mut num = 0;
                            while num < vars.len() {
                                vars[num] = crate::dense_poly::ALL_VARS[num];
                                num += 1;
                            }
                            vars
                        };
                        self.with_vars(&VARS).fmt(f)
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
                        let mut res: [<NewtonPoly $x>]<Z64<P>> = Default::default();
                        let rand: Z64<P> = rng.gen();
                        for shift in 0..P {
                            let x1 = rand + Z64::new_unchecked(shift);
                            trace!("x1 = {x1}");
                            let prefact = res.coeffs.iter().fold(
                                Z64::<P>::one(),
                                |p, pt| p * (x1 - pt.1)
                            );
                            let Some(prefact) = prefact.try_inv() else {
                                continue
                            };
                            let mut poly_rest = |xs: [Z64<P>; $y]| {
                                let mut args = [x1; $x];
                                args[1..].copy_from_slice(&xs);
                                prefact * ((self)(args) - res.eval(&args))
                            };
                            let a = Rec::<NewtonRec, [Z64<P>; $y]>::rec_with_ran(
                                &mut poly_rest,
                                rec,
                                &mut rng
                            );
                            let Some(a) = a  else {
                                return None
                            };
                            trace!("Reconstructed coefficient {}", a.with_vars(&["x2"]));
                            res.coeffs.push((a, x1));
                            trace!("Intermediate poly: {res}");
                            let trailing_zeroes = res.coeffs.iter().rev().take_while(|(c, _)| c.is_zero()).count();
                            if trailing_zeroes > rec.extra_pts {
                                res.coeffs.truncate(res.coeffs.len() - trailing_zeroes);
                                if let Some((a, _y)) = res.coeffs.pop() {
                                    res.a_last = a;
                                }
                                return Some(res)
                            }
                        }
                        None
                    }
                }

            }
        )*
    };
}

impl_newton_poly!(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
impl_newton_poly_recursive!(16,15,15,14,14,13,13,12,12,11,11,10,10,9,9,8,8,7,7,6,6,5,5,4,4,3,3,2,2,1);

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
            let coeff = repeat_with(|| rng.gen::<Z64<P>>()).take(nterms).collect();
            let poly = DensePoly::from_coeff(coeff);
            eprintln!("testing {poly}");
            let reconstructed = (|x: Z64<P>| poly.eval(&x)).rec_with_ran(rec, &mut rng).unwrap();
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
            let coeff = repeat_with(|| rng.gen::<Z64<P>>()).take(nterms).collect();
            let poly = DensePoly::from_coeff(coeff);
            let reconstructed = (|x: Z64<P>| poly.eval(&x)).rec_with_ran(rec, &mut rng).unwrap();
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
                let coeff = repeat_with(|| rng.gen::<Z64<P>>()).take(nterms).collect();
                let poly = DensePoly::from_coeff(coeff);
                coeffs.push(poly);
            }
            let poly = DensePoly::from_coeff(coeffs);
            eprintln!("original: {poly}");
            let reconstructed = (|x| poly.eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
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
                let coeff = repeat_with(|| rng.gen::<Z64<P>>()).take(nterms).collect();
                let poly = DensePoly::from_coeff(coeff);
                coeffs.push(poly);
            }
            let poly = DensePoly::from_coeff(coeffs);
            let reconstructed = (|x| poly.eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            let reconstructed: DensePoly2<Z64<P>> = reconstructed.into();
            assert_eq!(poly, reconstructed)
        }
    }

    #[test]
    fn rec3_poly_small() {
        log_init();

        const P: u64 = 5;
        let x1: DensePoly3<Z64<P>> = DensePoly3::from_coeff(
            vec![DensePoly2::zero(), DensePoly2::one()]
        );
        let x2: DensePoly2<Z64<P>> = DensePoly::from_coeff(
            vec![DensePoly::zero(), DensePoly::one()]
        );
        let x2: DensePoly3<Z64<P>> = DensePoly::from(x2);
        let x3 = DensePoly::from_coeff(vec![Z64::zero(), Z64::one()]);
        let x3: DensePoly2<_> = DensePoly::from(x3);
        let x3: DensePoly3<Z64<P>> = DensePoly::from(x3);
        let poly = x1 + x2 + x3;

        eprintln!("original: {poly}");
        let rec = NewtonRec::new(1);
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let reconstructed = (|x| poly.eval(&x))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        eprintln!("{reconstructed}");
        let reconstructed: DensePoly3<Z64<P>> = reconstructed.into();
        eprintln!("{reconstructed}");
        assert_eq!(poly, reconstructed)
    }

}
