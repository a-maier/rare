use std::{ops::{Index, IndexMut}, fmt::{Display, self}, cmp::max};

use galois_fields::Z64;
use log::{debug, trace};
use rand::{Rng, thread_rng};
use paste::paste;

use crate::{traits::{Zero, One, WithVars, Rec}, rand::{UniqueRand, UniqueRandIter}};
use crate::{poly::UniPolynomial, traits::Eval};

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

impl NewtonRec {
    pub fn new(extra_pts: usize) -> Self {
        Self { extra_pts }
    }

    pub fn rec_from_seq<I, const P: u64>(
        &self,
        pts: I
    ) -> Option<NewtonPolynomial<Z64<P>>>
    where I: IntoIterator<Item = (Z64<P>, Z64<P>)>
    {
        debug!("1d reconstruction");
        let mut a = TriangleMatrix::with_row_capacity(self.extra_pts + 1);
        let mut y = Vec::with_capacity(self.extra_pts + 1);
        for (pt, val) in pts.into_iter() {
            trace!("Point ({pt}, {val})");
            let i = a.nrows();
            y.push(pt);
            a.add_zero_row();
            a[(i, 0)] = val;
            for j in 1..=i {
                a[(i,j)] = (a[(i,j-1)] - a[(j-1,j-1)]) / (y[i] - y[j-1]);
            }
            if y.len() > self.extra_pts {
                let mut last_coeffs = (i - self.extra_pts)..=i;
                if last_coeffs.all(|d| a[(d, d)] == Zero::zero()) {
                    let coeff = Vec::from_iter((0..i - self.extra_pts).map(|c| a[(c, c)]));
                    y.truncate(max(coeff.len(), 1) - 1);
                    return Some(NewtonPolynomial {
                        coeff,
                        val: y
                    })
                }
            }
        }
        None
    }

    pub fn rec_univariate_with_ran<F, const P: u64>(
        &self,
        poly: F,
        rng: impl Rng,
    ) -> Option<NewtonPolynomial<Z64<P>>>
    where F: FnMut(Z64<P>) -> Z64<P>
    {
        self.rec1_with_ran(poly, rng)
    }

    pub fn rec1_with_ran<F, const P: u64>(
        &self,
        mut poly: F,
        rng: impl Rng,
    ) -> Option<NewtonPolynomial<Z64<P>>>
    where F: FnMut(Z64<P>) -> Z64<P>
    {
        self.rec_from_seq(
            UniqueRandIter::new(rng).map(|pt| (pt, poly(pt)))
        )
    }

    pub fn rec_univariate<F, const P: u64>(
        &self,
        poly: F
    ) -> Option<NewtonPolynomial<Z64<P>>>
    where F: FnMut(Z64<P>) -> Z64<P>
    {
        self.rec_univariate_with_ran(poly, thread_rng())
    }
}

impl<F, const P: u64> Rec<NewtonRec, (Z64<P>,)> for F
where F: FnMut(Z64<P>) -> Z64<P> {
    type Output = Option<NewtonPolynomial1<Z64<P>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: NewtonRec,
        rng: impl Rng
    ) -> Self::Output {
        reconstructor.rec_from_seq(
            UniqueRandIter::new(rng).map(|pt| (pt, (self)(pt)))
        )
    }
}

impl<F, const P: u64> Rec<NewtonRec, [Z64<P>; 1]> for F
where F: FnMut([Z64<P>; 1]) -> Z64<P> {
    type Output = Option<NewtonPolynomial1<Z64<P>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: NewtonRec,
        rng: impl Rng
    ) -> Self::Output {
        reconstructor.rec_from_seq(
            UniqueRandIter::new(rng).map(|pt| (pt, (self)([pt])))
        )
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct TriangleMatrix<T> {
    elem: Vec<T>,
    nrows: usize
}

impl<T> TriangleMatrix<T> {
    fn with_row_capacity(nrow: usize) -> Self {
        let capacity = row_start(nrow);
        Self {
            elem: Vec::with_capacity(capacity),
            nrows: 0
        }
    }

    fn nrows(&self) -> usize {
        self.nrows
    }
}

impl<T: Clone + Zero> TriangleMatrix<T> {
    fn add_zero_row(&mut self) {
        self.nrows += 1;
        self.elem.resize(self.elem.len() + self.nrows, T::zero())
    }
}

impl<T> Index<(usize, usize)> for TriangleMatrix<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.elem[row_start(index.0) + index.1]
    }
}

impl<T> IndexMut<(usize, usize)> for TriangleMatrix<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        &mut self.elem[row_start(index.0) + index.1]
    }
}

fn row_start(nrow: usize) -> usize {
    (nrow * (nrow + 1)) / 2
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct NewtonPolynomial<T, U = T> {
    coeff: Vec<T>,
    val: Vec<U>,
}

impl<T, U> NewtonPolynomial<T, U> {
    fn new() -> Self {
        NewtonPolynomial {coeff: Vec::new(), val: Vec::new()}
    }
}

impl<T, U> NewtonPolynomial<T, U> {
    fn add_term(&mut self, coeff: T, val: U) {
        self.coeff.push(coeff);
        self.val.push(val);
    }
}

impl<T: Zero, U> Zero for NewtonPolynomial<T, U> {
    fn zero() -> Self {
        Self::new()
    }

    fn is_zero(&self) -> bool {
        self.coeff.is_empty()
    }
}

impl<T: One, U> One for NewtonPolynomial<T, U> {
    fn one() -> Self {
        Self::new()
    }

    fn is_one(&self) -> bool {
        self.coeff.len() == 1
            && self.coeff[0].is_one()
    }
}

impl<const P: u64> From<&NewtonPolynomial<Z64<P>>> for UniPolynomial<Z64<P>> {
    fn from(p: &NewtonPolynomial<Z64<P>>) -> Self {
        if p.coeff.is_empty() {
            return Self::new()
        }
        debug_assert_eq!(p.coeff.len(), p.val.len() + 1);
        let mut prod = Self::one();
        let mut res = Self::new();
        for (a, y) in p.coeff.iter().zip(p.val.iter()) {
            res += &prod * a;
            prod *= Self::from_coeff_unchecked(vec![-*y, One::one()]);
        }
        if let Some(last) = p.coeff.last() {
            res += &prod * last;
        }
        res
    }
}

impl<const P: u64> Eval<Z64<P>> for NewtonPolynomial<Z64<P>> {
    type Output = Z64<P>;

    fn eval(&self, x: &Z64<P>) -> Z64<P> {
        if let Some((coeff, coeffs)) = self.coeff.split_first() {
            let mut res = *coeff;
            let mut prod = One::one();
            for (coeff, val) in coeffs.into_iter().zip(self.val.iter()) {
                prod *= x - val;
                res += coeff * &prod;
            }
            res
        } else {
            Zero::zero()
        }
    }
}

impl<const P: u64> Eval<[Z64<P>; 1]> for NewtonPolynomial<Z64<P>> {
    type Output = Z64<P>;

    fn eval(&self, x: &[Z64<P>; 1]) -> Z64<P> {
        self.eval(&x[0])
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct FmtNewtonPoly<'a, 'b, T: Display + One + Zero, U: Display + One + Zero, V: Display> {
    poly: &'a NewtonPolynomial<T, U>,
    var: &'b [V],
}

impl<'a, 'b, T: Display + One + Zero, U: Display + One + Zero, V: Display> FmtNewtonPoly<'a, 'b, T, U, V> {
    fn new(poly: &'a NewtonPolynomial<T, U>, var: &'b [V]) -> Self {
        Self { poly, var }
    }
}

impl<'a, 'b, V: Display, const P: u64> Display for FmtNewtonPoly<'a, 'b, Z64<P>, Z64<P>, V> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let coeff = &self.poly.coeff;
        let val = &self.poly.val;
        let var = &self.var[0];
        if let Some((first, rest)) = coeff.split_first() {
            debug_assert_eq!(coeff.len(), val.len() + 1);
            write!(f, "{first}")?;
            if let Some((last_coeff, rest)) = rest.split_last() {
                write!(f, " + ")?;
                let (last_pt, pts) = val.split_last().unwrap();
                for (coeff, pt) in rest.iter().zip(pts.iter()) {
                    write!(f, "({var} - {pt})*({coeff} + ")?;
                }
                write!(f, "({var} - {last_pt})*{last_coeff}")?;
                for _ in 0..rest.len() {
                    write!(f, ")")?;
                }
            }
            Ok(())
        } else {
            write!(f, "0")
        }
    }
}

pub type NewtonPolynomial1<T> = NewtonPolynomial<T>;

impl<'a, 'b, T, U, S: Display> WithVars<'a, &'b [S; 1]> for NewtonPolynomial<T, U>
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
                use crate::poly::[<UniPolynomial $x>];

                impl<const P: u64> Display for [<NewtonPolynomial $x>]<Z64<P>> {
                    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                        const VARS: [&str; $x] = {
                            let mut vars = [""; $x];
                            let mut num = 0;
                            while num < vars.len() {
                                vars[num] = crate::poly::ALL_VARS[num];
                                num += 1;
                            }
                            vars
                        };
                        self.with_vars(&VARS).fmt(f)
                    }
                }

                impl<'a, 'b, V: Display, const P: u64> Display for FmtNewtonPoly<'a, 'b, [<NewtonPolynomial $x>]<Z64<P>>, Z64<P>, V> {
                    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                        let coeff = &self.poly.coeff;
                        let val = &self.poly.val;
                        let var = &self.var[0];
                        if let Some((first, rest)) = coeff.split_first() {
                            // debug_assert_eq!(coeff.len(), val.len() + 1);
                            write!(f, "{first}")?;
                            if let Some((last_coeff, rest)) = rest.split_last() {
                                write!(f, " + ")?;
                                let (last_pt, val) = val[..=rest.len()].split_last().unwrap();
                                for (coeff, pt) in rest.iter().zip(val.iter()) {
                                    let coeff = FmtNewtonPoly::new(coeff, &self.var[1..]);
                                    write!(f, "({var} - {pt})*(({coeff}) + ")?;
                                }
                                let last_coeff = FmtNewtonPoly::new(last_coeff, &self.var[1..]);
                                write!(f, "({var} - {last_pt})*({last_coeff})")?;
                                for _ in 0..rest.len() {
                                    write!(f, ")")?;
                                }
                            }
                            Ok(())
                        } else {
                            write!(f, "0")
                        }
                    }
                }

                impl<const P: u64> From<[<NewtonPolynomial $x>]<Z64<P>>> for [<UniPolynomial $x>]<Z64<P>> {
                    fn from(p: [<NewtonPolynomial $x>]<Z64<P>>) -> Self {
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
                pub type [<NewtonPolynomial $x>]<T> = NewtonPolynomial<[<NewtonPolynomial $y >]<T>, T>;

                impl<'a, 'b, S: Display, const P: u64> WithVars<'a, &'b [S; $x]> for [<NewtonPolynomial $x>]<Z64<P>>
                {
                    type Output = FmtNewtonPoly<'a, 'b, [<NewtonPolynomial $y >]<Z64<P>>, Z64<P>, S>;

                    fn with_vars(&'a self, vars: &'b[S; $x]) -> Self::Output {
                        FmtNewtonPoly::new(self, vars)
                    }
                }

                impl<const P: u64> From<&[<NewtonPolynomial $x>]<Z64<P>>> for [<UniPolynomial $x>]<Z64<P>> {
                    fn from(p: &[<NewtonPolynomial $x>]<Z64<P>>) -> Self {
                        debug!("convert to poly: {p}");
                        if p.coeff.is_empty() {
                            return Self::new()
                        }
                        debug_assert_eq!(p.coeff.len(), p.val.len() + 1);
                        let mut prod = Self::one();
                        let mut res = Self::new();
                        for (a, y) in p.coeff.iter().zip(p.val.iter()) {
                            let a = UniPolynomial::from(a);
                            let term = &prod * &a;
                            trace!("({prod}) * ({}) = {term}", Self::from(a));
                            res += term;
                            let my = [<UniPolynomial $y>]::from(-*y);
                            prod *= &Self::from_coeff_unchecked(vec![my, One::one()]);
                        }
                        if let Some(last) = p.coeff.last() {
                            let last = UniPolynomial::from(last);
                            let term = &prod * &last;
                            trace!("({prod}) * ({}) = {term}", Self::from(last));
                            res += term;
                        }
                        res
                    }
                }

                impl<const P: u64> Eval<[Z64<P>; $x]> for [<NewtonPolynomial $x>]<Z64<P>> {
                    type Output = Z64<P>;

                    fn eval(&self, x: &[Z64<P>; $x]) -> Z64<P> {
                        let (x0, rest) = x.split_first().unwrap();
                        if let Some((coeff, coeffs)) = self.coeff.split_first() {
                            let mut xs = [Z64::zero(); $y];
                            xs.copy_from_slice(&rest);
                            let mut res = coeff.eval(&xs);
                            let mut prod: Z64<P> = One::one();
                            for (coeff, val) in coeffs.into_iter().zip(self.val.iter()) {
                                prod *= x0 - val;
                                res += coeff.eval(&xs) * prod;
                            }
                            res
                        } else {
                            Zero::zero()
                        }
                    }
                }

                impl<F, const P: u64> Rec<NewtonRec, [Z64<P>; $x]> for F
                where F: FnMut([Z64<P>; $x]) -> Z64<P> {
                    type Output = Option<[<NewtonPolynomial $x>]<Z64<P>>>;

                    fn rec_with_ran(
                        &mut self,
                        rec: NewtonRec,
                        mut rng: impl Rng
                    ) -> Self::Output {
                        debug!("{}d reconstruction", $x);
                        let mut res: [<NewtonPolynomial $x>]<Z64<P>> = Default::default();
                        let mut rand = UniqueRand::new();
                        while let Some(x1) = rand.try_gen(&mut rng) {
                            trace!("x1 = {x1}");
                            let prefact = res.val.iter().fold(
                                Z64::<P>::one(),
                                |p, pt| p * (x1 - pt)
                            );
                            let prefact = prefact.inv();
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
                            if let Some(a) = a {
                                //trace!("Reconstructed coefficient {}", a.with_vars(&["x2"]));
                                res.add_term(a, x1);
                                trace!("Intermediate poly: {res}");
                                if res.coeff.iter().rev().take_while(|c| c.is_zero()).count() > rec.extra_pts {
                                    res.coeff.truncate(res.coeff.len() - rec.extra_pts - 1);
                                    res.val.truncate(max(res.coeff.len(), 1) - 1);
                                    return Some(res)
                                }
                            } else {
                                return None
                            };
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
            let poly = UniPolynomial::from_coeff(coeff);
            eprintln!("{poly}");
            let reconstructed = (|x: Z64<P>| poly.eval(&x)).rec_with_ran(rec, &mut rng).unwrap();
            eprintln!("{reconstructed}");
            let reconstructed: UniPolynomial<Z64<P>> = reconstructed.into();
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
            let poly = UniPolynomial::from_coeff(coeff);
            let reconstructed = (|x: Z64<P>| poly.eval(&x)).rec_with_ran(rec, &mut rng).unwrap();
            let reconstructed: UniPolynomial<Z64<P>> = reconstructed.into();
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
                let poly = UniPolynomial::from_coeff(coeff);
                coeffs.push(poly);
            }
            let poly = UniPolynomial::from_coeff(coeffs);
            eprintln!("original: {poly}");
            let reconstructed = (|x| poly.eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            eprintln!("{reconstructed}");
            let reconstructed: UniPolynomial2<Z64<P>> = reconstructed.into();
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
                let poly = UniPolynomial::from_coeff(coeff);
                coeffs.push(poly);
            }
            let poly = UniPolynomial::from_coeff(coeffs);
            let reconstructed = (|x| poly.eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            let reconstructed: UniPolynomial2<Z64<P>> = reconstructed.into();
            assert_eq!(poly, reconstructed)
        }
    }

    #[test]
    fn rec3_poly_small() {
        log_init();

        const P: u64 = 5;
        let x1: UniPolynomial3<Z64<P>> = UniPolynomial3::from_coeff(
            vec![UniPolynomial2::zero(), UniPolynomial2::one()]
        );
        let x2: UniPolynomial2<Z64<P>> = UniPolynomial::from_coeff(
            vec![UniPolynomial::zero(), UniPolynomial::one()]
        );
        let x2: UniPolynomial3<Z64<P>> = UniPolynomial::from(x2);
        let x3 = UniPolynomial::from_coeff(vec![Z64::zero(), Z64::one()]);
        let x3: UniPolynomial2<_> = UniPolynomial::from(x3);
        let x3: UniPolynomial3<Z64<P>> = UniPolynomial::from(x3);
        let poly = x1 + x2 + x3;

        eprintln!("original: {poly}");
        let rec = NewtonRec::new(1);
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let reconstructed = (|x| poly.eval(&x))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        eprintln!("{reconstructed}");
        let reconstructed: UniPolynomial3<Z64<P>> = reconstructed.into();
        eprintln!("{reconstructed}");
        assert_eq!(poly, reconstructed)
    }

}
