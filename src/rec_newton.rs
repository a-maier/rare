use std::{ops::{Index, IndexMut}, fmt::{Display, self}, cmp::max};

use galois_fields::Z64;
use log::{debug, trace};
use rand::{Rng, thread_rng};

use crate::{traits::{Zero, One}, rand::{UniqueRand, UniqueRandIter}, poly::{UniPolynomial2, UniPolynomial1}};
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

    pub fn rec2_with_ran<F, const P: u64>(
        &self,
        mut poly: F,
        mut rng: impl Rng,
    ) -> Option<NewtonPolynomial2<Z64<P>>>
    where F: FnMut(Z64<P>, Z64<P>) -> Z64<P>
    {
        debug!("2d reconstruction");
        let mut res: NewtonPolynomial2<Z64<P>> = Default::default();
        let mut rand = UniqueRand::new();
        while let Some(x1) = rand.try_gen(&mut rng) {
            trace!("x = {x1}");
            let prefact = res.val.iter().fold(
                Z64::<P>::one(),
                |p, pt| p * (x1 - pt)
            );
            let prefact = prefact.inv();
            let a = self.rec1_with_ran(
                |x2| prefact * (poly(x1, x2) - res.eval(&[x1, x2])),
                &mut rng
            );
            if let Some(a) = a {
                trace!("Reconstructed coefficient {}", FmtNewtonPoly::new(&a, &["y"]));
                res.add_term(a, x1);
                trace!("Intermediate poly: {res}");
                if res.coeff.iter().rev().take_while(|c| c.is_zero()).count() > self.extra_pts {
                    res.coeff.truncate(res.coeff.len() - self.extra_pts - 1);
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

pub type NewtonPolynomial1<T> = NewtonPolynomial<T>;
pub type NewtonPolynomial2<T> = NewtonPolynomial<NewtonPolynomial<T>, T>;
pub type NewtonPolynomial3<T> = NewtonPolynomial<NewtonPolynomial2<T>, T>;
pub type NewtonPolynomial4<T> = NewtonPolynomial<NewtonPolynomial3<T>, T>;
pub type NewtonPolynomial5<T> = NewtonPolynomial<NewtonPolynomial4<T>, T>;

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

impl<const P: u64> From<NewtonPolynomial<Z64<P>>> for UniPolynomial<Z64<P>> {
    fn from(p: NewtonPolynomial<Z64<P>>) -> Self {
        Self::from(&p)
    }
}

impl<const P: u64> From<&NewtonPolynomial2<Z64<P>>> for UniPolynomial2<Z64<P>> {
    fn from(p: &NewtonPolynomial2<Z64<P>>) -> Self {
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
            trace!("({prod}) * ({a}) = {term}");
            res += term;
            let my = UniPolynomial1::from_coeff_unchecked(vec![-*y]);
            prod *= &Self::from_coeff_unchecked(vec![my, One::one()]);
        }
        if let Some(last) = p.coeff.last() {
            let last = UniPolynomial::from(last);
            let term = &prod * &last;
            trace!("({prod}) * ({last}) = {term}");
            res += term;
        }
        res
    }
}

impl<const P: u64> From<NewtonPolynomial2<Z64<P>>> for UniPolynomial2<Z64<P>> {
    fn from(p: NewtonPolynomial2<Z64<P>>) -> Self {
        Self::from(&p)
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
// TODO: macros?
impl<const P: u64> Eval<[Z64<P>; 1]> for NewtonPolynomial<Z64<P>> {
    type Output = Z64<P>;

    fn eval(&self, x: &[Z64<P>; 1]) -> Z64<P> {
        self.eval(&x[0])
    }
}

impl<const P: u64> Eval<[Z64<P>; 2]> for NewtonPolynomial2<Z64<P>> {
    type Output = Z64<P>;

    fn eval(&self, x: &[Z64<P>; 2]) -> Z64<P> {
        if let Some((coeff, coeffs)) = self.coeff.split_first() {
            let mut res = coeff.eval(&x[1]);
            let mut prod: Z64<P> = One::one();
            for (coeff, val) in coeffs.into_iter().zip(self.val.iter()) {
                prod *= x[0] - val;
                res += coeff.eval(&x[1]) * prod;
            }
            res
        } else {
            Zero::zero()
        }
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

impl<const P: u64> Display for NewtonPolynomial<Z64<P>> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        FmtNewtonPoly::new(self, &["x"]).fmt(f)
    }
}

impl<const P: u64> Display for NewtonPolynomial2<Z64<P>> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        FmtNewtonPoly::new(self, &["x", "y"]).fmt(f)
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

impl<'a, 'b, V: Display, const P: u64> Display for FmtNewtonPoly<'a, 'b, NewtonPolynomial<Z64<P>>, Z64<P>, V> {
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
            let reconstructed = rec.rec_univariate_with_ran(|x| poly.eval(&x), &mut rng).unwrap();
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
            let reconstructed = rec.rec_univariate_with_ran(|x| poly.eval(&x), &mut rng).unwrap();
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
            let reconstructed = rec.rec2_with_ran(|x, y| poly.eval(&[x, y]), &mut rng).unwrap();
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
            eprintln!("original: {poly}");
            let reconstructed = rec.rec2_with_ran(|x, y| poly.eval(&[x, y]), &mut rng).unwrap();
            eprintln!("{reconstructed}");
            let reconstructed: UniPolynomial2<Z64<P>> = reconstructed.into();
            eprintln!("{reconstructed}");
            assert_eq!(poly, reconstructed)
        }
    }

}
