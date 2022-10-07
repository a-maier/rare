use std::{ops::{AddAssign, SubAssign, Add, Sub, Mul, MulAssign, DivAssign, Div}, fmt::{Display, self}, cmp::min};

use galois_fields::Z64;
use itertools::Itertools;

use crate::traits::{Zero, One, Eval};

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct UniPolynomial<T> {
    coeff: Vec<T>
}

impl<T: Zero> UniPolynomial<T> {
    pub fn new() -> Self {
        Self { coeff: Vec::new() }
    }

    pub fn from_coeff(coeff: Vec<T>) -> Self {
        let mut res = Self { coeff };
        res.delete_trailing_zeroes();
        res
    }

    pub fn from_coeff_unchecked(coeff: Vec<T>) -> Self {
        Self { coeff }
    }

    pub fn into_coeff(self) -> Vec<T> {
        self.coeff
    }

    pub fn degree(&self) -> Option<usize> {
        if self.coeff.is_empty() {
            None
        } else {
            Some(self.coeff.len() - 1)
        }
    }

    pub fn nterms(&self) -> usize {
        self.coeff.iter().filter(|c| !c.is_zero()).count()
    }

    // pub fn eval<'a, 'b, X>(&'a self, x: &'b X) -> X
    // where
    //     X: Zero + Mul<&'b X, Output = X> + Add<&'a T, Output = X>
    // {
    //     self.coeff.iter().rev().fold(
    //         X::zero(),
    //         |acc, c| acc * x + c
    //     )
    // }

    fn delete_trailing_zeroes(&mut self) {
        let last_nonzero = self.coeff.iter()
            .rposition(|c| !c.is_zero())
            .map(|pos| pos + 1)
            .unwrap_or_default();
        self.coeff.truncate(last_nonzero);
    }
}

impl<T: One + Eq> UniPolynomial<T> {
    pub fn is_one(&self) -> bool {
        match self.coeff.first() {
            Some(c) if c.is_one() => true,
            _ => false
        }
    }
}

impl<T: AddAssign + Zero> AddAssign for UniPolynomial<T> {
    fn add_assign(&mut self, mut rhs: Self) {
        if self.coeff.len() < rhs.coeff.len() {
            std::mem::swap(&mut self.coeff, &mut rhs.coeff)
        }
        for (lhs, rhs) in self.coeff.iter_mut().zip(rhs.coeff.into_iter()) {
            *lhs += rhs;
        }
        self.delete_trailing_zeroes()
    }
}

impl<T: Zero> AddAssign<&UniPolynomial<T>> for UniPolynomial<T>
where for<'a> T: AddAssign<&'a T>
{
    fn add_assign(&mut self, rhs: &Self) {
        if self.coeff.len() < rhs.coeff.len() {
            self.coeff.resize_with(rhs.coeff.len(), || T::zero())
        }
        for (lhs, rhs) in self.coeff.iter_mut().zip(rhs.coeff.iter()) {
            *lhs += rhs;
        }
        self.delete_trailing_zeroes()
    }
}

impl<T: SubAssign + Zero> SubAssign for UniPolynomial<T> {
    fn sub_assign(&mut self, mut rhs: Self) {
        if self.coeff.len() < rhs.coeff.len() {
            std::mem::swap(&mut self.coeff, &mut rhs.coeff)
        }
        for (lhs, rhs) in self.coeff.iter_mut().zip(rhs.coeff.into_iter()) {
            *lhs -= rhs;
        }
        self.delete_trailing_zeroes()
    }
}

impl<T: Zero> SubAssign<&UniPolynomial<T>> for UniPolynomial<T>
where for<'a> T: SubAssign<&'a T>
{
    fn sub_assign(&mut self, rhs: &Self) {
        if self.coeff.len() < rhs.coeff.len() {
            self.coeff.resize_with(rhs.coeff.len(), || T::zero())
        }
        for (lhs, rhs) in self.coeff.iter_mut().zip(rhs.coeff.iter()) {
            *lhs -= rhs;
        }
        self.delete_trailing_zeroes()
    }
}

impl<T: AddAssign + Zero> Add for UniPolynomial<T> {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<T: Zero> Add<&UniPolynomial<T>> for UniPolynomial<T>
where for<'a> T: AddAssign<&'a T>
{
    type Output = Self;

    fn add(mut self, rhs: &Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<T: Zero> Add<UniPolynomial<T>> for &UniPolynomial<T>
where for<'a> T: AddAssign<&'a T>
{
    type Output = Self;

    fn add(self, mut rhs: UniPolynomial<T>) -> Self::Output {
        rhs  += self;
        self
    }
}

impl<T: SubAssign + Zero> Sub for UniPolynomial<T> {
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<T: Zero> Sub<&UniPolynomial<T>> for UniPolynomial<T>
where for<'a> T: SubAssign<&'a T>
{
    type Output = Self;

    fn sub(mut self, rhs: &Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<'a, 'b, T: Zero> Mul<&'b UniPolynomial<T>> for &'a UniPolynomial<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>
{
    type Output = UniPolynomial<T>;

    fn mul(self, rhs: &'b UniPolynomial<T>) -> Self::Output {
        let res_len = self.coeff.len() + rhs.coeff.len();
        let mut res = Vec::with_capacity(res_len);
        for n in 0..res_len {
            let start = (n + 1).saturating_sub(rhs.coeff.len());
            let end = min(n + 1, self.coeff.len());
            let mut c = T::zero();
            for i in start..end {
                c += &self.coeff[i] * &rhs.coeff[n - i]
            }
            res.push(c)
        }
        Self::Output::from_coeff(res)
    }
}

impl<T: Zero> Mul<&UniPolynomial<T>> for UniPolynomial<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>
{
    type Output = UniPolynomial<T>;

    fn mul(self, rhs: &UniPolynomial<T>) -> Self::Output {
        &self * rhs
    }
}

impl<T: Zero> Mul<UniPolynomial<T>> for &UniPolynomial<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>
{
    type Output = UniPolynomial<T>;

    fn mul(self, rhs: UniPolynomial<T>) -> Self::Output {
        self * &rhs
    }
}

impl<T: Zero> Mul<UniPolynomial<T>> for UniPolynomial<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>
{
    type Output = UniPolynomial<T>;

    fn mul(self, rhs: UniPolynomial<T>) -> Self::Output {
        &self * &rhs
    }
}

impl<T: Zero> MulAssign for UniPolynomial<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>
{
    fn mul_assign(&mut self, rhs: Self) {
        *self = &*self * rhs;
    }
}

impl<T: Zero> MulAssign<&UniPolynomial<T>> for UniPolynomial<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>
{
    fn mul_assign(&mut self, rhs: &Self) {
        *self = &*self * rhs;
    }
}

impl<T: Zero> MulAssign<&T> for UniPolynomial<T>
where for<'a> T: MulAssign<&'a T>
{
    fn mul_assign(&mut self, rhs: &T) {
        for c in &mut self.coeff {
            *c *= rhs
        }
        self.delete_trailing_zeroes()
    }
}

impl<T: Zero> DivAssign<&T> for UniPolynomial<T>
where for<'a> T: DivAssign<&'a T>
{
    fn div_assign(&mut self, rhs: &T) {
        for c in &mut self.coeff {
            *c /= rhs
        }
        self.delete_trailing_zeroes()
    }
}

impl<T: Zero> Mul<&T> for UniPolynomial<T>
where for<'a> T: MulAssign<&'a T>
{
    type Output = Self;

    fn mul(mut self, rhs: &T) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<'a, 'b, T: Zero + Clone> Mul<&'a T> for &'b UniPolynomial<T>
where for<'c> T: MulAssign<&'c T>
{
    type Output = UniPolynomial<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        self.to_owned() * rhs
    }
}

impl<T: Zero> Div<&T> for UniPolynomial<T>
where for<'a> T: DivAssign<&'a T>
{
    type Output = Self;

    fn div(mut self, rhs: &T) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<T: Zero> Zero for UniPolynomial<T> {
    fn zero() -> Self {
        Self::new()
    }

    fn is_zero(&self) -> bool {
        self.coeff.is_empty()
    }
}

impl<T: One> One for UniPolynomial<T> {
    fn one() -> Self {
        Self { coeff: vec![T::one()] }
    }

    fn is_one(&self) -> bool {
        self.coeff.len() == 1
            && self.coeff[0].is_one()
    }
}

pub type UniPolynomial1<T> = UniPolynomial<T>;
pub type UniPolynomial2<T> = UniPolynomial<UniPolynomial<T>>;
pub type UniPolynomial3<T> = UniPolynomial<UniPolynomial2<T>>;
pub type UniPolynomial4<T> = UniPolynomial<UniPolynomial3<T>>;
pub type UniPolynomial5<T> = UniPolynomial<UniPolynomial4<T>>;

// TODO: macros?
impl<const P: u64> Eval<Z64<P>> for UniPolynomial<Z64<P>> {
    type Output = Z64<P>;

    fn eval(&self, x: &Z64<P>) -> Z64<P> {
        self.coeff.iter().rev().fold(
            Zero::zero(),
            |acc, c| acc * x + c
        )
    }
}

impl<const P: u64> Eval<[Z64<P>; 1]> for UniPolynomial<Z64<P>> {
    type Output = Z64<P>;

    fn eval(&self, x: &[Z64<P>; 1]) -> Z64<P> {
        self.eval(&x[0])
    }
}

impl<const P: u64> Eval<[Z64<P>; 2]> for UniPolynomial2<Z64<P>> {
    type Output = Z64<P>;

    fn eval(&self, x: &[Z64<P>; 2]) -> Z64<P> {
        let (x0, rest) = x.split_first().unwrap();
        let mut xs = [Zero::zero(); 1];
        xs.copy_from_slice(rest);
        self.coeff.iter().rev().fold(
            Zero::zero(),
            |acc, c| acc * x0 + c.eval(&xs)
        )
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct FmtUniPoly<'a, 'b, T: Display + One + Zero, V: Display> {
    poly: &'a UniPolynomial<T>,
    var: &'b[V],
}

impl<'a, 'b, T: Display + One + Zero, V: Display> FmtUniPoly<'a, 'b, T, V> {
    fn new(poly: &'a UniPolynomial<T>, var: &'b [V]) -> Self {
        Self { poly, var }
    }
}

impl<const P: u64> Display for UniPolynomial<Z64<P>> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        FmtUniPoly::new(self, &["x"]).fmt(f)
    }
}

impl<const P: u64> Display for UniPolynomial<UniPolynomial<Z64<P>>> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        FmtUniPoly::new(self, &["x", "y"]).fmt(f)
    }
}

impl<'a, 'b, V: Display, const P: u64> Display for FmtUniPoly<'a, 'b, Z64<P>, V> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.poly.coeff.is_empty() {
            return write!(f, "0");
        }
        let coeff = self.poly.coeff.iter().enumerate().filter(|(_, c)| !c.is_zero());
        let var = &self.var[0];
        write!(f, "{}", coeff.map(|(pow, coeff)| match pow {
            0 => format!("{coeff}"),
            1 => if coeff.is_one() {
                format!("{var}")
            } else {
                format!("{coeff}*{var}")
            },
            _ => if coeff.is_one() {
                format!("{var}^{pow}")
            } else {
                format!("{coeff}*{var}^{pow}")
            }
        }).join(" + "))
    }
}

impl<'a, 'b, V: Display, const P: u64> Display for FmtUniPoly<'a, 'b, UniPolynomial<Z64<P>>, V> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.poly.coeff.is_empty() {
            return write!(f, "0");
        }
        let coeff = self.poly.coeff.iter().enumerate().filter(|(_, c)| !c.is_zero());
        let var = &self.var[0];
        write!(f, "{}", coeff.map(|(pow, coeff)| {
            let coeff = FmtUniPoly::new(coeff, &self.var[1..]);
            match pow {
                0 => format!("{coeff}"),
                1 => if coeff.poly.is_one() {
                    format!("{var}")
                } else if coeff.poly.nterms() == 1 {
                    format!("{coeff}*{var}")
                } else {
                    format!("({coeff})*{var}")
                },
                _ => if coeff.poly.is_one() {
                    format!("{var}^{pow}")
                } else if coeff.poly.nterms() == 1 {
                    format!("{coeff}*{var}^{pow}")
                } else {
                    format!("({coeff})*{var}^{pow}")
                }
            }
        }).join(" + ")
        )
    }
}
