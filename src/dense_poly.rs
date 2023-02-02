use std::{ops::{AddAssign, SubAssign, Add, Sub, Mul, MulAssign, DivAssign, Div}, fmt::{Display, self}, cmp::min};

use galois_fields::Z64;
use itertools::Itertools;
use paste::paste;

use crate::traits::{Zero, One, Eval, WithVars};

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct DensePoly<T> {
    coeff: Vec<T>
}

impl<T> DensePoly<T> {
    pub fn new() -> Self {
        Self { coeff: Vec::new() }
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
}

impl<T: Zero> DensePoly<T> {
    pub fn from_coeff(coeff: Vec<T>) -> Self {
        let mut res = Self { coeff };
        res.delete_trailing_zeroes();
        res
    }

    pub fn from_coeff_unchecked(coeff: Vec<T>) -> Self {
        debug_assert!(!matches!(coeff.last(), Some(c) if c.is_zero()));
        Self { coeff }
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

impl<T: AddAssign + Zero> AddAssign for DensePoly<T> {
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

impl<T: Zero> AddAssign<&DensePoly<T>> for DensePoly<T>
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

impl<T: AddAssign + Zero> AddAssign<T> for DensePoly<T> {
    fn add_assign(&mut self, rhs: T) {
        if rhs.is_zero() {
            return;
        }
        if self.coeff.is_empty() {
            self.coeff = vec![rhs];
        } else {
            self.coeff[0] += rhs;
            if self.coeff.len() == 1 && self.coeff[0].is_zero() {
                self.coeff.clear()
            }
        }
    }
}

impl<T: Zero> AddAssign<&T> for DensePoly<T>
where for<'a> T: AddAssign<&'a T>
{
    fn add_assign(&mut self, rhs: &T) {
        if rhs.is_zero() {
            return;
        }
        if self.coeff.is_empty() {
            self.coeff.push(Zero::zero())
        }
        self.coeff[0] += rhs;
        if self.coeff.len() == 1 && self.coeff[0].is_zero() {
            self.coeff.clear()
        }
    }
}

impl<T: SubAssign + Zero> SubAssign for DensePoly<T> {
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

impl<T: Zero> SubAssign<&DensePoly<T>> for DensePoly<T>
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

impl<T: SubAssign + Zero> SubAssign<T> for DensePoly<T> {
    fn sub_assign(&mut self, rhs: T) {
        if rhs.is_zero() {
            return;
        }
        if self.coeff.is_empty() {
            self.coeff.push(Zero::zero())
        }
        self.coeff[0] -= rhs;
        if self.coeff.len() == 1 && self.coeff[0].is_zero() {
            self.coeff.clear()
        }
    }
}

impl<T: Zero> SubAssign<&T> for DensePoly<T>
where for<'a> T: SubAssign<&'a T>
{
    fn sub_assign(&mut self, rhs: &T) {
        if rhs.is_zero() {
            return;
        }
        if self.coeff.is_empty() {
            self.coeff.push(Zero::zero())
        }
        self.coeff[0] -= rhs;
        if self.coeff.len() == 1 && self.coeff[0].is_zero() {
            self.coeff.clear()
        }
    }
}

impl<T: AddAssign + Zero> Add for DensePoly<T> {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<T: Zero> Add<&DensePoly<T>> for DensePoly<T>
where for<'a> T: AddAssign<&'a T>
{
    type Output = Self;

    fn add(mut self, rhs: &Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<T: Zero> Add<DensePoly<T>> for &DensePoly<T>
where for<'a> T: AddAssign<&'a T>
{
    type Output = Self;

    fn add(self, mut rhs: DensePoly<T>) -> Self::Output {
        rhs  += self;
        self
    }
}

impl<T: AddAssign + Zero> Add<T> for DensePoly<T> {
    type Output = Self;

    fn add(mut self, rhs: T) -> Self::Output {
        self += rhs;
        self
    }
}

impl<T: Zero> Add<&T> for DensePoly<T>
where for<'a> T: AddAssign<&'a T>
{
    type Output = Self;

    fn add(mut self, rhs: &T) -> Self::Output {
        self += rhs;
        self
    }
}

impl<T: Zero> Sub<&DensePoly<T>> for DensePoly<T>
where for<'a> T: SubAssign<&'a T>
{
    type Output = Self;

    fn sub(mut self, rhs: &Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<T: Zero> Sub<&T> for DensePoly<T>
where for<'a> T: SubAssign<&'a T>
{
    type Output = Self;

    fn sub(mut self, rhs: &T) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<T: SubAssign + Zero> Sub<T> for DensePoly<T> {
    type Output = Self;

    fn sub(mut self, rhs: T) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<'a, 'b, T: Zero> Mul<&'b DensePoly<T>> for &'a DensePoly<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>
{
    type Output = DensePoly<T>;

    fn mul(self, rhs: &'b DensePoly<T>) -> Self::Output {
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

impl<T: Zero> Mul<&DensePoly<T>> for DensePoly<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>
{
    type Output = DensePoly<T>;

    fn mul(self, rhs: &DensePoly<T>) -> Self::Output {
        &self * rhs
    }
}

impl<T: Zero> Mul<DensePoly<T>> for &DensePoly<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>
{
    type Output = DensePoly<T>;

    fn mul(self, rhs: DensePoly<T>) -> Self::Output {
        self * &rhs
    }
}

impl<T: Zero> Mul<DensePoly<T>> for DensePoly<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>
{
    type Output = DensePoly<T>;

    fn mul(self, rhs: DensePoly<T>) -> Self::Output {
        &self * &rhs
    }
}

impl<T: Zero> MulAssign for DensePoly<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>
{
    fn mul_assign(&mut self, rhs: Self) {
        *self = &*self * rhs;
    }
}

impl<T: Zero> MulAssign<&DensePoly<T>> for DensePoly<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>
{
    fn mul_assign(&mut self, rhs: &Self) {
        *self = &*self * rhs;
    }
}

impl<T: Zero> MulAssign<&T> for DensePoly<T>
where for<'a> T: MulAssign<&'a T>
{
    fn mul_assign(&mut self, rhs: &T) {
        for c in &mut self.coeff {
            *c *= rhs
        }
        self.delete_trailing_zeroes()
    }
}

impl<T: Zero> DivAssign<&T> for DensePoly<T>
where for<'a> T: DivAssign<&'a T>
{
    fn div_assign(&mut self, rhs: &T) {
        for c in &mut self.coeff {
            *c /= rhs
        }
        self.delete_trailing_zeroes()
    }
}

impl<T: Zero> Mul<&T> for DensePoly<T>
where for<'a> T: MulAssign<&'a T>
{
    type Output = Self;

    fn mul(mut self, rhs: &T) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<'a, 'b, T: Zero + Clone> Mul<&'a T> for &'b DensePoly<T>
where for<'c> T: MulAssign<&'c T>
{
    type Output = DensePoly<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        self.to_owned() * rhs
    }
}

impl<T: Zero> Div<&T> for DensePoly<T>
where for<'a> T: DivAssign<&'a T>
{
    type Output = Self;

    fn div(mut self, rhs: &T) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<T> Zero for DensePoly<T> {
    fn zero() -> Self {
        Self::new()
    }

    fn is_zero(&self) -> bool {
        self.coeff.is_empty()
    }
}

impl<T: One> One for DensePoly<T> {
    fn one() -> Self {
        Self { coeff: vec![T::one()] }
    }

    fn is_one(&self) -> bool {
        self.coeff.len() == 1
            && self.coeff[0].is_one()
    }
}

impl<T: Zero> From<T> for DensePoly<T> {
    fn from(z: T) -> Self {
        if z.is_zero() {
            Self::zero()
        } else {
            Self::from_coeff_unchecked(vec![z])
        }
    }
}

impl<const P: u64> Eval<Z64<P>> for DensePoly<Z64<P>> {
    type Output = Z64<P>;

    fn eval(&self, x: &Z64<P>) -> Z64<P> {
        self.coeff.iter().rev().fold(
            Zero::zero(),
            |acc, c| acc * x + c
        )
    }
}

impl<const P: u64> Eval<[Z64<P>; 1]> for DensePoly<Z64<P>> {
    type Output = Z64<P>;

    fn eval(&self, x: &[Z64<P>; 1]) -> Z64<P> {
        self.eval(&x[0])
    }
}

impl<'a, 'b, T, S: Display> WithVars<'a, &'b [S; 1]> for DensePoly<T>
where
    T: Display + One + Zero + 'a,
{
    type Output = FmtUniPoly<'a, 'b, T, S>;

    fn with_vars(&'a self, vars: &'b[S; 1]) -> Self::Output {
        FmtUniPoly::new(self, vars)
    }
}

pub type DensePoly1<T> = DensePoly<T>;

macro_rules! count {
    () => {0usize};
    ($_head:tt $($tail:tt)*) => {1usize + count!($($tail)*)};
}

macro_rules! impl_dense_poly {
    ( $($x:literal), *) => {
        pub(crate) const ALL_VARS: [&str; count!($($x)*)] = [$(
            concat!("x", stringify!($x)),
        )*];

        $(
            paste! {
                impl<const P: u64> Display for [<DensePoly $x>]<Z64<P>> {
                    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                        const VARS: [&str; $x] = {
                            let mut vars = [""; $x];
                            let mut num = 0;
                            while num < vars.len() {
                                vars[num] = ALL_VARS[num];
                                num += 1;
                            }
                            vars
                        };
                        self.with_vars(&VARS).fmt(f)
                    }
                }

                impl<'a, 'b, V: Display, const P: u64> Display for FmtUniPoly<'a, 'b, [<DensePoly $x>]<Z64<P>>, V> {
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
            }
        )*
    };
}

macro_rules! impl_dense_poly_recursive {
    ( $($x:literal, $y:literal), * ) => {
        $(
            paste! {
                pub type [<DensePoly $x>]<T> = DensePoly<[<DensePoly $y>]<T>>;

                impl<const P: u64> From<Z64<P>> for [<DensePoly $x>]<Z64<P>> {
                    fn from(z: Z64<P>) -> Self {
                        if z.is_zero() {
                            Self::zero()
                        } else {
                            Self::from_coeff_unchecked(vec![
                                [<DensePoly $y>]::from(z)
                            ])
                        }
                    }
                }

                impl<const P: u64> Eval<[Z64<P>; $x]> for [<DensePoly $x>]<Z64<P>> {
                    type Output = Z64<P>;

                    fn eval(&self, x: &[Z64<P>; $x]) -> Z64<P> {
                        let (x0, rest) = x.split_first().unwrap();
                        let mut xs = [Zero::zero(); $y];
                        xs.copy_from_slice(rest);
                        self.coeff.iter().rev().fold(
                            Zero::zero(),
                            |acc, c| acc * x0 + c.eval(&xs)
                        )
                    }
                }

                impl<'a, 'b, S: Display, const P: u64> WithVars<'a, &'b [S; $x]> for [<DensePoly $x>]<Z64<P>>
                {
                    type Output = FmtUniPoly<'a, 'b, [<DensePoly $y >]<Z64<P>>, S>;

                    fn with_vars(&'a self, vars: &'b[S; $x]) -> Self::Output {
                        FmtUniPoly::new(self, vars)
                    }
                }
            }
        )*
    };
}

impl_dense_poly!(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
impl_dense_poly_recursive!(16,15,15,14,14,13,13,12,12,11,11,10,10,9,9,8,8,7,7,6,6,5,5,4,4,3,3,2,2,1);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct FmtUniPoly<'a, 'b, T: Display + One + Zero, V: Display> {
    poly: &'a DensePoly<T>,
    var: &'b[V],
}

impl<'a, 'b, T: Display + One + Zero, V: Display> FmtUniPoly<'a, 'b, T, V> {
    fn new(poly: &'a DensePoly<T>, var: &'b [V]) -> Self {
        Self { poly, var }
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
