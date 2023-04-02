use std::{
    cmp::min,
    fmt::{self, Display},
    iter::repeat_with,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign},
};

use galois_fields::Z64;
use itertools::Itertools;
use paste::paste;

use crate::{
    sparse_poly::{SparseMono, SparsePoly},
    traits::{Eval, One, Shift, TryEval, WithVars, Zero},
    util::count,
};

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct DensePoly<T> {
    coeff: Vec<T>,
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

    pub fn len(&self) -> usize {
        self.coeff.len()
    }

    pub fn is_empty(&self) -> bool {
        self.coeff.is_empty()
    }

    pub fn coeffs(&self) -> &[T] {
        &self.coeff
    }

    pub fn coeff(&self, i: usize) -> &T {
        &self.coeff[i]
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

    fn delete_trailing_zeroes(&mut self) {
        let last_nonzero = self
            .coeff
            .iter()
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
where
    for<'a> T: AddAssign<&'a T>,
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
where
    for<'a> T: AddAssign<&'a T>,
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
where
    for<'a> T: SubAssign<&'a T>,
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
where
    for<'a> T: SubAssign<&'a T>,
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
where
    for<'a> T: AddAssign<&'a T>,
{
    type Output = Self;

    fn add(mut self, rhs: &Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<T: Zero> Add<DensePoly<T>> for &DensePoly<T>
where
    for<'a> T: AddAssign<&'a T>,
{
    type Output = Self;

    fn add(self, mut rhs: DensePoly<T>) -> Self::Output {
        rhs += self;
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
where
    for<'a> T: AddAssign<&'a T>,
{
    type Output = Self;

    fn add(mut self, rhs: &T) -> Self::Output {
        self += rhs;
        self
    }
}

impl<T: Zero> Sub<&DensePoly<T>> for DensePoly<T>
where
    for<'a> T: SubAssign<&'a T>,
{
    type Output = Self;

    fn sub(mut self, rhs: &Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<T: Zero> Sub<&T> for DensePoly<T>
where
    for<'a> T: SubAssign<&'a T>,
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
    for<'c> T: AddAssign<<&'c T as Mul>::Output>,
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
    for<'c> T: AddAssign<<&'c T as Mul>::Output>,
{
    type Output = DensePoly<T>;

    fn mul(self, rhs: &DensePoly<T>) -> Self::Output {
        &self * rhs
    }
}

impl<T: Zero> Mul<DensePoly<T>> for &DensePoly<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>,
{
    type Output = DensePoly<T>;

    fn mul(self, rhs: DensePoly<T>) -> Self::Output {
        self * &rhs
    }
}

impl<T: Zero> Mul<DensePoly<T>> for DensePoly<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>,
{
    type Output = DensePoly<T>;

    fn mul(self, rhs: DensePoly<T>) -> Self::Output {
        &self * &rhs
    }
}

impl<T: Zero> MulAssign for DensePoly<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>,
{
    fn mul_assign(&mut self, rhs: Self) {
        *self = &*self * rhs;
    }
}

impl<T: Zero> MulAssign<&DensePoly<T>> for DensePoly<T>
where
    for<'c> &'c T: Mul,
    for<'c> T: AddAssign<<&'c T as Mul>::Output>,
{
    fn mul_assign(&mut self, rhs: &Self) {
        *self = &*self * rhs;
    }
}

impl<T: Zero> MulAssign<&T> for DensePoly<T>
where
    for<'a> T: MulAssign<&'a T>,
{
    fn mul_assign(&mut self, rhs: &T) {
        for c in &mut self.coeff {
            *c *= rhs
        }
        self.delete_trailing_zeroes()
    }
}

impl<T: Zero> DivAssign<&T> for DensePoly<T>
where
    for<'a> T: DivAssign<&'a T>,
{
    fn div_assign(&mut self, rhs: &T) {
        for c in &mut self.coeff {
            *c /= rhs
        }
        self.delete_trailing_zeroes()
    }
}

impl<T: Zero> Mul<&T> for DensePoly<T>
where
    for<'a> T: MulAssign<&'a T>,
{
    type Output = Self;

    fn mul(mut self, rhs: &T) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<'a, 'b, T: Zero + Clone> Mul<&'a T> for &'b DensePoly<T>
where
    for<'c> T: MulAssign<&'c T>,
{
    type Output = DensePoly<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        self.to_owned() * rhs
    }
}

impl<T: Zero> Div<&T> for DensePoly<T>
where
    for<'a> T: DivAssign<&'a T>,
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
        Self {
            coeff: vec![T::one()],
        }
    }

    fn is_one(&self) -> bool {
        self.coeff.len() == 1 && self.coeff[0].is_one()
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

impl<T: Zero> From<SparsePoly<T, 1>> for DensePoly<T> {
    fn from(p: SparsePoly<T, 1>) -> Self {
        let orig_coeff = p.into_terms();
        if let Some(max_pow) = orig_coeff.last().map(|c| c.powers[0]) {
            let mut coeff = Vec::from_iter(
                std::iter::repeat_with(|| T::zero()).take(1 + max_pow as usize),
            );
            for c in orig_coeff {
                coeff[c.powers[0] as usize] = c.coeff;
            }
            Self::from_coeff_unchecked(coeff)
        } else {
            Self::zero()
        }
    }
}

impl<T: Zero> From<SparseMono<T, 1>> for DensePoly<T> {
    fn from(m: SparseMono<T, 1>) -> Self {
        let mut coeff = Vec::from_iter(
            std::iter::repeat_with(|| T::zero()).take(1 + m.powers[0] as usize),
        );
        coeff[m.powers[0] as usize] = m.coeff;
        Self::from_coeff_unchecked(coeff)
    }
}

impl<const P: u64> TryEval<Z64<P>> for DensePoly<Z64<P>> {
    type Output = Z64<P>;

    fn try_eval(&self, x: &Z64<P>) -> Option<Z64<P>> {
        Some(
            self.coeff
                .iter()
                .rev()
                .fold(Zero::zero(), |acc, c| acc * x + c),
        )
    }
}

impl<const P: u64> TryEval<[Z64<P>; 1]> for DensePoly<Z64<P>> {
    type Output = Z64<P>;

    fn try_eval(&self, x: &[Z64<P>; 1]) -> Option<Z64<P>> {
        self.try_eval(&x[0])
    }
}

impl<const P: u64> Eval<Z64<P>> for DensePoly<Z64<P>> {}

impl<const P: u64> Eval<[Z64<P>; 1]> for DensePoly<Z64<P>> {}

// TODO: would be more logical to have consistent input and return types
fn binom<const P: u64>(n: usize, k: usize) -> Z64<P> {
    assert!((n as u64) < P);
    assert!((k as u64) < P);
    let num = ((n - k + 1)..=n)
        .fold(Z64::one(), |acc, i| acc * Z64::new_unchecked(i as u64));
    let den =
        (1..=k).fold(Z64::one(), |acc, i| acc * Z64::new_unchecked(i as u64));
    num / den
}

impl<const P: u64> Shift<Z64<P>> for DensePoly<Z64<P>> {
    fn shift(self, shift: Z64<P>) -> Self {
        if shift.is_zero() {
            return self;
        }
        let mut res_coeff = vec![Z64::zero(); self.len()];
        for (pow, coeff) in self.coeff.into_iter().enumerate().rev() {
            if coeff.is_zero() {
                continue;
            }
            assert!(pow < res_coeff.len());
            let len = pow + 1;
            let mid = len / 2;
            let is_even_pow = (len % 2) == 1;
            for i in 0..mid {
                let binom = binom(pow, i);
                res_coeff[i] += coeff * binom * shift.powu((pow - i) as u64);
                res_coeff[pow - i] += coeff * binom * shift.powu(i as u64);
            }
            if is_even_pow {
                res_coeff[mid] +=
                    coeff * binom(pow, mid) * shift.powu(mid as u64);
            }
        }
        DensePoly::from_coeff(res_coeff)
    }
}

impl<const P: u64> Shift<[Z64<P>; 1]> for DensePoly<Z64<P>> {
    fn shift(self, shift: [Z64<P>; 1]) -> Self {
        self.shift(shift[0])
    }
}

impl<'a, 'b, T, S: Display> WithVars<'a, &'b [S; 1]> for DensePoly<T>
where
    T: Display + One + Zero + 'a,
{
    type Output = FmtUniPoly<'a, 'b, T, S>;

    fn with_vars(&'a self, vars: &'b [S; 1]) -> Self::Output {
        FmtUniPoly::new(self, vars)
    }
}

pub type DensePoly1<T> = DensePoly<T>;

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

                impl<T: AddAssign + Zero> From<SparsePoly<T, $x>> for [<DensePoly $x>]<T> {
                    fn from(p: SparsePoly<T, $x>) -> Self {
                        let orig_coeff = p.into_terms();
                        if let Some(max_pow) = orig_coeff.last().map(|c| c.powers[0]) {
                            let mut res_coeff = Vec::from_iter(
                                std::iter::repeat_with(|| Zero::zero()).take(1 + max_pow as usize)
                            );
                            // TODO: this is not very efficient (lots of allocations)
                            // it's probably better to try and convert all monomials
                            // with the same power[0] at once.
                            for SparseMono{coeff, powers} in orig_coeff {
                                let rest = SparseMono::new(
                                    coeff,
                                    powers[1..].try_into().unwrap()
                                );
                                res_coeff[powers[0] as usize] += [<DensePoly $y>]::<T>::from(rest)
                            }
                            Self::from_coeff_unchecked(res_coeff)
                        } else {
                            Self::zero()
                        }
                    }
                }

                impl<T: Zero> From<SparseMono<T, $x>> for [<DensePoly $x>]<T> {
                    fn from(m: SparseMono<T, $x>) -> Self {
                        let (&pow, rest) = m.powers.split_first().unwrap();
                        let rest: [u32; $y] = rest.try_into().unwrap();
                        let rest = SparseMono::new(m.coeff, rest);
                        let mut res_coeff = Vec::from_iter(
                            std::iter::repeat_with(|| Zero::zero()).take(1 + pow as usize)
                        );
                        res_coeff[pow as usize] = [<DensePoly $y>]::from(rest);
                        Self::from_coeff_unchecked(res_coeff)
                    }
                }

                impl<const P: u64> TryEval<[Z64<P>; $x]> for [<DensePoly $x>]<Z64<P>> {
                    type Output = Z64<P>;

                    fn try_eval(&self, x: &[Z64<P>; $x]) -> Option<Z64<P>> {
                        let (x0, rest) = x.split_first().unwrap();
                        let mut xs = [Zero::zero(); $y];
                        xs.copy_from_slice(rest);
                        Some(self.coeff.iter().rev().fold(
                            Zero::zero(),
                            |acc, c| acc * x0 + c.eval(&xs)
                        ))
                    }
                }

                impl<const P: u64> Eval<[Z64<P>; $x]> for [<DensePoly $x>]<Z64<P>> {
                }

                impl<'a, const P: u64> MulAssign<&'a Z64<P>> for [<DensePoly $x>]<Z64<P>> {
                    fn mul_assign(&mut self, rhs: &'a Z64<P>) {
                        for c in &mut self.coeff {
                            c.mul_assign(rhs);
                        }
                    }
                }

                impl<'a, const P: u64> Mul<&'a Z64<P>> for [<DensePoly $x>]<Z64<P>> {
                    type Output = Self;

                    fn mul(mut self, rhs: &'a Z64<P>) -> Self {
                        self *= rhs;
                        self
                    }
                }

                impl<const P: u64> Shift<[Z64<P>; $x]> for [<DensePoly $x>]<Z64<P>> {
                    fn shift(mut self, shift: [Z64<P>; $x]) -> Self {
                        if shift.is_zero() {
                            return self;
                        }
                        let (this_shift, rest) = shift.split_first().unwrap();
                        let rest: [Z64<P>; $y] = rest.try_into().unwrap();
                        self.coeff = self.coeff.into_iter()
                            .map(|c| c.shift(rest))
                            .collect();
                        if this_shift.is_zero() {
                            return self;
                        }
                        let mut res_coeff = Vec::from_iter(
                            repeat_with(|| [<DensePoly $y>]::zero()).take(self.len())
                        );
                        for (pow, coeff) in self.coeff.into_iter().enumerate().rev() {
                            if coeff.is_zero() {
                                continue;
                            }
                            assert!(pow < res_coeff.len());
                            let len = pow + 1;
                            let mid = len / 2;
                            let is_even_pow = (len % 2) == 1;
                            for i in 0..mid {
                                let binom = binom(pow, i);
                                res_coeff[i] += coeff.clone() * &(
                                    binom * this_shift.powu((pow - i) as u64)
                                );
                                res_coeff[pow - i] += coeff.clone() * &(
                                    binom * this_shift.powu(i as u64)
                                );
                            }
                            if is_even_pow {
                                res_coeff[mid] += coeff * &(
                                    binom(pow, mid) * this_shift.powu(mid as u64)
                                )
                            }
                        }
                        DensePoly::from_coeff(res_coeff)
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

impl_dense_poly!(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
impl_dense_poly_recursive!(
    16, 15, 15, 14, 14, 13, 13, 12, 12, 11, 11, 10, 10, 9, 9, 8, 8, 7, 7, 6, 6,
    5, 5, 4, 4, 3, 3, 2, 2, 1
);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct FmtUniPoly<'a, 'b, T: Display + One + Zero, V: Display> {
    poly: &'a DensePoly<T>,
    var: &'b [V],
}

impl<'a, 'b, T: Display + One + Zero, V: Display> FmtUniPoly<'a, 'b, T, V> {
    fn new(poly: &'a DensePoly<T>, var: &'b [V]) -> Self {
        Self { poly, var }
    }
}

impl<'a, 'b, V: Display, const P: u64> Display
    for FmtUniPoly<'a, 'b, Z64<P>, V>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.poly.coeff.is_empty() {
            return write!(f, "0");
        }
        let coeff = self
            .poly
            .coeff
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero());
        let var = &self.var[0];
        write!(
            f,
            "{}",
            coeff
                .map(|(pow, coeff)| match pow {
                    0 => format!("{coeff}"),
                    1 =>
                        if coeff.is_one() {
                            format!("{var}")
                        } else {
                            format!("{coeff}*{var}")
                        },
                    _ =>
                        if coeff.is_one() {
                            format!("{var}^{pow}")
                        } else {
                            format!("{coeff}*{var}^{pow}")
                        },
                })
                .join(" + ")
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn shift_const() {
        log_init();
        const P: u64 = 29;

        let orig: DensePoly<Z64<P>> = DensePoly::from_coeff(vec![Z64::one()]);

        let unshifted = orig.clone().shift(Z64::zero());
        assert_eq!(orig, unshifted);

        let shifted = orig.clone().shift(Z64::one());
        assert_eq!(orig, shifted);
    }

    #[test]
    fn shift_linear() {
        log_init();
        const P: u64 = 29;

        let orig: DensePoly<Z64<P>> =
            DensePoly::from_coeff(vec![Z64::zero(), Z64::one()]);

        let unshifted = orig.clone().shift(Z64::zero());
        assert_eq!(orig, unshifted);

        let shifted = orig.clone().shift(Z64::one());
        assert_eq!(
            shifted,
            DensePoly::from_coeff(vec![Z64::one(), Z64::one()])
        );
    }

    #[test]
    fn shift_quadratic() {
        log_init();
        const P: u64 = 29;

        let orig: DensePoly<Z64<P>> =
            DensePoly::from_coeff(vec![Z64::zero(), Z64::zero(), Z64::one()]);

        let unshifted = orig.clone().shift(Z64::zero());
        assert_eq!(orig, unshifted);

        let shifted = orig.clone().shift(Z64::one());
        assert_eq!(
            shifted,
            DensePoly::from_coeff(vec![1.into(), 2.into(), 1.into()])
        );

        let shifted = orig.clone().shift(Z64::new(2));
        assert_eq!(
            shifted,
            DensePoly::from_coeff(vec![4.into(), 4.into(), 1.into()])
        );
    }

    #[test]
    fn shift_cubic() {
        log_init();
        const P: u64 = 29;

        let orig: DensePoly<Z64<P>> = DensePoly::from_coeff(vec![
            Z64::zero(),
            Z64::zero(),
            Z64::zero(),
            Z64::one(),
        ]);

        let unshifted = orig.clone().shift(Z64::zero());
        assert_eq!(orig, unshifted);

        let shifted = orig.clone().shift(Z64::one());
        assert_eq!(
            shifted,
            DensePoly::from_coeff(vec![1.into(), 3.into(), 3.into(), 1.into()])
        );
    }

    #[test]
    fn shift_quartic() {
        log_init();
        const P: u64 = 29;

        let mut coeff = vec![Z64::zero(); 5];
        coeff[4] = Z64::one();
        let orig: DensePoly<Z64<P>> = DensePoly::from_coeff(coeff);

        let unshifted = orig.clone().shift(Z64::zero());
        assert_eq!(orig, unshifted);

        let shifted = orig.clone().shift(Z64::one());
        assert_eq!(
            shifted,
            DensePoly::from_coeff(vec![
                1.into(),
                4.into(),
                6.into(),
                4.into(),
                1.into()
            ])
        );
    }

    #[test]
    fn shift_2d() {
        log_init();
        const P: u64 = 29;

        let poly: DensePoly<Z64<P>> =
            DensePoly::from_coeff(vec![0.into(), 1.into()]);
        let orig = DensePoly::from_coeff(vec![Zero::zero(), poly]);

        eprintln!("orig: {orig}");

        let shifted = orig.clone().shift([1.into(), 1.into()]);
        eprintln!("shifted by 1: {shifted}");
        let poly: DensePoly<Z64<P>> =
            DensePoly::from_coeff(vec![1.into(), 1.into()]);
        let ref_shifted = DensePoly::from_coeff(vec![poly.clone(), poly]);
        assert_eq!(shifted, ref_shifted);

        let shifted = orig.clone().shift([2.into(), 2.into()]);
        eprintln!("shifted by 2: {shifted}");
        let poly: DensePoly<Z64<P>> =
            DensePoly::from_coeff(vec![2.into(), 1.into()]);
        let ref_shifted =
            DensePoly::from_coeff(vec![poly.clone() * &Z64::new(2), poly]);
        assert_eq!(shifted, ref_shifted)
    }
}
