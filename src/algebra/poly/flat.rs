use std::{
    cmp::Ordering,
    fmt::{self, Display},
    ops::{
        Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
    },
};

use ffnt::Z64;
use num_traits::Pow;
use paste::paste;
use rug::{integer::IntegerExt64, Complete, Integer};

use crate::{
    traits::{Eval, One, TryEval, WithVars, Zero},
    util::{slice_start, ALL_VARS},
};

use super::dense::DensePoly;

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct FlatPoly<T, const Z: usize> {
    terms: Vec<FlatMono<T, Z>>,
}

impl<T, const Z: usize> FlatPoly<T, Z> {
    pub fn new() -> Self {
        Self { terms: Vec::new() }
    }

    pub fn terms(&self) -> &[FlatMono<T, Z>] {
        self.terms.as_ref()
    }

    pub fn len(&self) -> usize {
        self.terms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn into_terms(self) -> Vec<FlatMono<T, Z>> {
        self.terms
    }

    pub fn term(&self, i: usize) -> &FlatMono<T, Z> {
        &self.terms[i]
    }
}

impl<T, const Z: usize> FlatPoly<T, Z>
where
    FlatPoly<T, Z>: One + Clone,
    FlatMono<T, Z>: One + Clone,
    for<'a> FlatMono<T, Z>: MulAssign<&'a FlatMono<T, Z>> + MulAssign<&'a T>,
    for<'a> &'a T: Pow<u32, Output = T>,
    T: From<u64> + Zero,
{
    pub fn powu(&self, n: u32) -> Self {
        if n.is_zero() {
            return Self::one();
        }
        if n.is_one() {
            return self.clone();
        }
        match self.len() {
            0 => self.clone(),
            1 => self.term(0).powu(n).into(),
            2 => {
                let a = self.term(0);
                let b = self.term(1);
                let mut res_terms = Vec::from_iter((0..=n).scan(
                    FlatMono::<T, Z>::one(),
                    |b_to_k, _| {
                        let res = b_to_k.clone();
                        *b_to_k *= b;
                        Some(res)
                    },
                ));
                let mut binom_coeff = 1;
                let mut a_to_l = FlatMono::one();
                let n = n as u64;
                for (k, term) in res_terms.iter_mut().rev().enumerate() {
                    let k = k as u64;
                    *term *= &a_to_l;
                    *term *= &T::from(binom_coeff);
                    a_to_l *= a;
                    binom_coeff = (binom_coeff * (n - k)) / (k + 1);
                }
                Self::from_raw_terms(res_terms)
            }
            _ => todo!("integer power for polynomial with more than two terms"),
        }
    }
}

impl<T: Zero, const Z: usize> FlatPoly<T, Z> {
    pub fn from_raw_terms(terms: Vec<FlatMono<T, Z>>) -> Self {
        debug_assert!(!terms.iter().any(|c| c.is_zero()));
        Self { terms }
    }
}

impl<T: AddAssign + Zero, const Z: usize> FlatPoly<T, Z> {
    pub fn from_terms(mut terms: Vec<FlatMono<T, Z>>) -> Self {
        terms.sort_unstable_by_key(|t| t.powers);
        let mut terms = terms.into_iter();
        let mut res_terms = Vec::with_capacity(terms.len());
        if let Some(first) = terms.next() {
            res_terms.push(first);
            for t in terms {
                let last_idx = res_terms.len() - 1;
                if t.powers == res_terms[last_idx].powers {
                    res_terms[last_idx].coeff += t.coeff
                } else if res_terms[last_idx].is_zero() {
                    res_terms[last_idx] = t;
                } else {
                    res_terms.push(t);
                }
            }
            if res_terms[res_terms.len() - 1].is_zero() {
                res_terms.pop();
            }
        }
        Self::from_raw_terms(res_terms)
    }
}

impl<T: Zero, const Z: usize> From<FlatMono<T, Z>> for FlatPoly<T, Z> {
    fn from(source: FlatMono<T, Z>) -> Self {
        if source.is_zero() {
            Self::zero()
        } else {
            Self::from_raw_terms(vec![source])
        }
    }
}

impl<T, const Z: usize> AddAssign for FlatPoly<T, Z>
where
    FlatPoly<T, Z>: Add<Output = Self>,
{
    fn add_assign(&mut self, rhs: FlatPoly<T, Z>) {
        *self = std::mem::replace(self, Self::new()) + rhs;
    }
}

impl<T, const Z: usize> Add for FlatPoly<T, Z>
where
    T: Zero + AddAssign,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut res = Vec::with_capacity(self.len() + rhs.len());
        let mut lhs_terms = self.terms.into_iter().peekable();
        let mut rhs_terms = rhs.terms.into_iter().peekable();
        loop {
            match (lhs_terms.peek(), rhs_terms.peek()) {
                (None, _) => {
                    res.extend(rhs_terms);
                    break;
                }
                (_, None) => {
                    res.extend(lhs_terms);
                    break;
                }
                (Some(lhs), Some(rhs)) => match lhs.powers.cmp(&rhs.powers) {
                    Ordering::Less => res.push(lhs_terms.next().unwrap()),
                    Ordering::Greater => res.push(rhs_terms.next().unwrap()),
                    Ordering::Equal => {
                        let mut lhs = lhs_terms.next().unwrap();
                        let rhs = rhs_terms.next().unwrap();
                        lhs.coeff += rhs.coeff;
                        if !lhs.coeff.is_zero() {
                            res.push(lhs)
                        }
                    }
                },
            }
        }
        Self::from_raw_terms(res)
    }
}

impl<T, const Z: usize> SubAssign for FlatPoly<T, Z>
where
    FlatPoly<T, Z>: AddAssign,
    for<'a> &'a T: Neg<Output = T>,
{
    #[allow(clippy::suspicious_op_assign_impl)]
    fn sub_assign(&mut self, mut rhs: FlatPoly<T, Z>) {
        for term in &mut rhs.terms {
            term.coeff = term.coeff.neg();
        }
        *self += rhs;
    }
}

impl<T, const Z: usize> Sub for FlatPoly<T, Z>
where
    FlatPoly<T, Z>: SubAssign,
{
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<'a, 'b, T, const Z: usize> Mul<&'b FlatPoly<T, Z>> for &'a FlatPoly<T, Z>
where
    T: Zero + AddAssign,
    &'a T: Mul<&'b T, Output = T>,
{
    type Output = FlatPoly<T, Z>;

    fn mul(self, rhs: &'b FlatPoly<T, Z>) -> Self::Output {
        let mut res = Vec::with_capacity(self.len() * rhs.len());
        for lhs in self.terms() {
            for rhs in rhs.terms() {
                res.push(lhs * rhs)
            }
        }
        Self::Output::from_terms(res)
    }
}

impl<T, const Z: usize> MulAssign<T> for FlatPoly<T, Z>
where
    T: Copy + MulAssign + Zero,
{
    fn mul_assign(&mut self, rhs: T) {
        if rhs.is_zero() {
            self.terms.clear();
        } else {
            for term in &mut self.terms {
                *term *= rhs;
                debug_assert!(!term.is_zero())
            }
        }
    }
}

impl<'a, T, const Z: usize> MulAssign<&'a T> for FlatPoly<T, Z>
where
    T: MulAssign<&'a T> + Zero,
{
    fn mul_assign(&mut self, rhs: &'a T) {
        if rhs.is_zero() {
            self.terms.clear();
        } else {
            for term in &mut self.terms {
                *term *= rhs;
                debug_assert!(!term.is_zero())
            }
        }
    }
}

impl<T, const Z: usize> Mul<T> for FlatPoly<T, Z>
where
    T: Copy + Mul<T> + Zero,
    <T as Mul<T>>::Output: Zero,
{
    type Output = FlatPoly<<T as Mul<T>>::Output, Z>;

    fn mul(self, rhs: T) -> Self::Output {
        if rhs.is_zero() {
            Zero::zero()
        } else {
            Self::Output::from_raw_terms(
                self.into_terms().into_iter().map(|t| t * rhs).collect(),
            )
        }
    }
}

impl<'a, T, const Z: usize> Mul<&'a T> for FlatPoly<T, Z>
where
    T: Mul<&'a T> + Zero,
    <T as Mul<&'a T>>::Output: Zero,
{
    type Output = FlatPoly<<T as Mul<&'a T>>::Output, Z>;

    fn mul(self, rhs: &'a T) -> Self::Output {
        if rhs.is_zero() {
            Zero::zero()
        } else {
            Self::Output::from_raw_terms(
                self.into_terms().into_iter().map(|t| t * rhs).collect(),
            )
        }
    }
}

impl<T, const Z: usize> DivAssign<T> for FlatPoly<T, Z>
where
    T: Copy + DivAssign + Zero,
{
    fn div_assign(&mut self, rhs: T) {
        for term in &mut self.terms {
            *term /= rhs;
            debug_assert!(!term.is_zero())
        }
    }
}

impl<'a, T, const Z: usize> DivAssign<&'a T> for FlatPoly<T, Z>
where
    T: DivAssign<&'a T> + Zero,
{
    fn div_assign(&mut self, rhs: &'a T) {
        for term in &mut self.terms {
            *term /= rhs;
            debug_assert!(!term.is_zero())
        }
    }
}

impl<T, const Z: usize> Div<T> for FlatPoly<T, Z>
where
    T: Copy + Div<T>,
    <T as Div<T>>::Output: Zero,
{
    type Output = FlatPoly<<T as Div<T>>::Output, Z>;

    fn div(self, rhs: T) -> Self::Output {
        Self::Output::from_raw_terms(
            self.into_terms().into_iter().map(|t| t / rhs).collect(),
        )
    }
}

impl<'a, T, const Z: usize> Div<&'a T> for FlatPoly<T, Z>
where
    T: Div<&'a T>,
    <T as Div<&'a T>>::Output: Zero,
{
    type Output = FlatPoly<<T as Div<&'a T>>::Output, Z>;

    fn div(self, rhs: &'a T) -> Self::Output {
        Self::Output::from_raw_terms(
            self.into_terms().into_iter().map(|t| t / rhs).collect(),
        )
    }
}

impl<T, const Z: usize> Zero for FlatPoly<T, Z> {
    fn zero() -> Self {
        Self::new()
    }

    fn is_zero(&self) -> bool {
        self.terms().is_empty()
    }
}

impl<T: Zero + One, const Z: usize> One for FlatPoly<T, Z> {
    fn one() -> Self {
        Self::from_raw_terms(vec![FlatMono::<T, Z>::one()])
    }

    fn is_one(&self) -> bool {
        if let Some(first) = self.terms().first() {
            self.len() == 1 && first.is_one()
        } else {
            false
        }
    }
}

impl<T, const Z: usize> TryEval<[T; Z]> for FlatPoly<T, Z>
where
    FlatMono<T, Z>: Eval<[T; Z], Output = T>,
    T: Add<Output = T> + Zero,
{
    type Output = T;

    fn try_eval(&self, x: &[T; Z]) -> Option<Self::Output> {
        Some(
            self.terms()
                .iter()
                .fold(Zero::zero(), |acc, t| acc + t.eval(x)),
        )
    }
}

impl<T, const Z: usize> Eval<[T; Z]> for FlatPoly<T, Z> where
    FlatPoly<T, Z>: TryEval<[T; Z]>
{
}

impl<const P: u64, const Z: usize> TryEval<[Z64<P>; Z]>
    for FlatPoly<Integer, Z>
{
    type Output = Z64<P>;

    fn try_eval(&self, x: &[Z64<P>; Z]) -> Option<Self::Output> {
        Some(
            self.terms()
                .iter()
                .fold(Zero::zero(), |acc, t| acc + t.eval(x)),
        )
    }
}

impl<const P: u64, const Z: usize> Eval<[Z64<P>; Z]> for FlatPoly<Integer, Z> {}

impl<T: Zero> From<DensePoly<T>> for FlatPoly<T, 1> {
    fn from(source: DensePoly<T>) -> Self {
        let terms = source
            .into_coeff()
            .into_iter()
            .enumerate()
            .filter_map(|(n, c)| {
                if c.is_zero() {
                    None
                } else {
                    Some(FlatMono::new(c, [n as u32]))
                }
            })
            .collect();
        Self { terms }
    }
}
macro_rules! impl_poly_conv_rec {
    ( $($x:literal, $y:literal), * ) => {
        $(
            impl<T: Zero> From<DensePoly<FlatPoly<T, $y>>> for FlatPoly<T, $x> {
                fn from(source: DensePoly<FlatPoly<T, $y>>) -> Self {
                    let terms = source
                        .into_coeff()
                        .into_iter()
                        .enumerate()
                        .filter(|(_, c)| !c.is_zero())
                        .flat_map(|(n, c)| c.into_terms().into_iter().map(move |c| (n, c)))
                        .map(|(n, c)| {
                            let mut pow = [0; $x];
                            pow[1..].copy_from_slice(&c.powers);
                            pow[0] = n as u32;
                            FlatMono::new(c.coeff, pow)
                        })
                        .collect();
                    Self { terms }
                }
            }

            paste! {
                use crate::algebra::poly::dense::[<DensePoly $x>];

                impl<T: Zero> From<[<DensePoly $x>]<T>> for FlatPoly<T, $x> {
                    fn from(source: [<DensePoly $x>]<T>) -> Self {
                        let tmp: DensePoly<FlatPoly<T, $y>> = DensePoly::from_coeff_unchecked(
                            source.into_coeff().into_iter()
                                .map(|c| c.into())
                                .collect()
                        );
                        tmp.into()
                    }
                }
            }
        )*
    };
}

impl_poly_conv_rec!(
    16, 15, 15, 14, 14, 13, 13, 12, 12, 11, 11, 10, 10, 9, 9, 8, 8, 7, 7, 6, 6,
    5, 5, 4, 4, 3, 3, 2, 2, 1
);

pub struct FmtFlatPoly<'a, 'b, V, T, const Z: usize> {
    p: &'a FlatPoly<T, Z>,
    vars: &'b [V; Z],
}

impl<'a, 'b, V: Display + 'b, T: 'a, const Z: usize> WithVars<'a, &'b [V; Z]>
    for FlatPoly<T, Z>
{
    type Output = FmtFlatPoly<'a, 'b, V, T, Z>;

    fn with_vars(&'a self, vars: &'b [V; Z]) -> Self::Output {
        FmtFlatPoly { p: self, vars }
    }
}

impl<'a, 'b, V: Display, const P: u64, const Z: usize> Display
    for FmtFlatPoly<'a, 'b, V, Z64<P>, Z>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some((first, rest)) = self.p.terms().split_first() {
            first.with_vars(self.vars).fmt(f)?;
            for term in rest {
                write!(f, " + {}", term.with_vars(self.vars))?
            }
            Ok(())
        } else {
            write!(f, "0")
        }
    }
}

impl<const P: u64, const Z: usize> Display for FlatPoly<Z64<P>, Z> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let vars: &[_; Z] = slice_start(&ALL_VARS);
        self.with_vars(vars).fmt(f)
    }
}

impl<'a, 'b, V: Display, const Z: usize> Display
    for FmtFlatPoly<'a, 'b, V, Integer, Z>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some((first, rest)) = self.p.terms().split_first() {
            first.with_vars(self.vars).fmt(f)?;
            for term in rest {
                if term.coeff.is_positive() {
                    write!(f, " + {}", term.with_vars(self.vars))?
                } else {
                    let term =
                        FlatMono::new((-&term.coeff).complete(), term.powers);
                    write!(f, " - {}", term.with_vars(self.vars))?
                }
            }
            Ok(())
        } else {
            write!(f, "0")
        }
    }
}

impl<const Z: usize> Display for FlatPoly<Integer, Z> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut vars = [""; Z];
        vars.copy_from_slice(&ALL_VARS[..Z]);
        self.with_vars(&vars).fmt(f)
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct FlatMono<T, const Z: usize> {
    pub powers: [u32; Z],
    pub coeff: T,
}

impl<T, const Z: usize> FlatMono<T, Z> {
    pub fn new(coeff: T, powers: [u32; Z]) -> Self {
        Self { coeff, powers }
    }
}

impl<T, const Z: usize> FlatMono<T, Z>
where
    for<'a> &'a T: Pow<u32, Output = T>,
{
    fn powu(&self, k: u32) -> Self {
        Self {
            powers: self.powers.map(|p| p * k),
            coeff: self.coeff.pow(k),
        }
    }
}

impl<T: Default, const Z: usize> Default for FlatMono<T, Z> {
    fn default() -> Self {
        Self {
            coeff: Default::default(),
            powers: [Default::default(); Z],
        }
    }
}

impl<T: AddAssign + Zero, const Z: usize> Add for FlatMono<T, Z> {
    type Output = FlatPoly<T, Z>;

    fn add(self, rhs: Self) -> Self::Output {
        Self::Output::from_terms(vec![self, rhs])
    }
}

impl<T: AddAssign + Zero, const Z: usize> Sub for FlatMono<T, Z>
where
    FlatMono<T, Z>: Neg<Output = Self>,
{
    type Output = FlatPoly<T, Z>;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output::from_terms(vec![self, rhs.neg()])
    }
}

impl<T, const Z: usize> MulAssign for FlatMono<T, Z>
where
    T: MulAssign,
{
    fn mul_assign(&mut self, rhs: FlatMono<T, Z>) {
        self.coeff *= rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
    }
}

impl<'a, T, const Z: usize> MulAssign<&'a FlatMono<T, Z>> for FlatMono<T, Z>
where
    T: MulAssign<&'a T>,
{
    fn mul_assign(&mut self, rhs: &'a FlatMono<T, Z>) {
        self.coeff *= &rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.iter()) {
            *a += *b;
        }
    }
}

impl<T, const Z: usize> Mul<FlatMono<T, Z>> for FlatMono<T, Z>
where
    T: Mul,
{
    type Output = FlatMono<<T as Mul>::Output, Z>;

    fn mul(mut self, rhs: FlatMono<T, Z>) -> Self::Output {
        let coeff = self.coeff * rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
        Self::Output::new(coeff, self.powers)
    }
}

impl<'a, T, const Z: usize> Mul<FlatMono<T, Z>> for &'a FlatMono<T, Z>
where
    &'a T: Mul<T>,
{
    type Output = FlatMono<<&'a T as Mul<T>>::Output, Z>;

    fn mul(self, rhs: FlatMono<T, Z>) -> Self::Output {
        let coeff = &self.coeff * rhs.coeff;
        let mut powers = self.powers;
        for (a, b) in powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
        Self::Output::new(coeff, powers)
    }
}

impl<'a, T, const Z: usize> Mul<&'a FlatMono<T, Z>> for FlatMono<T, Z>
where
    T: Mul<&'a T>,
{
    type Output = FlatMono<<T as Mul<&'a T>>::Output, Z>;

    fn mul(mut self, rhs: &'a FlatMono<T, Z>) -> Self::Output {
        let coeff = self.coeff * &rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
        Self::Output::new(coeff, self.powers)
    }
}

impl<'a, 'b, T, const Z: usize> Mul<&'b FlatMono<T, Z>> for &'a FlatMono<T, Z>
where
    &'a T: Mul<&'b T>,
{
    type Output = FlatMono<<&'a T as Mul<&'b T>>::Output, Z>;

    fn mul(self, rhs: &'b FlatMono<T, Z>) -> Self::Output {
        let coeff = &self.coeff * &rhs.coeff;
        let mut powers = self.powers;
        for (a, b) in powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
        Self::Output::new(coeff, powers)
    }
}

impl<T, const Z: usize> DivAssign<FlatMono<T, Z>> for FlatMono<T, Z>
where
    T: DivAssign<T>,
{
    fn div_assign(&mut self, rhs: FlatMono<T, Z>) {
        self.coeff /= rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
    }
}

impl<T, const Z: usize> Div<FlatMono<T, Z>> for FlatMono<T, Z>
where
    T: Div<T>,
{
    type Output = FlatMono<<T as Div<T>>::Output, Z>;

    fn div(mut self, rhs: FlatMono<T, Z>) -> Self::Output {
        let coeff = self.coeff / rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
        Self::Output::new(coeff, self.powers)
    }
}

impl<'a, T, const Z: usize> Div<FlatMono<T, Z>> for &'a FlatMono<T, Z>
where
    &'a T: Div<T>,
{
    type Output = FlatMono<<&'a T as Div<T>>::Output, Z>;

    fn div(self, rhs: FlatMono<T, Z>) -> Self::Output {
        let coeff = &self.coeff / rhs.coeff;
        let mut powers = self.powers;
        for (a, b) in powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
        Self::Output::new(coeff, powers)
    }
}

impl<'a, T, const Z: usize> Div<&'a FlatMono<T, Z>> for FlatMono<T, Z>
where
    T: Div<&'a T>,
{
    type Output = FlatMono<<T as Div<&'a T>>::Output, Z>;

    fn div(mut self, rhs: &'a FlatMono<T, Z>) -> Self::Output {
        let coeff = self.coeff / &rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
        Self::Output::new(coeff, self.powers)
    }
}

impl<'a, 'b, T, const Z: usize> Div<&'b FlatMono<T, Z>> for &'a FlatMono<T, Z>
where
    &'a T: Div<&'b T>,
{
    type Output = FlatMono<<&'a T as Div<&'b T>>::Output, Z>;

    fn div(self, rhs: &'b FlatMono<T, Z>) -> Self::Output {
        let coeff = &self.coeff / &rhs.coeff;
        let mut powers = self.powers;
        for (a, b) in powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
        Self::Output::new(coeff, powers)
    }
}

impl<T, const Z: usize> MulAssign<T> for FlatMono<T, Z>
where
    T: MulAssign,
{
    fn mul_assign(&mut self, rhs: T) {
        self.coeff *= rhs;
    }
}

impl<'a, T, const Z: usize> MulAssign<&'a T> for FlatMono<T, Z>
where
    T: MulAssign<&'a T>,
{
    fn mul_assign(&mut self, rhs: &'a T) {
        self.coeff *= rhs;
    }
}

impl<T, const Z: usize> Mul<T> for FlatMono<T, Z>
where
    T: Mul<T>,
{
    type Output = FlatMono<<T as Mul<T>>::Output, Z>;

    fn mul(self, rhs: T) -> Self::Output {
        Self::Output::new(self.coeff * rhs, self.powers)
    }
}

impl<'a, T, const Z: usize> Mul<&'a T> for FlatMono<T, Z>
where
    T: Mul<&'a T>,
{
    type Output = FlatMono<<T as Mul<&'a T>>::Output, Z>;

    fn mul(self, rhs: &'a T) -> Self::Output {
        Self::Output::new(self.coeff * rhs, self.powers)
    }
}

impl<'a, T, const Z: usize> Mul<T> for &'a FlatMono<T, Z>
where
    &'a T: Mul<T>,
{
    type Output = FlatMono<<&'a T as Mul<T>>::Output, Z>;

    fn mul(self, rhs: T) -> Self::Output {
        Self::Output::new((&self.coeff) * rhs, self.powers)
    }
}

impl<'a, 'b, T, const Z: usize> Mul<&'a T> for &'b FlatMono<T, Z>
where
    &'b T: Mul<&'a T>,
{
    type Output = FlatMono<<&'b T as Mul<&'a T>>::Output, Z>;

    fn mul(self, rhs: &'a T) -> Self::Output {
        Self::Output::new((&self.coeff) * rhs, self.powers)
    }
}

impl<T, const Z: usize> DivAssign<T> for FlatMono<T, Z>
where
    T: DivAssign,
{
    fn div_assign(&mut self, rhs: T) {
        self.coeff /= rhs;
    }
}

impl<'a, T, const Z: usize> DivAssign<&'a T> for FlatMono<T, Z>
where
    T: DivAssign<&'a T>,
{
    fn div_assign(&mut self, rhs: &'a T) {
        self.coeff /= rhs;
    }
}

impl<T, const Z: usize> Div<T> for FlatMono<T, Z>
where
    T: Div<T>,
{
    type Output = FlatMono<<T as Div<T>>::Output, Z>;

    fn div(self, rhs: T) -> Self::Output {
        Self::Output::new(self.coeff / rhs, self.powers)
    }
}

impl<'a, T, const Z: usize> Div<&'a T> for FlatMono<T, Z>
where
    T: Div<&'a T>,
{
    type Output = FlatMono<<T as Div<&'a T>>::Output, Z>;

    fn div(self, rhs: &'a T) -> Self::Output {
        Self::Output::new(self.coeff / rhs, self.powers)
    }
}

impl<'a, T, const Z: usize> Div<T> for &'a FlatMono<T, Z>
where
    &'a T: Div<T>,
{
    type Output = FlatMono<<&'a T as Div<T>>::Output, Z>;

    fn div(self, rhs: T) -> Self::Output {
        Self::Output::new((&self.coeff) / rhs, self.powers)
    }
}

impl<'a, 'b, T, const Z: usize> Div<&'a T> for &'b FlatMono<T, Z>
where
    &'b T: Div<&'a T>,
{
    type Output = FlatMono<<&'b T as Div<&'a T>>::Output, Z>;

    fn div(self, rhs: &'a T) -> Self::Output {
        Self::Output::new((&self.coeff) / rhs, self.powers)
    }
}

impl<T: Neg, const Z: usize> Neg for FlatMono<T, Z> {
    type Output = FlatMono<<T as Neg>::Output, Z>;

    fn neg(self) -> Self::Output {
        Self::Output::new(self.coeff.neg(), self.powers)
    }
}

impl<T: Zero, const Z: usize> Zero for FlatMono<T, Z> {
    fn zero() -> Self {
        Self {
            coeff: T::zero(),
            powers: [0; Z],
        }
    }

    fn is_zero(&self) -> bool {
        self.coeff.is_zero()
    }
}

impl<T: One, const Z: usize> One for FlatMono<T, Z> {
    fn one() -> Self {
        Self {
            coeff: T::one(),
            powers: [0; Z],
        }
    }

    fn is_one(&self) -> bool {
        self.coeff.is_one() && self.powers.iter().all(|&p| p == 0)
    }
}

impl<T: One, const Z: usize> From<T> for FlatMono<T, Z> {
    fn from(source: T) -> Self {
        Self::new(source, [0; Z])
    }
}

impl<const P: u64, const Z: usize> TryEval<[Z64<P>; Z]>
    for FlatMono<Z64<P>, Z>
{
    type Output = Z64<P>;

    fn try_eval(&self, x: &[Z64<P>; Z]) -> Option<Self::Output> {
        Some(
            x.iter()
                .zip(self.powers.iter().copied())
                .fold(self.coeff, |acc, (x, y)| acc * x.pow(y)),
        )
    }
}

impl<const P: u64, const Z: usize> Eval<[Z64<P>; Z]> for FlatMono<Z64<P>, Z> {}

impl<const P: u64, const Z: usize> TryEval<[Z64<P>; Z]>
    for FlatMono<Integer, Z>
{
    type Output = Z64<P>;

    fn try_eval(&self, x: &[Z64<P>; Z]) -> Option<Self::Output> {
        let coeff = unsafe { Z64::new_unchecked(self.coeff.mod_u64(P)) };
        Some(
            x.iter()
                .zip(self.powers.iter().copied())
                .fold(coeff, |acc, (x, y)| acc * x.pow(y)),
        )
    }
}

impl<const P: u64, const Z: usize> Eval<[Z64<P>; Z]> for FlatMono<Integer, Z> {}

impl<const Z: usize> TryEval<[Integer; Z]> for FlatMono<Integer, Z> {
    type Output = Integer;

    fn try_eval(&self, x: &[Integer; Z]) -> Option<Self::Output> {
        Some(
            &self.coeff
                * x.iter()
                    .zip(self.powers.iter().copied())
                    .fold(Integer::one(), |acc, (x, y)| {
                        acc * x.pow(y).complete()
                    }),
        )
    }
}

impl<const Z: usize> Eval<[Integer; Z]> for FlatMono<Integer, Z> {}

pub struct FmtFlatMono<'a, 'b, V, T, const Z: usize> {
    m: &'a FlatMono<T, Z>,
    vars: &'b [V],
}

impl<'a, 'b, V: Display, T: 'a, const Z: usize> WithVars<'a, &'b [V; Z]>
    for FlatMono<T, Z>
{
    type Output = FmtFlatMono<'a, 'b, V, T, Z>;

    fn with_vars(&'a self, vars: &'b [V; Z]) -> Self::Output {
        FmtFlatMono { m: self, vars }
    }
}

impl<'a, 'b, V: Display, const P: u64, const Z: usize> Display
    for FmtFlatMono<'a, 'b, V, Z64<P>, Z>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.m.coeff)?;
        if !self.m.is_one() && !self.m.is_zero() {
            debug_assert_eq!(self.m.powers.len(), self.vars.len());
            for (p, v) in self.m.powers.iter().zip(self.vars.iter()) {
                match p {
                    0 => {}
                    1 => write!(f, "*{v}")?,
                    _ => write!(f, "*{v}^{p}")?,
                }
            }
        }
        Ok(())
    }
}

impl<const P: u64, const Z: usize> Display for FlatMono<Z64<P>, Z> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut vars = [""; Z];
        vars.copy_from_slice(&ALL_VARS);
        self.with_vars(&vars).fmt(f)
    }
}

impl<'a, 'b, V: Display, const Z: usize> Display
    for FmtFlatMono<'a, 'b, V, Integer, Z>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        debug_assert_eq!(self.m.powers.len(), self.vars.len());
        let var_pows = self
            .vars
            .iter()
            .zip(self.m.powers.iter())
            .filter(|(_, p)| **p > 0);
        // omit coefficient if it's 1 and there are variables.
        // TODO: similar for -1
        if self.m.coeff.is_one() && self.m.powers.iter().any(|p| !p.is_zero()) {
            let mut var_pows = var_pows;
            let (v, p) = var_pows.next().unwrap();
            write!(f, "{v}")?;
            if *p != 1 {
                write!(f, "^{p}")?;
            }
            for (v, p) in var_pows {
                write!(f, "*{v}")?;
                if *p != 1 {
                    write!(f, "^{p}")?;
                }
            }
        } else {
            write!(f, "{}", self.m.coeff)?;
            if !self.m.is_one() && !self.m.is_zero() {
                for (v, p) in var_pows {
                    write!(f, "*{v}")?;
                    if *p != 1 {
                        write!(f, "^{p}")?;
                    }
                }
            }
        }
        Ok(())
    }
}

impl<const Z: usize> Display for FlatMono<Integer, Z> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut vars = [""; Z];
        vars.copy_from_slice(&ALL_VARS);
        self.with_vars(&vars).fmt(f)
    }
}

impl<const P: u64, const N: usize> From<FlatPoly<Z64<P>, N>>
    for FlatPoly<Integer, N>
{
    fn from(p: FlatPoly<Z64<P>, N>) -> Self {
        let terms = p
            .into_terms()
            .into_iter()
            .map(|c| FlatMono::new(i64::from(c.coeff).into(), c.powers))
            .collect();
        FlatPoly::from_raw_terms(terms)
    }
}

#[cfg(test)]
mod tests {
    use rand::Rng;
    use rand_xoshiro::rand_core::SeedableRng;

    use super::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn conv_1d() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_LEN: usize = 10;
        const MAX_POW: u32 = 10 * MAX_LEN as u32;
        const P: u64 = 61;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        for _ in 0..NTESTS {
            let len = rng.gen_range(0..=MAX_LEN);
            let terms: Vec<FlatMono<Z64<P>, 1>> =
                std::iter::repeat_with(|| {
                    FlatMono::new(rng.gen(), [rng.gen_range(0..=MAX_POW)])
                })
                .take(len)
                .collect();
            let poly = FlatPoly::from_terms(terms);
            eprintln!("original: {poly}");
            let as_dense = DensePoly::from(poly.clone());
            eprintln!("collected: {as_dense}");
            let reexpanded: FlatPoly<Z64<P>, 1> = as_dense.into();
            eprintln!("re-expanded: {reexpanded}");
            assert_eq!(poly, reexpanded);
        }
    }

    #[test]
    fn conv_2d() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_LEN: usize = 10;
        const MAX_POW: u32 = 10;
        const P: u64 = 61;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        for _ in 0..NTESTS {
            let len = rng.gen_range(0..=MAX_LEN);
            let terms: Vec<FlatMono<Z64<P>, 2>> =
                std::iter::repeat_with(|| {
                    FlatMono::new(
                        rng.gen(),
                        [
                            rng.gen_range(0..=MAX_POW),
                            rng.gen_range(0..=MAX_POW),
                        ],
                    )
                })
                .take(len)
                .collect();
            let poly = FlatPoly::from_terms(terms);
            eprintln!("original: {poly}");
            let as_dense = DensePoly2::from(poly.clone());
            eprintln!("collected: {as_dense}");
            let reexpanded: FlatPoly<Z64<P>, 2> = as_dense.into();
            eprintln!("re-expanded: {reexpanded}");
            assert_eq!(poly, reexpanded);
        }
    }
}
