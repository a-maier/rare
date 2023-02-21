use std::{ops::{Add, Neg, Sub, MulAssign, DivAssign, Mul, Div, AddAssign, SubAssign}, fmt::{Display, self}, cmp::Ordering};

use galois_fields::Z64;

use crate::{traits::{Zero, One, WithVars}, dense_poly::ALL_VARS};

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct SparsePoly<T, const Z: usize> {
    terms: Vec<SparseMono<T, Z>>
}

impl<T, const Z: usize> SparsePoly<T, Z> {
    pub fn new() -> Self {
        Self{ terms: Vec::new() }
    }

    pub fn terms(&self) -> &[SparseMono<T, Z>] {
        self.terms.as_ref()
    }

    pub fn len(&self) -> usize {
        self.terms.len()
    }

    pub fn into_terms(self) -> Vec<SparseMono<T, Z>> {
        self.terms
    }
}

impl<T: Zero, const Z: usize> SparsePoly<T, Z> {
    pub fn from_raw_terms(terms: Vec<SparseMono<T, Z>>) -> Self {
        debug_assert!(! terms.iter().any(|c| c.is_zero()));
        Self { terms }
    }
}

impl<T: AddAssign + Zero, const Z: usize> SparsePoly<T, Z> {
    pub fn from_terms(mut terms: Vec<SparseMono<T, Z>>) -> Self {
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
                    res_terms.pop();
                }
            }
            if res_terms[res_terms.len() - 1].is_zero() {
                res_terms.pop();
            }
        }
        Self::from_raw_terms(res_terms)
    }
}

impl<T: Zero, const Z: usize> From<SparseMono<T, Z>> for SparsePoly<T, Z> {
    fn from(source: SparseMono<T, Z>) -> Self {
        if source.is_zero() {
            Self::zero()
        } else {
            Self::from_raw_terms(vec![source])
        }
    }
}

impl<T, const Z: usize> AddAssign for SparsePoly<T, Z>
where SparsePoly<T, Z>: Add<Output = Self>
{
    fn add_assign(&mut self, rhs: SparsePoly<T, Z>) {
        *self = std::mem::replace(self, Self::new()) + rhs;
    }
}

impl<T, const Z: usize> Add for SparsePoly<T, Z>
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
                },
                (_, None) => {
                    res.extend(lhs_terms);
                    break;
                },
                (Some(lhs), Some(rhs)) => {
                    match lhs.powers.cmp(&rhs.powers) {
                        Ordering::Less => res.push(lhs_terms.next().unwrap()),
                        Ordering::Greater => res.push(rhs_terms.next().unwrap()),
                        Ordering::Equal => {
                            let mut lhs = lhs_terms.next().unwrap();
                            let rhs = rhs_terms.next().unwrap();
                            lhs.coeff += rhs.coeff;
                            res.push(lhs)
                        }
                    }
                }
            }
        }
        Self::from_raw_terms(res)
    }
}

impl<T, const Z: usize> SubAssign for SparsePoly<T, Z>
where
    SparsePoly<T, Z>: AddAssign,
    for<'a> &'a T: Neg<Output = T>
{
    fn sub_assign(&mut self, mut rhs: SparsePoly<T, Z>) {
        for term in &mut rhs.terms {
            term.coeff = term.coeff.neg();
        }
        *self += rhs;
    }
}

impl<T, const Z: usize> Sub for SparsePoly<T, Z>
where SparsePoly<T, Z>: SubAssign
{
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<'a, 'b, T, const Z: usize> Mul<&'b SparsePoly<T, Z>> for &'a SparsePoly<T, Z>
where
    T: Zero + AddAssign,
&'a T: Mul<&'b T, Output = T>
{
    type Output = SparsePoly<T, Z>;

    fn mul(self, rhs: &'b SparsePoly<T, Z>) -> Self::Output {
        let mut res = Vec::with_capacity(self.len() * rhs.len());
        for lhs in self.terms() {
            for rhs in rhs.terms() {
                res.push(lhs * rhs)
            }
        }
        Self::Output::from_terms(res)
    }
}

impl<T, const Z: usize> MulAssign<T> for SparsePoly<T, Z>
where T: Copy + MulAssign + Zero
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

impl<'a, T, const Z: usize> MulAssign<&'a T> for SparsePoly<T, Z>
where T: MulAssign<&'a T> + Zero
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

impl<T, const Z: usize> Mul<T> for SparsePoly<T, Z>
where
    T: Copy + Mul<T> + Zero,
    <T as Mul<T>>::Output: Zero
{
    type Output = SparsePoly<<T as Mul<T>>::Output, Z>;

    fn mul(self, rhs: T) -> Self::Output {
        if rhs.is_zero() {
            Zero::zero()
        } else {
            Self::Output::from_raw_terms(
                self.into_terms()
                    .into_iter()
                    .map(|t| t * rhs)
                    .collect()
            )
        }
    }
}

impl<'a, T, const Z: usize> Mul<&'a T> for SparsePoly<T, Z>
where
    T: Mul<&'a T> + Zero,
    <T as Mul<&'a T>>::Output: Zero
{
    type Output = SparsePoly<<T as Mul<&'a T>>::Output, Z>;

    fn mul(self, rhs: &'a T) -> Self::Output {
        if rhs.is_zero() {
            Zero::zero()
        } else {
            Self::Output::from_raw_terms(
                self.into_terms()
                    .into_iter()
                    .map(|t| t * rhs)
                    .collect()
            )
        }
    }
}

impl<T, const Z: usize> DivAssign<T> for SparsePoly<T, Z>
where T: Copy + DivAssign + Zero
{
    fn div_assign(&mut self, rhs: T) {
        for term in &mut self.terms {
            *term /= rhs;
            debug_assert!(!term.is_zero())
        }
    }
}

impl<'a, T, const Z: usize> DivAssign<&'a T> for SparsePoly<T, Z>
where T: DivAssign<&'a T> + Zero
{
    fn div_assign(&mut self, rhs: &'a T) {
        for term in &mut self.terms {
            *term /= rhs;
            debug_assert!(!term.is_zero())
        }
    }
}

impl<T, const Z: usize> Div<T> for SparsePoly<T, Z>
where
    T: Copy + Div<T>,
    <T as Div<T>>::Output: Zero
{
    type Output = SparsePoly<<T as Div<T>>::Output, Z>;

    fn div(self, rhs: T) -> Self::Output {
        Self::Output::from_raw_terms(
            self.into_terms()
                .into_iter()
                .map(|t| t / rhs)
                .collect()
        )
    }
}

impl<'a, T, const Z: usize> Div<&'a T> for SparsePoly<T, Z>
where
    T: Div<&'a T>,
    <T as Div<&'a T>>::Output: Zero
{
    type Output = SparsePoly<<T as Div<&'a T>>::Output, Z>;

    fn div(self, rhs: &'a T) -> Self::Output {
        Self::Output::from_raw_terms(
            self.into_terms()
                .into_iter()
                .map(|t| t / rhs)
                .collect()
        )
    }
}

impl<T, const Z: usize> Zero for SparsePoly<T, Z> {
    fn zero() -> Self {
        Self::new()
    }

    fn is_zero(&self) -> bool {
        self.terms().is_empty()
    }
}

impl<T: Zero + One, const Z: usize> One for SparsePoly<T, Z> {
    fn one() -> Self {
        Self::from_raw_terms(vec![SparseMono::<T, Z>::one()])
    }

    fn is_one(&self) -> bool {
        if let Some(first) = self.terms().first() {
            self.len() == 1 && first.is_one()
        } else {
            false
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct SparseMono<T, const Z: usize> {
    pub powers: [u32; Z],
    pub coeff: T,
}

impl<T, const Z: usize> SparseMono<T, Z> {
    pub fn new(coeff: T, powers: [u32; Z]) -> Self {
        Self { coeff, powers }
    }
}

impl<T: Default, const Z: usize> Default for SparseMono<T, Z> {
    fn default() -> Self {
        Self {
            coeff: Default::default(),
            powers: [Default::default(); Z]
        }
    }
}

impl<T: AddAssign + Zero, const Z: usize> Add for SparseMono<T, Z> {
    type Output = SparsePoly<T, Z>;

    fn add(self, rhs: Self) -> Self::Output {
        Self::Output::from_terms(vec![self, rhs])
    }
}

impl<T: AddAssign + Zero, const Z: usize> Sub for SparseMono<T, Z>
where
    SparseMono<T, Z>: Neg<Output = Self>
{
    type Output = SparsePoly<T, Z>;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output::from_terms(vec![self, rhs.neg()])
    }
}

impl<T, const Z: usize> MulAssign for SparseMono<T, Z>
where T: MulAssign
{
    fn mul_assign(&mut self, rhs: SparseMono<T, Z>) {
        self.coeff *= rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
    }
}

impl<'a, T, const Z: usize> MulAssign<&'a SparseMono<T,Z>> for SparseMono<T, Z>
where T: MulAssign<&'a T>
{
    fn mul_assign(&mut self, rhs: &'a SparseMono<T, Z>) {
        self.coeff *= &rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.iter()) {
            *a += *b;
        }
    }
}

impl<T, const Z: usize> Mul<SparseMono<T, Z>> for SparseMono<T, Z>
where T: Mul
{
    type Output = SparseMono<<T as Mul>::Output, Z>;

    fn mul(mut self, rhs: SparseMono<T, Z>) -> Self::Output {
        let coeff = self.coeff * rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
        Self::Output::new(coeff, self.powers)
    }
}

impl<'a, T, const Z: usize> Mul<SparseMono<T, Z>> for &'a SparseMono<T, Z>
where &'a T: Mul<T>
{
    type Output = SparseMono<<&'a T as Mul<T>>::Output, Z>;

    fn mul(self, rhs: SparseMono<T, Z>) -> Self::Output {
        let coeff = &self.coeff * rhs.coeff;
        let mut powers = self.powers;
        for (a, b) in powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
        Self::Output::new(coeff, powers)
    }
}

impl<'a, T, const Z: usize> Mul<&'a SparseMono<T, Z>> for SparseMono<T, Z>
where T: Mul<&'a T>
{
    type Output = SparseMono<<T as Mul<&'a T>>::Output, Z>;

    fn mul(mut self, rhs: &'a SparseMono<T, Z>) -> Self::Output {
        let coeff = self.coeff * &rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
        Self::Output::new(coeff, self.powers)
    }
}

impl<'a, 'b, T, const Z: usize> Mul<&'b SparseMono<T, Z>> for &'a SparseMono<T, Z>
where &'a T: Mul<&'b T>
{
    type Output = SparseMono<<&'a T as Mul<&'b T>>::Output, Z>;

    fn mul(self, rhs: &'b SparseMono<T, Z>) -> Self::Output {
        let coeff = &self.coeff * &rhs.coeff;
        let mut powers = self.powers;
        for (a, b) in powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
        Self::Output::new(coeff, powers)
    }
}

impl<T, const Z: usize> DivAssign<SparseMono<T, Z>> for SparseMono<T, Z>
where T: DivAssign<T>
{
    fn div_assign(&mut self, rhs: SparseMono<T, Z>) {
        self.coeff /= rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
    }
}

impl<T, const Z: usize> Div<SparseMono<T, Z>> for SparseMono<T, Z>
where T: Div<T>
{
    type Output = SparseMono<<T as Div<T>>::Output, Z>;

    fn div(mut self, rhs: SparseMono<T, Z>) -> Self::Output {
        let coeff = self.coeff / rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
        Self::Output::new(coeff, self.powers)
    }
}

impl<'a, T, const Z: usize> Div<SparseMono<T, Z>> for &'a SparseMono<T, Z>
where &'a T: Div<T>
{
    type Output = SparseMono<<&'a T as Div<T>>::Output, Z>;

    fn div(self, rhs: SparseMono<T, Z>) -> Self::Output {
        let coeff = &self.coeff / rhs.coeff;
        let mut powers = self.powers;
        for (a, b) in powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
        Self::Output::new(coeff, powers)
    }
}

impl<'a, T, const Z: usize> Div<&'a SparseMono<T, Z>> for SparseMono<T, Z>
where T: Div<&'a T>
{
    type Output = SparseMono<<T as Div<&'a T>>::Output, Z>;

    fn div(mut self, rhs: &'a SparseMono<T, Z>) -> Self::Output {
        let coeff = self.coeff / &rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
        Self::Output::new(coeff, self.powers)
    }
}

impl<'a, 'b, T, const Z: usize> Div<&'b SparseMono<T, Z>> for &'a SparseMono<T, Z>
where &'a T: Div<&'b T>
{
    type Output = SparseMono<<&'a T as Div<&'b T>>::Output, Z>;

    fn div(self, rhs: &'b SparseMono<T, Z>) -> Self::Output {
        let coeff = &self.coeff / &rhs.coeff;
        let mut powers = self.powers;
        for (a, b) in powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
        Self::Output::new(coeff, powers)
    }
}

impl<T, const Z: usize> MulAssign<T> for SparseMono<T, Z>
where T: MulAssign
{
    fn mul_assign(&mut self, rhs: T) {
        self.coeff *= rhs;
    }
}

impl<'a, T, const Z: usize> MulAssign<&'a T> for SparseMono<T, Z>
where T: MulAssign<&'a T>
{
    fn mul_assign(&mut self, rhs: &'a T) {
        self.coeff *= rhs;
    }
}

impl<T, const Z: usize> Mul<T> for SparseMono<T, Z>
where T: Mul<T>
{
    type Output = SparseMono<<T as Mul<T>>::Output, Z>;

    fn mul(self, rhs: T) -> Self::Output {
        Self::Output::new(self.coeff * rhs, self.powers)
    }
}

impl<'a, T, const Z: usize> Mul<&'a T> for SparseMono<T, Z>
where T: Mul<&'a T>
{
    type Output = SparseMono<<T as Mul<&'a T>>::Output, Z>;

    fn mul(self, rhs: &'a T) -> Self::Output {
        Self::Output::new(self.coeff * rhs, self.powers)
    }
}

impl<'a, T, const Z: usize> Mul<T> for &'a SparseMono<T, Z>
where &'a T: Mul<T>
{
    type Output = SparseMono<<&'a T as Mul<T>>::Output, Z>;

    fn mul(self, rhs: T) -> Self::Output {
        Self::Output::new((&self.coeff) * rhs, self.powers)
    }
}

impl<'a, 'b, T, const Z: usize> Mul<&'a T> for &'b SparseMono<T, Z>
where &'b T: Mul<&'a T>
{
    type Output = SparseMono<<&'b T as Mul<&'a T>>::Output, Z>;

    fn mul(self, rhs: &'a T) -> Self::Output {
        Self::Output::new((&self.coeff) * rhs, self.powers)
    }
}

impl<T, const Z: usize> DivAssign<T> for SparseMono<T, Z>
where T: DivAssign
{
    fn div_assign(&mut self, rhs: T) {
        self.coeff /= rhs;
    }
}

impl<'a, T, const Z: usize> DivAssign<&'a T> for SparseMono<T, Z>
where T: DivAssign<&'a T>
{
    fn div_assign(&mut self, rhs: &'a T) {
        self.coeff /= rhs;
    }
}

impl<T, const Z: usize> Div<T> for SparseMono<T, Z>
where T: Div<T>
{
    type Output = SparseMono<<T as Div<T>>::Output, Z>;

    fn div(self, rhs: T) -> Self::Output {
        Self::Output::new(self.coeff / rhs, self.powers)
    }
}

impl<'a, T, const Z: usize> Div<&'a T> for SparseMono<T, Z>
where T: Div<&'a T>
{
    type Output = SparseMono<<T as Div<&'a T>>::Output, Z>;

    fn div(self, rhs: &'a T) -> Self::Output {
        Self::Output::new(self.coeff / rhs, self.powers)
    }
}

impl<'a, T, const Z: usize> Div<T> for &'a SparseMono<T, Z>
where &'a T: Div<T>
{
    type Output = SparseMono<<&'a T as Div<T>>::Output, Z>;

    fn div(self, rhs: T) -> Self::Output {
        Self::Output::new((&self.coeff) / rhs, self.powers)
    }
}

impl<'a, 'b, T, const Z: usize> Div<&'a T> for &'b SparseMono<T, Z>
where &'b T: Div<&'a T>
{
    type Output = SparseMono<<&'b T as Div<&'a T>>::Output, Z>;

    fn div(self, rhs: &'a T) -> Self::Output {
        Self::Output::new((&self.coeff) / rhs, self.powers)
    }
}

impl<T: Neg, const Z: usize> Neg for SparseMono<T, Z> {
    type Output = SparseMono<<T as Neg>::Output, Z>;

    fn neg(self) -> Self::Output {
        Self::Output::new(self.coeff.neg(), self.powers)
    }
}

impl<T: Zero, const Z: usize> Zero for SparseMono<T, Z> {
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

impl<T: One, const Z: usize> One for SparseMono<T, Z> {
    fn one() -> Self {
        Self {
            coeff: T::one(),
            powers: [0; Z],
        }
    }

    fn is_one(&self) -> bool {
        self.coeff.is_one()
            && self.powers.iter().all(|&p| p == 0)
    }
}

impl<T: One, const Z: usize> From<T> for SparseMono<T, Z> {
    fn from(source: T) -> Self {
        Self::new(source, [0; Z])
    }
}

pub struct FmtSMonomial<'a, V, const P: u64, const Z: usize> {
    m: SparseMono<Z64<P>, Z>,
    vars: &'a [V],
}

impl<'a, 'b, V: Display, const P: u64, const Z: usize> WithVars<'a, &'b [V; Z]> for SparseMono<Z64<P>, Z> {
    type Output = FmtSMonomial<'b, V, P, Z>;

    fn with_vars(&'a self, vars: &'b [V; Z]) -> Self::Output {
        FmtSMonomial{ m: *self, vars }
    }
}


impl<'a, V: Display, const P: u64, const Z: usize> Display for FmtSMonomial<'a, V, P, Z> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.m.coeff)?;
        if !self.m.is_one() && !self.m.is_zero() {
            debug_assert_eq!(self.m.powers.len(), self.vars.len());
            for (p, v) in self.m.powers.iter().zip(self.vars.iter()) {
                match p {
                    0 => {},
                    1 => write!(f, "*{v}")?,
                    _ => write!(f, "*{v}^{p}")?,
                }
            }
        }
        Ok(())
    }
}

impl<const P: u64, const Z: usize> Display for SparseMono<Z64<P>, Z> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut vars = [""; Z];
        vars.copy_from_slice(&ALL_VARS);
        self.with_vars(&vars).fmt(f)
    }
}
