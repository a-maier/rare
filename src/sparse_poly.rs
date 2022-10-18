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

impl<T, const Z: usize> Zero for SparsePoly<T, Z> {
    fn zero() -> Self {
        Self::new()
    }

    fn is_zero(&self) -> bool {
        self.terms().is_empty()
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

impl<T, U, const Z: usize> MulAssign<SparseMono<U, Z>> for SparseMono<T, Z>
where T: MulAssign<U>
{
    fn mul_assign(&mut self, rhs: SparseMono<U, Z>) {
        self.coeff *= rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
    }
}

impl<T, U, const Z: usize> Mul<SparseMono<U, Z>> for SparseMono<T, Z>
where T: Mul<U>
{
    type Output = SparseMono<<T as Mul<U>>::Output, Z>;

    fn mul(mut self, rhs: SparseMono<U, Z>) -> Self::Output {
        let coeff = self.coeff * rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
        Self::Output::new(coeff, self.powers)
    }
}

impl<'a, T, U, const Z: usize> Mul<SparseMono<U, Z>> for &'a SparseMono<T, Z>
where &'a T: Mul<U>
{
    type Output = SparseMono<<&'a T as Mul<U>>::Output, Z>;

    fn mul(self, rhs: SparseMono<U, Z>) -> Self::Output {
        let coeff = &self.coeff * rhs.coeff;
        let mut powers = self.powers;
        for (a, b) in powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
        Self::Output::new(coeff, powers)
    }
}

impl<'a, T, U, const Z: usize> Mul<&'a SparseMono<U, Z>> for SparseMono<T, Z>
where T: Mul<&'a U>
{
    type Output = SparseMono<<T as Mul<&'a U>>::Output, Z>;

    fn mul(mut self, rhs: &'a SparseMono<U, Z>) -> Self::Output {
        let coeff = self.coeff * &rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
        Self::Output::new(coeff, self.powers)
    }
}

impl<'a, 'b, T, U, const Z: usize> Mul<&'b SparseMono<U, Z>> for &'a SparseMono<T, Z>
where &'a T: Mul<&'b U>
{
    type Output = SparseMono<<&'a T as Mul<&'b U>>::Output, Z>;

    fn mul(self, rhs: &'b SparseMono<U, Z>) -> Self::Output {
        let coeff = &self.coeff * &rhs.coeff;
        let mut powers = self.powers;
        for (a, b) in powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a += b;
        }
        Self::Output::new(coeff, powers)
    }
}

impl<T, U, const Z: usize> DivAssign<SparseMono<U, Z>> for SparseMono<T, Z>
where T: DivAssign<U>
{
    fn div_assign(&mut self, rhs: SparseMono<U, Z>) {
        self.coeff /= rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
    }
}

impl<T, U, const Z: usize> Div<SparseMono<U, Z>> for SparseMono<T, Z>
where T: Div<U>
{
    type Output = SparseMono<<T as Div<U>>::Output, Z>;

    fn div(mut self, rhs: SparseMono<U, Z>) -> Self::Output {
        let coeff = self.coeff / rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
        Self::Output::new(coeff, self.powers)
    }
}

impl<'a, T, U, const Z: usize> Div<SparseMono<U, Z>> for &'a SparseMono<T, Z>
where &'a T: Div<U>
{
    type Output = SparseMono<<&'a T as Div<U>>::Output, Z>;

    fn div(self, rhs: SparseMono<U, Z>) -> Self::Output {
        let coeff = &self.coeff / rhs.coeff;
        let mut powers = self.powers;
        for (a, b) in powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
        Self::Output::new(coeff, powers)
    }
}

impl<'a, T, U, const Z: usize> Div<&'a SparseMono<U, Z>> for SparseMono<T, Z>
where T: Div<&'a U>
{
    type Output = SparseMono<<T as Div<&'a U>>::Output, Z>;

    fn div(mut self, rhs: &'a SparseMono<U, Z>) -> Self::Output {
        let coeff = self.coeff / &rhs.coeff;
        for (a, b) in self.powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
        Self::Output::new(coeff, self.powers)
    }
}

impl<'a, 'b, T, U, const Z: usize> Div<&'b SparseMono<U, Z>> for &'a SparseMono<T, Z>
where &'a T: Div<&'b U>
{
    type Output = SparseMono<<&'a T as Div<&'b U>>::Output, Z>;

    fn div(self, rhs: &'b SparseMono<U, Z>) -> Self::Output {
        let coeff = &self.coeff / &rhs.coeff;
        let mut powers = self.powers;
        for (a, b) in powers.iter_mut().zip(rhs.powers.into_iter()) {
            *a -= b;
        }
        Self::Output::new(coeff, powers)
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
