use std::{ops::{AddAssign, SubAssign, Add, Sub, Mul, MulAssign, DivAssign, Div}, fmt::{Display, self}, cmp::min};

use num_traits::{Zero, One};

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

    pub fn degree(&self) -> Option<usize> {
        if self.coeff.is_empty() {
            None
        } else {
            Some(self.coeff.len() - 1)
        }
    }

    pub fn eval<'a, 'b, X>(&'a self, x: &'b X) -> X
    where
        X: Zero + Mul<&'b X, Output = X> + Add<&'a T, Output = X>
    {
        self.coeff.iter().rev().fold(
            X::zero(),
            |acc, c| acc * x + c
        )
    }

    fn delete_trailing_zeroes(&mut self) {
        let last_nonzero = self.coeff.iter()
            .rposition(|c| !c.is_zero())
            .map(|pos| pos + 1)
            .unwrap_or_default();
        self.coeff.truncate(last_nonzero);
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

impl<T: AddAssign + Zero> Zero for UniPolynomial<T> {
    fn zero() -> Self {
        Self::new()
    }

    fn is_zero(&self) -> bool {
        self.coeff.is_empty()
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct FmtUniPoly<'a, T: Display + One + Zero, V: Display> {
    poly: &'a UniPolynomial<T>,
    var: V,
}

impl<'a, T: Display + One + Zero, V: Display> FmtUniPoly<'a, T, V> {
    fn new(poly: &'a UniPolynomial<T>, var: V) -> Self {
        Self { poly, var }
    }
}

impl<T: Display + One + Zero + Eq> Display for UniPolynomial<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        FmtUniPoly::new(self, "x").fmt(f)
    }
}

impl<'a, T: Display + One + Zero + Eq, V: Display> Display for FmtUniPoly<'a, T, V> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.poly.coeff.is_empty() {
            return write!(f, "0");
        }
        let mut coeff = self.poly.coeff.iter().enumerate().filter(|(_, c)| !c.is_zero());
        let var = &self.var;
        if let Some((pow, first_coeff)) = coeff.next() {
            match pow {
                0 => write!(f, "{first_coeff}"),
                1 => if first_coeff.is_one() {
                    write!(f, "{var}")
                } else {
                    write!(f, "{first_coeff}*{var}")
                },
                _ => if first_coeff.is_one() {
                    write!(f, "{var}^{pow}")
                } else {
                    write!(f, "{first_coeff}*{var}^{pow}")
                }
            }?;
        }
        for (pow, coeff) in coeff {
            let coeff_str = coeff.to_string();
            let coeff_str = coeff_str.trim();
            let coeff_str = if let Some(coeff_str) = coeff_str.strip_prefix('-') {
                write!(f, " - ")?;
                coeff_str.trim_start()
            } else {
                write!(f, " + ")?;
                coeff_str
            };
            // TODO: code duplication
            match pow {
                0 => write!(f, "{coeff_str}"),
                1 => if coeff_str == "1" {
                    write!(f, "{var}")
                } else {
                    write!(f, "{coeff_str}*{var}")
                },
                _ => if coeff_str == "1" {
                    write!(f, "{var}^{pow}")
                } else {
                    write!(f, "{coeff_str}*{var}^{pow}")
                }
            }?;
        }
        Ok(())
    }
}
