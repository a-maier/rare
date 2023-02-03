use std::{fmt::{Display, self}, ops::Div};

use galois_fields::Z64;

use crate::{traits::{Zero, One, WithVars, TryEval}, dense_poly::DensePoly};

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct DenseRat<T> {
    num: DensePoly<T>,
    den: DensePoly<T>
}

impl<T> DenseRat<T> {
    pub fn from_num_den_unchecked(
        num: DensePoly<T>,
        den: DensePoly<T>
    ) -> Self {
        debug_assert!(!den.is_zero());
        Self { num, den }
    }

    pub fn num(&self) -> &DensePoly<T> {
        &self.num
    }

    pub fn den(&self) -> &DensePoly<T> {
        &self.den
    }
}

impl<T: One> DenseRat<T> {
    pub fn new() -> Self {
        Self { num: Zero::zero(), den: One::one() }
    }
}

impl<T: One> Zero for DenseRat<T> {
    fn zero() -> Self {
        Self::new()
    }

    fn is_zero(&self) -> bool {
        self.num().is_zero()
    }
}

impl<T: One> One for DenseRat<T> {
    fn one() -> Self {
        Self{ num: One::one(), den: One::one() }
    }

    fn is_one(&self) -> bool {
        self.num().is_one() && self.den().is_one()
    }
}

impl<T, U> TryEval<U> for DenseRat<T>
where
    DensePoly<T>: TryEval<U>,
    <DensePoly<T> as TryEval<U>>::Output: Div + Zero
{
    type Output = <<DensePoly<T> as TryEval<U>>::Output as Div>::Output;

    fn try_eval(&self, pt: &U) -> Option<Self::Output> {
        let den = self.den().try_eval(pt)?;
        // TODO: might be better to check if division suceeds
        //       but TryDiv is not a frequently implemented trait
        //       and we would always have to calculate num
        if den.is_zero() {
            return None;
        }
        let num = self.num().try_eval(pt)?;
        Some(num / den)
    }
}


impl<'a, 'b, T, S: Display> WithVars<'a, &'b [S; 1]> for DenseRat<T>
where
    T: Display + One + Zero + 'a,
{
    type Output = FmtUniRat<'a, 'b, T, S>;

    fn with_vars(&'a self, vars: &'b[S; 1]) -> Self::Output {
        FmtUniRat::new(self, vars)
    }
}

impl<T> Display for DenseRat<T>
where DensePoly<T>: Display
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.num.is_zero() {
            write!(f, "0")
        } else {
            write!(f, "({}) / ({})", self.num(), self.den())
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct FmtUniRat<'a, 'b, T: Display + One + Zero, V: Display> {
    rat: &'a DenseRat<T>,
    var: &'b[V],
}

impl<'a, 'b, T: Display + One + Zero, V: Display> FmtUniRat<'a, 'b, T, V> {
    fn new(rat: &'a DenseRat<T>, var: &'b [V]) -> Self {
        Self { rat, var }
    }
}

impl<'a, 'b, V: Display, const P: u64> Display for FmtUniRat<'a, 'b, Z64<P>, V> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.rat.is_zero() {
            return write!(f, "0");
        }
        let var = &[&self.var[0]];
        let rat = &self.rat;
        write!(f, "({})/({})", rat.num().with_vars(var), rat.den().with_vars(var))
    }
}
