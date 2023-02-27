use std::{fmt::{Display, self}, ops::Div};

use galois_fields::Z64;

use crate::{traits::{Zero, One, WithVars, TryEval}, dense_poly::DensePoly};

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Rat<N, D = N> {
    num: N,
    den: D
}

impl<N, D> Rat<N, D> {
    pub fn from_num_den_unchecked(
        num: N,
        den: D
    ) -> Self {
        //debug_assert!(!den.is_zero());
        Self { num, den }
    }

    pub fn num(&self) -> &N {
        &self.num
    }

    pub fn den(&self) -> &D {
        &self.den
    }

    pub fn into_num_den(self) -> (N, D) {
        (self.num, self.den)
    }
}

impl<N: Zero, D: One> Rat<N, D> {
    pub fn new() -> Self {
        Self { num: Zero::zero(), den: One::one() }
    }
}

impl<N: Zero, D: One> Zero for Rat<N, D> {
    fn zero() -> Self {
        Self::new()
    }

    fn is_zero(&self) -> bool {
        self.num().is_zero()
    }
}

impl<N: One, D: One> One for Rat<N, D> {
    fn one() -> Self {
        Self{ num: One::one(), den: One::one() }
    }

    fn is_one(&self) -> bool {
        self.num().is_one() && self.den().is_one()
    }
}

impl<N, D, U> TryEval<U> for Rat<N, D>
where
    N: TryEval<U>,
    D: TryEval<U>,
   <D as TryEval<U>>::Output: Zero,
   <N as TryEval<U>>::Output: Div<<D as TryEval<U>>::Output>,
{
    type Output = <<N as TryEval<U>>::Output as Div<<D as TryEval<U>>::Output>>::Output;

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


impl<'a, 'b, N, D, S: Display> WithVars<'a, &'b [S; 1]> for Rat<N, D>
where
    N: Display + One + Zero + 'a,
    D: Display + One + Zero + 'a,
{
    type Output = FmtUniRat<'a, 'b, N, D, S>;

    fn with_vars(&'a self, vars: &'b[S; 1]) -> Self::Output {
        FmtUniRat::new(self, vars)
    }
}

impl<N, D> Display for Rat<N, D>
where
    N: Display + Zero,
    D: Display
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
pub struct FmtUniRat<'a, 'b, N, D, V: Display>
where
    N: Display + One + Zero,
    D: Display + One + Zero
{
    rat: &'a Rat<N, D>,
    var: &'b[V],
}

impl<'a, 'b, N, D, V: Display> FmtUniRat<'a, 'b, N, D, V>
where
    N: Display + One + Zero,
    D: Display + One + Zero
{
    fn new(rat: &'a Rat<N, D>, var: &'b [V]) -> Self {
        Self { rat, var }
    }
}

impl<'a, 'b, V: Display, const P: u64> Display for FmtUniRat<'a, 'b, DensePoly<Z64<P>>, DensePoly<Z64<P>>, V> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.rat.is_zero() {
            return write!(f, "0");
        }
        let var = &[&self.var[0]];
        let rat = &self.rat;
        write!(f, "({})/({})", rat.num().with_vars(var), rat.den().with_vars(var))
    }
}
