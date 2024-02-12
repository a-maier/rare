use std::{
    fmt::{self, Display},
    ops::Div,
};

use ffnt::Z64;
use rug::{ops::NegAssign, Integer, Rational};
use thiserror::Error;

use crate::{
    arr::Arr,
    dense_poly::DensePoly,
    sparse_poly::{SparseMono, SparsePoly},
    traits::{One, TryEval, WithVars, Zero},
};

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Rat<N, D = N> {
    num: N,
    den: D,
}

impl<N, D> Rat<N, D> {
    pub fn from_num_den_unchecked(num: N, den: D) -> Self {
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
        Self {
            num: Zero::zero(),
            den: One::one(),
        }
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
        Self {
            num: One::one(),
            den: One::one(),
        }
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
    type Output =
        <<N as TryEval<U>>::Output as Div<<D as TryEval<U>>::Output>>::Output;

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

    fn with_vars(&'a self, vars: &'b [S; 1]) -> Self::Output {
        FmtUniRat::new(self, vars)
    }
}

impl<N, D> Display for Rat<N, D>
where
    N: Display + Zero,
    D: Display,
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
    D: Display + One + Zero,
{
    rat: &'a Rat<N, D>,
    var: &'b [V],
}

impl<'a, 'b, N, D, V: Display> FmtUniRat<'a, 'b, N, D, V>
where
    N: Display + One + Zero,
    D: Display + One + Zero,
{
    fn new(rat: &'a Rat<N, D>, var: &'b [V]) -> Self {
        Self { rat, var }
    }
}

impl<'a, 'b, V: Display, const P: u64> Display
    for FmtUniRat<'a, 'b, DensePoly<Z64<P>>, DensePoly<Z64<P>>, V>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.rat.is_zero() {
            return write!(f, "0");
        }
        let var = &[&self.var[0]];
        let rat = &self.rat;
        write!(
            f,
            "({})/({})",
            rat.num().with_vars(var),
            rat.den().with_vars(var)
        )
    }
}

#[derive(
    Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash, Error,
)]
pub struct NoneError {}

impl Display for NoneError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "No result found")
    }
}

pub type Rat64 = Rat<i64, u64>;

impl<const P: u64> TryFrom<Z64<P>> for Rat64 {
    type Error = NoneError;

    fn try_from(value: Z64<P>) -> Result<Self, Self::Error> {
        let max_bound = (P / 2) as i64;
        let max_den = (max_bound as f64).powf(1. / 4.) as i64;
        wang_reconstruct(value, max_bound / max_den, max_den)
            .ok_or(NoneError {})
    }
}

fn wang_reconstruct<const P: u64>(
    value: Z64<P>,
    max_num: i64,
    max_den: i64,
) -> Option<Rat64> {
    let mut v = Arr([P as i64, 0]);
    let mut w = Arr([i64::from(value), 1]);
    while w[0] > max_num {
        let q = v[0] / w[0];
        let z = v - w * q;
        (v, w) = (w, z);
    }
    if w[1] < 0 {
        w = -w;
    }
    if w[1] < max_den && gcd(w[0], w[1]) == 1 {
        Some(Rat64::from_num_den_unchecked(w[0], w[1] as u64))
    } else {
        None
    }
}

fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        (a, b) = (b, a.rem_euclid(b))
    }
    a
}

impl<const P: u64, const N: usize> From<Rat<SparsePoly<Z64<P>, N>>>
    for Rat<SparsePoly<Integer, N>>
{
    fn from(rat: Rat<SparsePoly<Z64<P>, N>>) -> Self {
        let (num, den) = rat.into_num_den();
        Rat::from_num_den_unchecked(num.into(), den.into())
    }
}

impl<const N: usize> From<Rat<SparsePoly<Rational, N>>>
    for Rat<SparsePoly<Integer, N>>
{
    fn from(rat: Rat<SparsePoly<Rational, N>>) -> Self {
        let (num, den) = rat.into_num_den();

        let mut lcm = Integer::one();
        let num_terms = num.terms().iter().map(|t| &t.coeff);
        let den_terms = den.terms().iter().map(|t| &t.coeff);

        for c in num_terms.chain(den_terms) {
            lcm.lcm_mut(c.denom());
        }
        // ensure unique normalisation: first coefficient in denominator is positive
        if den.terms()[0].coeff.numer().is_negative() {
            lcm.neg_assign()
        }

        let to_int = |p: SparsePoly<Rational, N>, lcm| {
            let terms = p
                .into_terms()
                .into_iter()
                .map(|t| {
                    let coeff: Rational = t.coeff * lcm;
                    debug_assert!(coeff.denom().is_one());
                    SparseMono::new(coeff.into_numer_denom().0, t.powers)
                })
                .collect();
            SparsePoly::from_raw_terms(terms)
        };

        let num = to_int(num, &lcm);
        let den = to_int(den, &lcm);
        Rat::from_num_den_unchecked(num, den)
    }
}

#[cfg(test)]
mod tests {
    use log::debug;
    use rand::Rng;
    use rand_xoshiro::rand_core::SeedableRng;

    use super::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn rec_rat_num64() {
        log_init();

        const NTESTS: u32 = 10000;
        const P: u64 = 1152921504606846883;

        let max_num = 1000;
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

        for _ in 0..NTESTS {
            let mut num = rng.gen_range(-max_num..=max_num);
            let mut den = rng.gen_range(1..=max_num);
            let gcd = gcd(num, den);
            num /= gcd;
            den /= gcd;
            if num == 0 {
                den = 1;
            }
            let n: Z64<P> = Z64::new(num) / Z64::new(den);
            let rec = Rat64::try_from(n).unwrap();
            debug!("{num}/{den} == {rec}");
            assert_eq!(rec.num(), &num);
            assert_eq!(*rec.den(), den as u64);
        }
    }
}
