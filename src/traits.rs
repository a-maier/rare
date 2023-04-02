use std::num::NonZeroUsize;

use galois_fields::Z64;
use rand::{thread_rng, Rng};
use rug::{Integer, Rational};

pub trait Rec<R, Args> {
    type Output;

    fn rec_with_ran(&mut self, reconstructor: R, rng: impl Rng)
        -> Self::Output;

    fn rec(&mut self, reconstructor: R) -> Self::Output {
        self.rec_with_ran(reconstructor, thread_rng())
    }
}

pub trait Eval<T>: TryEval<T> {
    fn eval(&self, pt: &T) -> <Self as TryEval<T>>::Output {
        self.try_eval(pt).unwrap()
    }
}

pub trait TryEval<T> {
    type Output;

    fn try_eval(&self, pt: &T) -> Option<Self::Output>;
}

// custom Zero trait that doesn't require Add
pub trait Zero {
    fn zero() -> Self;
    fn is_zero(&self) -> bool;
}

impl<const P: u64> Zero for Z64<P> {
    fn zero() -> Self {
        <Z64<P> as num_traits::Zero>::zero()
    }

    fn is_zero(&self) -> bool {
        num_traits::Zero::is_zero(self)
    }
}

impl<T: Zero, const N: usize> Zero for [T; N] {
    fn zero() -> Self {
        array_init::array_init(|_| T::zero())
    }

    fn is_zero(&self) -> bool {
        self.iter().all(Zero::is_zero)
    }
}

// custom One trait that doesn't require Mul
pub trait One {
    fn one() -> Self;
    fn is_one(&self) -> bool;
}

impl<const P: u64> One for Z64<P> {
    fn one() -> Self {
        <Z64<P> as num_traits::One>::one()
    }

    fn is_one(&self) -> bool {
        num_traits::One::is_one(self)
    }
}

impl One for NonZeroUsize {
    fn one() -> Self {
        unsafe { Self::new_unchecked(1) }
    }

    fn is_one(&self) -> bool {
        *self == <Self as One>::one()
    }
}

macro_rules! impl_zero_one {
    ( $( $x:ty ),* ) => {
        $(
            impl Zero for $x {
                fn zero() -> Self {
                    num_traits::Zero::zero()
                }

                fn is_zero(&self) -> bool {
                    num_traits::Zero::is_zero(self)
                }
            }

            impl One for $x {
                fn one() -> Self {
                    num_traits::One::one()
                }

                fn is_one(&self) -> bool {
                    num_traits::One::is_one(self)
                }
            }
        )*
    };
}

impl_zero_one!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize, f32, f64, Integer, Rational
);

pub trait Cardinality {
    fn cardinality() -> Option<usize>;
}

macro_rules! impl_cardinality_from_bits {
    ( $( $x:ty ),* ) => {
        $(
            impl Cardinality for $x {
                fn cardinality() -> Option<usize> {
                    Some(1usize << <$x>::BITS)
                }
            }
        )*
    };
}

impl_cardinality_from_bits!(i8, i16, i32, u8, u16, u32);

impl<const P: u64> Cardinality for Z64<P> {
    fn cardinality() -> Option<usize> {
        Some(P as usize)
    }
}

pub trait WithVars<'a, V> {
    type Output;

    fn with_vars(&'a self, vars: V) -> Self::Output;
}

pub trait Shift<S> {
    fn shift(self, shift: S) -> Self;
}
