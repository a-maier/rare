use galois_fields::Z64;
use rand::{Rng, thread_rng};

pub trait Rec<R, Args> {
    type Output;

    fn rec_with_ran(
        &mut self,
        reconstructor: R,
        rng: impl Rng
    ) -> Self::Output;

    fn rec(&mut self, reconstructor: R) -> Self::Output {
        self.rec_with_ran(reconstructor, thread_rng())
    }
}

pub trait Eval<T> {
    type Output;

    fn eval(&self, pt: &T) -> Self::Output;
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

pub trait Cardinality {
    fn cardinality()-> Option<usize>;
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
    fn cardinality()-> Option<usize> {
        Some(P as usize)
    }
}

pub trait WithVars<'a, V> {
    type Output;

    fn with_vars(&'a self, vars: V) -> Self::Output;
}
