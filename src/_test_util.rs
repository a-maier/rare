use std::iter::repeat_with;

use ffnt::Z64;
use paste::paste;
use rand::{distributions::Standard, prelude::Distribution, Rng};
use seq_macro::seq;

use crate::{
    algebra::{
        poly::{dense::{DensePoly, DensePoly1}, flat::{FlatPoly, FlatMono}},
        rat::Rat
    },
    traits::{One, TryEval, Zero},
};

// generate a random dense polynomial
// the maximum power is 2^m - 1, where m is sampled uniformly between 0 and `n[0]`
pub fn gen_dense_poly1<const P: u64>(
    n: &[u32; 1],
    mut rng: impl Rng,
) -> DensePoly<Z64<P>> {
    let max_pow = rng.gen_range(0..=n[0]);
    let nterms = 2usize.pow(max_pow);
    let coeff = repeat_with(|| rng.gen::<Z64<P>>()).take(nterms).collect();
    DensePoly::from_coeff(coeff)
}

macro_rules! impl_gen_recursive {
    ( $($x_minus_one:literal, $x:literal), * ) => {
        $(
            paste! {
                use crate::algebra::poly::dense::[<DensePoly $x>];

                pub fn [<gen_dense_poly $x>]<const P: u64>(
                    n: &[u32; $x],
                    mut rng: impl Rng
                ) -> [<DensePoly $x>]<Z64<P>> {
                    let (n, rest) = n.split_first().unwrap();
                    let max_pow = rng.gen_range(0..=*n);
                    let nterms = 2usize.pow(max_pow);
                    let rest = rest.try_into().unwrap();
                    let coeff = repeat_with(
                        || [<gen_dense_poly $x_minus_one>](&rest, &mut rng)
                    ).take(nterms).collect();
                    DensePoly::from_coeff(coeff)
                }
            }
        )*
    }
}

impl_gen_recursive!(1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10);

seq!( N in 1..=10 {
    paste! {
        pub fn [<gen_dense_rat N>]<const P: u64>(
            n: &[u32; N],
            mut rng: impl Rng
        ) -> Rat<[<DensePoly N>]<Z64<P>>> {
            let num = [<gen_dense_poly N>](n, &mut rng);
            let den = if num.is_zero() {
                One::one()
            } else {
                loop {
                    let den = [<gen_dense_poly N>](n, &mut rng);
                    if !den.is_zero() {
                        break den;
                    }
                }
            };

            Rat::from_num_den_unchecked(num, den)
        }
    }
});

// generate a random sparse polynomial
// the number of terms is 2^m, where m is sampled uniformly between 0 and `n`
pub fn gen_sparse_poly<const P: u64, const N: usize>(
    n: u32,
    max_pow: u32,
    mut rng: impl Rng,
) -> FlatPoly<Z64<P>, N> {
    let nterms = 2usize.pow(rng.gen_range(0..=n));
    let terms = repeat_with(|| gen_sparse_mono(max_pow, &mut rng))
        .take(nterms)
        .collect();
    FlatPoly::from_terms(terms)
}

pub fn gen_sparse_mono<const P: u64, const N: usize>(
    max_pow: u32,
    mut rng: impl Rng,
) -> FlatMono<Z64<P>, N> {
    let coeff = rng.gen_range(1..P);
    let coeff = unsafe { Z64::new_unchecked(coeff) };
    let powers = [(); N].map(|_| rng.gen_range(0..=max_pow));
    FlatMono { powers, coeff }
}

// generate a random sparse rational function
// the coefficient of the first term in the denominator is normalised to 1
pub fn gen_sparse_rat<const P: u64, const N: usize>(
    n: u32,
    max_pow: u32,
    mut rng: impl Rng,
) -> Rat<FlatPoly<Z64<P>, N>> {
    let num = gen_sparse_poly(n, max_pow, &mut rng);
    let den = if num.is_zero() {
        One::one()
    } else {
        loop {
            let den = gen_sparse_poly(n, max_pow, &mut rng);
            if !den.is_zero() {
                break den;
            }
        }
    };
    Rat::from_num_den_unchecked(num, den)
}

pub fn sample_eq<const N: usize, const P: u64>(
    orig: &impl TryEval<[Z64<P>; N], Output = Z64<P>>,
    rec: &impl TryEval<[Z64<P>; N], Output = Z64<P>>,
    mut rng: impl Rng,
) -> bool
where
    Standard: Distribution<[Z64<P>; N]>,
{
    const TESTS: usize = 10;
    for _ in 0..TESTS {
        let pt: [Z64<P>; N] = rng.gen();
        if orig.try_eval(&pt) != rec.try_eval(&pt) {
            return false;
        }
    }
    true
}
