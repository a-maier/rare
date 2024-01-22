use std::iter::repeat_with;

use galois_fields::Z64;
use rand::Rng;
use paste::paste;
use seq_macro::seq;

use crate::{dense_poly::{DensePoly, DensePoly1}, rat::Rat, traits::{One, Zero}, sparse_poly::{SparsePoly, SparseMono}};

// generate a random dense polynomial
// the maximum power is 2^m - 1, where m is sampled uniformly between 0 and `n[0]`
pub fn gen_dense_poly1<const P: u64>(
    n: &[u32; 1],
    mut rng: impl Rng
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
                use crate::dense_poly::[<DensePoly $x>];

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

impl_gen_recursive!(1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10);

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
                let mut den = [<gen_dense_poly N>](n, &mut rng).into_coeff();
                den[0] = One::one();
                DensePoly::from_coeff(den)
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
    mut rng: impl Rng
) -> SparsePoly<Z64<P>, N> {
    let nterms = 2usize.pow(rng.gen_range(0..=n));
    let terms = repeat_with(|| gen_sparse_mono(max_pow, &mut rng)).take(nterms).collect();
    SparsePoly::from_terms(terms)
}

pub fn gen_sparse_mono<const P: u64, const N: usize>(
    max_pow: u32,
    mut rng: impl Rng
) -> SparseMono<Z64<P>, N> {
    let coeff = rng.gen_range(1..P);
    let coeff = unsafe{ Z64::new_unchecked(coeff) };
    let powers = [(); N].map(|_| rng.gen_range(0..=max_pow));
    SparseMono{ powers, coeff }
}

// generate a random sparse rational function
// the coefficient of the first term in the denominator is normalised to 1
pub fn gen_sparse_rat<const P: u64, const N: usize>(
    n: u32,
    max_pow: u32,
    mut rng: impl Rng
) -> Rat<SparsePoly<Z64<P>, N>> {
    let num = gen_sparse_poly(n, max_pow, &mut rng);
    let den = if num.is_zero() {
        One::one()
    } else {
        let mut den: SparsePoly<_, N> = Zero::zero();
        while den.is_zero() {
            den = gen_sparse_poly(n, max_pow, &mut rng);
        }
        // normalise
        let mut terms = den.into_terms();
        terms[0].coeff = One::one();
        SparsePoly::from_raw_terms(terms)
    };
    Rat::from_num_den_unchecked(num, den)
}
