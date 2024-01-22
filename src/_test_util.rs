use std::iter::repeat_with;

use galois_fields::Z64;
use rand::Rng;
use paste::paste;
use seq_macro::seq;

use crate::{dense_poly::{DensePoly, DensePoly1}, rat::Rat, traits::{One, Zero}};

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
