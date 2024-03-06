use ffnt::Z64;
use rug::{integer::IntegerExt64, ops::RemRounding, Integer, Rational};

use crate::{arr::Arr, traits::{One, Zero}};

pub(crate) fn merge_crt<const P: u64>(c: &mut Integer, d: Z64<P>, modulus: &Integer) {
    // Terms [a, b] in Bézout's identity a * N + b * M = 1 for coprime N, M
    let ExtendedGCDResult { gcd, bezout } =
        extended_gcd(modulus.clone(), Integer::from(P));
    debug_assert!(gcd.is_one());
    let [a, b] = bezout;
    debug_assert!((Integer::from(&a * modulus) + &b * P).is_one());
    //  x % N = c
    //  x % M = d
    // has the solution x = c - (c - d) * a * N
    let shift = Integer::from(&*c - u64::from(d));
    *c -= shift * a * modulus;
    if c.is_negative() {
        let new_mod = Integer::from(P * modulus);
        *c = std::mem::take(c).rem_euc(new_mod);
    }
    debug_assert!(!c.is_negative());
    debug_assert_eq!(c.mod_u64(P), u64::from(d));
}

pub(crate) fn rat_reconstruct(
    coeff: &Integer,
    modulus: &Integer
) -> Option<Rational> {
    // TODO: code duplication
    let max_bound = Integer::from(modulus / 2);
    // TODO: make configurable
    let max_den = Integer::from(max_bound.root_64_ref(5));
    wang_reconstruct(
        coeff.to_owned(),
        modulus.to_owned(),
        &(max_bound / &max_den),
        &max_den,
    )
}

fn wang_reconstruct(
    value: Integer,
    modulus: Integer,
    max_num: &Integer,
    max_den: &Integer,
) -> Option<Rational> {
    let mut v = Arr([modulus, Integer::zero()]);
    let mut w = Arr([value, Integer::one()]);
    while &w[0] > max_num {
        let q = Integer::from(&v[0] / &w[0]);
        let z = v - w.clone() * &q;
        (v, w) = (w, z);
    }
    if w[1] < 0 {
        w = -w;
    }
    if &w[1] < max_den && w[0].clone().gcd(&w[1]) == 1 {
        let [num, den] = w.0;
        Some(unsafe { Rational::from_canonical(num, den) })
    } else {
        None
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct ExtendedGCDResult {
    gcd: Integer,
    bezout: [Integer; 2],
}

fn extended_gcd(a: Integer, b: Integer) -> ExtendedGCDResult {
    let mut old_r = a;
    let mut r = b;
    let mut old_s = Integer::one();
    let mut s = Integer::zero();
    let mut old_t = Integer::zero();
    let mut t = Integer::one();

    while !r.is_zero() {
        let quotient = Integer::from(&old_r / &r);
        let sub = Integer::from(&quotient * &r);
        (old_r, r) = (r, old_r - sub);
        let sub = Integer::from(&quotient * &s);
        (old_s, s) = (s, old_s - sub);
        let mut sub = quotient;
        sub *= &t;
        (old_t, t) = (t, old_t - sub);
    }
    ExtendedGCDResult {
        gcd: old_r,
        bezout: [old_s, old_t],
    }
}
