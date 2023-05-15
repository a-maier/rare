use std::iter::repeat_with;

use galois_fields::Z64;
use lazy_static::lazy_static;
use log::{debug, trace};
use num_traits::Signed;
use rand::Rng;
use rug::{Integer, Rational, integer::IntegerExt64, ops::RemRounding};
use seq_macro::seq;
use paste::paste;

use crate::{traits::{Rec, TryEval, Zero, One}, rat::{Rat, NoneError}, sparse_poly::{SparsePoly, SparseMono}, rec_rat_mod::RatRecMod, arr::Arr};

const LARGE_PRIMES: [u64; 12] = [
    1152921504606846883, 1152921504606846869, 1152921504606846803,
    1152921504606846797, 1152921504606846719, 1152921504606846697,
    1152921504606846607, 1152921504606846581, 1152921504606846577,
    1152921504606846523, 1152921504606846419, 1152921504606846397
];

seq!(N in 0..12 {
    paste! { const [<P N>]: u64 = LARGE_PRIMES[N]; }
});

const fn large_prime_idx(p: u64) -> usize {
    let mut i = 0;
    while LARGE_PRIMES[i] != p {
        i += 1;
        if i >= LARGE_PRIMES.len() {
            panic!("Argument is not a large prime");
        }
    }
    i
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct RatRec {
    extra_pts: usize
}

impl RatRec {
    pub fn new(extra_pts: usize) -> Self {
        Self { extra_pts }
    }
}

impl<F, const N: usize> Rec<RatRec, [Integer; N]> for F
where
    F: Rec<RatRecMod, [[Z64<P0>; N]; 1], Output = Option<Rat<SparsePoly<Z64<P0>, N>>>>,
    F: Rec<RatRecMod, [[Z64<P1>; N]; 1], Output = Option<Rat<SparsePoly<Z64<P1>, N>>>>,
    F: Rec<RatRecMod, [[Z64<P2>; N]; 1], Output = Option<Rat<SparsePoly<Z64<P2>, N>>>>,
    F: Rec<RatRecMod, [[Z64<P3>; N]; 1], Output = Option<Rat<SparsePoly<Z64<P3>, N>>>>,
    F: Rec<RatRecMod, [[Z64<P4>; N]; 1], Output = Option<Rat<SparsePoly<Z64<P4>, N>>>>,
    F: Rec<RatRecMod, [[Z64<P5>; N]; 1], Output = Option<Rat<SparsePoly<Z64<P5>, N>>>>,
    F: Rec<RatRecMod, [[Z64<P6>; N]; 1], Output = Option<Rat<SparsePoly<Z64<P6>, N>>>>,
    F: Rec<RatRecMod, [[Z64<P7>; N]; 1], Output = Option<Rat<SparsePoly<Z64<P7>, N>>>>,
    F: Rec<RatRecMod, [[Z64<P8>; N]; 1], Output = Option<Rat<SparsePoly<Z64<P8>, N>>>>,
    F: Rec<RatRecMod, [[Z64<P9>; N]; 1], Output = Option<Rat<SparsePoly<Z64<P9>, N>>>>,
    F: Rec<RatRecMod, [[Z64<P10>; N]; 1], Output = Option<Rat<SparsePoly<Z64<P10>, N>>>>,
    F: Rec<RatRecMod, [[Z64<P11>; N]; 1], Output = Option<Rat<SparsePoly<Z64<P11>, N>>>>,
{
    type Output = Option<Rat<SparsePoly<Integer, N>>>;

    fn rec_with_ran(
        &mut self,
        rec: RatRec,
        mut rng: impl Rng
    ) -> Self::Output {
        let rec_mod = RatRecMod::new(rec.extra_pts);

        debug!("Trying rational reconstruction over characteristic {P0}");
        let mod_rec = Rec::<RatRecMod, [[Z64<P0>; N]; 1]>::rec_with_ran(
            self, rec_mod, &mut rng
        )?;
        let mut mod_rec = FFRat::from(mod_rec);
        let mut res: Result<Rat<SparsePoly<Integer, N>>, _> = (&mod_rec).try_into();

        seq!( M in 1..12 {{
            const P: u64 = paste!{ [<P M>] };
            debug!("Trying rational reconstruction over characteristic {P}");
            let next_mod_rec = Rec::<RatRecMod, [[Z64<P>; N]; 1]>::rec_with_ran(
                self, rec_mod, &mut rng
            )?;
            // TODO: compare already during reconstruction of `next_mod_rec`
            if let Ok(res_ref) = res.as_ref() {
                let sample_same = repeat_with(|| [(); N].map(|_| rng.gen()))
                    .take(rec.extra_pts)
                    .all(|pt| next_mod_rec.try_eval(&pt) == res_ref.try_eval(&pt));
                if sample_same {
                    return Some(res.unwrap());
                }
            }
            mod_rec = combine_crt_rat(mod_rec, next_mod_rec);
            res  = (&mod_rec).try_into();
        }});

        debug!("Rational reconstruction failed");
        trace!("Final value: {res:?}");
        None
    }
}

fn combine_crt_rat<const P: u64, const N: usize>(
    rat: FFRat<N>,
    new_rat: Rat<SparsePoly<Z64<P>, N>>
) -> FFRat<N> {
    let (num, den) = rat.rat.into_num_den();
    let modulus = rat.modulus;
    let (new_num, new_den) = new_rat.into_num_den();
    debug_assert_eq!(num.len(), new_num.len());
    debug_assert_eq!(den.len(), new_den.len());
    let mut num = num.into_terms();
    let mut den = den.into_terms();
    let new_num = new_num.into_terms();
    let new_den = new_den.into_terms();

    let terms = num.iter_mut().chain(den.iter_mut());
    let new_terms = new_num.into_iter().chain(new_den.into_iter());
    for (term, new_term) in terms.zip(new_terms) {
        debug_assert_eq!(term.powers, new_term.powers);
        merge_crt(&mut term.coeff, new_term.coeff, &modulus);
    }
    let num = SparsePoly::from_raw_terms(num);
    let den = SparsePoly::from_raw_terms(den);
    let rat = Rat::from_num_den_unchecked(num, den);
    FFRat {
        rat,
        modulus: modulus * P,
    }
}

fn merge_crt<const P: u64>(
    c: &mut Integer,
    d: Z64<P>,
    modulus: &Integer
) {
    // TODO: check that `p_idx` gets optimised out
    let p_idx = large_prime_idx(P);
    debug_assert!(p_idx >= 1);
    let shift = Integer::from(&*c - u64::from(d));
    *c -= shift * &BEZOUT[p_idx - 1];
    if c.is_negative() {
        let new_mod = Integer::from(P * modulus);
        *c = std::mem::take(c).rem_euc(new_mod);
    }
    debug_assert!(!c.is_negative());
    debug_assert_eq!(c.mod_u64(P), u64::from(d));
}

// Terms a * N in Bézout's identity a * N + b * M = 1 for coprime N, M
// in BEZOUT[i] the variables are
// M = LARGE_PRIMES[i + 1]
// N = LARGE_PRIMES[0] * ... * LARGE_PRIMES[i]
// These are exactly the terms required for the Chinese Remainder Theorem
// namely
//  x % N = c
//  x % M = d
// has the solution x = c - (c - d) * a * N
lazy_static! {
    static ref BEZOUT: [Integer; 11] = {
        let mut res: [Integer; 11] = Default::default();
        let mut n = Integer::one();
        for i in 0..res.len() {
            n *= LARGE_PRIMES[i];
            let m = LARGE_PRIMES[i + 1];
            let ExtendedGCDResult { gcd, bezout } = extended_gcd(
                n.clone(),
                Integer::from(m)
            );
            debug_assert!(gcd.is_one());
            let [a, b] = bezout;
            debug_assert!((Integer::from(&a * &n) + &b * m).is_one());
            res[i] = a * &n;
        }
        res
    };
}

// rational function over finite characteristic that does not necessarily fit in a `Z64`
#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct FFRat<const N: usize> {
    rat: Rat<SparsePoly<Integer, N>>,
    modulus: Integer
}

impl<const P: u64, const N: usize> From<Rat<SparsePoly<Z64<P>, N>>> for FFRat<N> {
    fn from(source: Rat<SparsePoly<Z64<P>, N>>) -> Self {
        let rat = source.into();
        Self { rat, modulus: P.into()}
    }
}

impl<'a, const N: usize> TryFrom<&'a FFRat<N>> for Rat<SparsePoly<Rational, N>> {
    type Error = NoneError;

    fn try_from(source: &'a FFRat<N>) -> Result<Self, Self::Error> {
        let mut num = Vec::with_capacity(source.rat.num().len());
        for term in source.rat.num().terms() {
            let rat_coeff = rat_reconstruct(&term.coeff, &source.modulus).ok_or(NoneError{})?;
            num.push(SparseMono::new(rat_coeff, term.powers));
        }
        let num = SparsePoly::from_raw_terms(num);

        let mut den = Vec::with_capacity(source.rat.den().len());
        for term in source.rat.den().terms() {
            let rat_coeff = rat_reconstruct(&term.coeff, &source.modulus).ok_or(NoneError{})?;
            den.push(SparseMono::new(rat_coeff, term.powers));
        }
        let den = SparsePoly::from_raw_terms(den);

        Ok(Rat::from_num_den_unchecked(num, den))
    }
}

impl<'a, const N: usize> TryFrom<&'a FFRat<N>> for Rat<SparsePoly<Integer, N>> {
    type Error = NoneError;

    fn try_from(source: &'a FFRat<N>) -> Result<Self, Self::Error> {
        Rat::<SparsePoly<Rational, N>>::try_from(source).map(
            |r| r.into()
        )
    }
}

fn rat_reconstruct(coeff: &Integer, modulus: &Integer) -> Option<Rational> {
    // TODO: code duplication
    let max_bound = Integer::from(modulus / 2);
    let max_den = Integer::from_f64(max_bound.to_f64().powf(1. / 4.)).unwrap();
    wang_reconstruct(coeff.to_owned(), modulus.to_owned(), &(max_bound / &max_den), &max_den)
}

fn wang_reconstruct(
    value: Integer,
    modulus: Integer,
    max_num: &Integer,
    max_den: &Integer
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
        Some(unsafe{ Rational::from_canonical(num, den) })
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
