use std::ops::ControlFlow::{self, *};

use log::{debug, trace};
use thiserror::Error;

use crate::{
    algebra::{poly::flat::{FlatMono, FlatPoly}, rat::Rat}, rec::rat::{
        finite::degree_rec::{self, DegreeRec}, thiele_univar, util,
    }, traits::Zero, Integer, Z64
};

pub type IntRat<const N: usize> = Rat<FlatPoly<Integer, N>>;

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Rec<const P: u64, const N: usize> {
    extra_pts: usize,
    rec: DegreeOrScaledRec<P, N>,
    scalings: Scalings<N>,
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
enum DegreeOrScaledRec<const P: u64, const N: usize> {
    DegreeRec(DegreeRec<P, N>),
    ScaledRec(thiele_univar::Rec),
}

impl<const P: u64, const N: usize> Rec<P, N> {
    pub fn new(extra_pts: usize) -> Self {
        let rec = if N == 1 {
            DegreeOrScaledRec::ScaledRec(thiele_univar::Rec::new(extra_pts))
        } else {
            DegreeOrScaledRec::DegreeRec(DegreeRec::new(extra_pts))
        };
        Self {
            extra_pts,
            rec,
            scalings: Scalings::default(),
        }
    }

    pub fn add_pt<const Q: u64>(
        &mut self,
        z: [Z64<Q>; N],
        q_z: Z64<Q>
    ) -> Result<ControlFlow<IntRat<N>, Status<N>>, Error<P, N>> {
        use DegreeOrScaledRec::*;
        use Status::*;
        trace!("Adding: f({z:?}) = {q_z}");
        match self.rec {
            DegreeRec(ref mut rec) => {
                if P != Q {
                    return Err(util::RecError::Mod{
                        expected: P,
                        found: Q,
                    }.into());
                }
                // SAFETY: always safe thanks to the assert
                let (z, q_z) = unsafe {
                    assert_eq!(P, Q);
                    let z = z.map(|z| Z64::new_unchecked(u64::from(z)));
                    (z, std::mem::transmute(q_z))
                };
                match rec.add_pt(z, q_z)? {
                    Continue(n) => Ok(Continue(Varying(n))),
                    Break(powers) => {
                        self.set_scalings(powers.map(|p| p.into_iter().max().unwrap()));
                        Ok(Continue(Scaling))
                    },
                }
            },
            ScaledRec(ref mut rec) => {
                self.scalings.check_pt(z)?;
                rec.add_pt([z[0]], q_z)?;
                let res = match rec.status() {
                    thiele_univar::Status::NeedNextPt => Continue(Scaling),
                    thiele_univar::Status::NeedNextMod => Continue(NextMod),
                    thiele_univar::Status::Done => {
                        let rat = std::mem::replace(
                            rec,
                            thiele_univar::Rec::new(self.extra_pts)
                        ).into_rat().unwrap();
                        debug!("Finished reconstruction: {rat}");
                        Break(self.scalings.undo(rat))
                    },
                };
                Ok(res)
            }
        }
    }

    pub fn set_scalings(&mut self, scalings: [u32; N]) {
        self.scalings = Scalings(scalings);
        debug!("Set scalings: {:?}", self.scalings);
        self.rec = DegreeOrScaledRec::ScaledRec(
            thiele_univar::Rec::new(self.extra_pts)
        );
    }

    pub fn scale<const Q: u64>(&self, z0: Z64<Q>) -> [Z64<Q>; N] {
        self.scalings.scale(z0)
    }

}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct Scalings<const N: usize>([u32; N]);

impl<const N: usize> Scalings<N> {
    fn scale<const P: u64>(&self, z0: Z64<P>) -> [Z64<P>; N] {
        let mut res = [z0; N];
        for n in 1..N {
            res[n] = res[n - 1].powu(self.0[n - 1] as u64)
        }
        res
    }

    fn check_pt<const P: u64, const Q: u64>(
        &self,
        z: [Z64<P>; N],
    ) -> Result<(), Error<Q, N>> {
        let expected = self.scale(z[0]);
        if z != expected {
            Err(Error::Scaling {
                expected: expected.map(|z| u64::from(z)),
                got: z.map(|z| u64::from(z))
            })
        } else {
            Ok(())
        }
    }

    fn undo(
        &self,
        rat: Rat<FlatPoly<Integer, 1>>
    ) -> Rat<FlatPoly<Integer, N>> {
        if rat.is_zero() {
            return Rat::zero();
        }
        let (num, den) = rat.into_num_den();
        let num = FlatPoly::from_terms(
            num.into_terms()
                .into_iter()
                .map(|t| term_from_scalings(t, &self.0))
                .collect()
        );
        let den = FlatPoly::from_terms(
            den.into_terms()
                .into_iter()
                .map(|t| term_from_scalings(t, &self.0))
                .collect()
        );
        Rat::from_num_den_unchecked(num, den)
    }
}

impl<const N: usize> Default for Scalings<N> {
    fn default() -> Self {
        Self([0; N])
    }
}

fn term_from_scalings<const N: usize>(
    t: FlatMono<Integer, 1>,
    scalings: &[u32; N]
) -> FlatMono<Integer, N> {
    let FlatMono{
        coeff,
        powers,
    } = t;
    let [mut pow] = powers;
    if N == 1 {
        return FlatMono { powers: [pow; N], coeff }
    }
    let mut res_powers = [0; N];
    for (n, s) in scalings.iter().enumerate() {
        res_powers[n] = pow % s;
        pow /= s;
    }
    FlatMono { powers: res_powers, coeff }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Status<const N: usize> {
    NextMod,
    Scaling,
    Varying(usize),
}

#[derive(Debug, Error)]
pub enum Error<const P: u64, const N: usize> {
    #[error("Error reconstructing rational function degrees")]
    DegreeErr(#[from] degree_rec::Error<P, N>),
    #[error("Error reconstructing rational function")]
    RecErr(#[from] util::RecError),
    #[error("Wrong Scaling: expected {expected:?}, got {got:?}")]
    Scaling{
        expected: [u64; N],
        got: [u64; N],
    },
}

#[cfg(test)]
mod tests {
    use super::*;
    use ::rand::{Rng, SeedableRng};

    use crate::rec::primes::LARGE_PRIMES;
    use crate::algebra::poly::flat::FlatMono;
    use crate::traits::{One, TryEval, Zero};
    use paste::paste;
    use rug::integer::Order;
    use seq_macro::seq;

    const NTESTS: usize = 100;
    const EXTRA_SAMPLES: usize = 10;
    const MAX_TERMS: usize = 5;
    const MAX_POW: u32 = 2;
    const MAX_COEFF_U64_POW: usize = 2;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    fn rand_int(mut rng: impl Rng) -> Integer {
        let len = rng.gen_range(1..=MAX_COEFF_U64_POW);
        let parts = Vec::from_iter((0..len).map(|_| rng.gen::<u64>()));
        let int = Integer::from_digits(&parts, Order::Lsf);
        if rng.gen() {
            int
        } else {
            -int
        }
    }

    fn rand_term<const N: usize>(mut rng: impl Rng) -> FlatMono<Integer, N> {
        let pow = [(); N].map(|_| rng.gen_range(0..=MAX_POW));
        FlatMono::new(rand_int(rng), pow)
    }

    fn rand_poly<const N: usize>(mut rng: impl Rng) -> FlatPoly<Integer, N> {
        let nterms = rng.gen_range(0..=MAX_TERMS);
        FlatPoly::from_terms((0..nterms).map(|_| rand_term(&mut rng)).collect())
    }

    fn rand_rat<const N: usize>(
        mut rng: impl Rng,
    ) -> Rat<FlatPoly<Integer, N>> {
        let mut den = FlatPoly::zero();
        // TODO: implement shift so that everything works with non-vanishing starting power
        while den.is_zero() || !den.term(0).powers.is_zero() {
            den = rand_poly(&mut rng);
        }
        let num = rand_poly(rng);
        Rat::from_num_den_unchecked(num, den)
    }

    seq!(NVARS in 1..=3 {
        paste! {
            #[test]
            fn [<rec_rat_explicit_ NVARS>]() {
                use ControlFlow::*;
                use Status::*;
                log_init();

                let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
                for _ in 0..NTESTS {
                    const P: u64 = LARGE_PRIMES[0];

                    let orig = rand_rat::<NVARS>(&mut rng);
                    eprintln!("trying to reconstruct {orig}");
                    let mut rec: Rec<P, NVARS> = Rec::new(1);
                    let (mut z, mut q_z) = std::iter::repeat_with(|| {
                        let z: [Z64<P>; NVARS] = [(); NVARS].map(|_| rng.gen());
                        orig.try_eval(&z).map(|q_z| (z, q_z))
                    })
                        .flatten()
                        .next()
                        .unwrap();
                    let mut status;
                    loop {
                        status = rec.add_pt(z, q_z).unwrap();
                        if let Continue(Varying(n)) = status {
                            z[n] += Z64::one();
                            q_z = loop {
                                if let Some(val) = orig.try_eval(&z) {
                                    break val;
                                }
                                z[n] += Z64::one() ;
                            };
                        } else {
                            break;
                        }
                    }
                    let rec = match status {
                        Break(rec) => rec,
                        Continue(Scaling) => 'rec: {
                            seq!{ N in 0..20 {{
                                const P: u64 = LARGE_PRIMES[N];
                                let mut z0: Z64<P> = rng.gen();

                                loop {
                                    z0 += Z64::one();
                                    let mut z = rec.scale(z0);
                                    let q_z = loop {
                                        if let Some(val) = orig.try_eval(&z) {
                                            break val;
                                        }
                                        z0 += Z64::one();
                                        z = rec.scale(z0);
                                    };
                                    let status = rec.add_pt(z, q_z).unwrap();
                                    match status {
                                        Break(rec) => break 'rec rec,
                                        Continue(NextMod) => break,
                                        Continue(Scaling) => {},
                                        _ => unreachable!("Unexpected reconstruction return value")
                                    }
                                }
                            }}}
                            panic!("Need more than 20 characteristics!");
                        }
                        _ => unreachable!("Unexpected reconstruction return value")
                    };
                    // assert_eq!(rec.status, Status::Done);
                    // let rec = rec.into_rat().unwrap();
                    // eprintln!("reconstructed {rec}");

                    for _ in 0..EXTRA_SAMPLES {
                        let n = [(); NVARS].map(|_| rand_int(&mut rng));
                        debug!("check at x = {n:?}");
                        let orig_val = orig.try_eval(&n);
                        let rec_val = rec.try_eval(&n);
                        debug!("{orig_val:?} == {rec_val:?}");
                        assert!((orig_val == rec_val) || orig_val.is_none())
                    }
                }
            }
        }
    });
}
