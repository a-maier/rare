use std::ops::ControlFlow;

use ffnt::Z64;
use log::debug;
use paste::paste;
use thiserror::Error;

use crate::{
    algebra::{
        poly::flat::FlatPoly,
        rat::{NoneError, Rat},
    },
    rec::{
        primes::LARGE_PRIMES,
        rat::{ffrat::FFRat, finite::cuyt_lee::Needed},
    },
    traits::TryEval,
    Integer,
};

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Status<const P: u64, const N: usize> {
    #[default]
    NextMod,
    Needed(Needed<P, N>),
    Done,
}

impl<const P: u64, const N: usize> From<Needed<P, N>> for Status<P, N> {
    fn from(value: Needed<P, N>) -> Self {
        Status::Needed(value)
    }
}

const P0: u64 = LARGE_PRIMES[0];

macro_rules! impl_rat_rec {
    ( $($n:literal), * ) => {
        $(
            paste! {
                use crate::rec::rat::finite::cuyt_lee::[<RatRecMod $n>];

                #[derive(Debug)]
                /// Rational function reconstruction using the method by Cuyt and Lee
                pub struct [<Rec $n>] {
                    extra_pts: usize,
                    rec: [<RatRecMod $n>]<P0>,
                    rat: FFRat<$n>,
                    res: Result<Rat<FlatPoly<Integer, $n>>, NoneError>,
                    status: Status<P0, $n>,
                    modulus: u64,
                    sample_agree: usize,
                    shift: [i64; $n],
                }

                impl [<Rec $n>] {
                    /// Construct a rational function reconstructor
                    ///
                    /// `extra_pts` is the number of additional redundant points used
                    /// to check the result. `shift` should be a shift of variables
                    /// ensuring that there is no pole when all arguments of the
                    /// rational function are zero.
                    pub fn with_shift(
                        extra_pts: usize,
                        shift: [i64; $n],
                    ) -> Self {
                        Self {
                            extra_pts,
                            rec: [<RatRecMod $n>]::new(0),
                            rat: Default::default(),
                            res: Err(NoneError {}),
                            status: Default::default(),
                            modulus: Default::default(),
                            sample_agree: 0,
                            shift,
                        }
                    }

                    pub fn add_pt<'a, const P: u64>(
                        &'a mut self,
                        z: [Z64<P>; $n],
                        q_z: Z64<P>,
                    ) -> Result<(), RecError> {
                        if self.status == Status::Done || self.sample_done(z, q_z) {
                            return Ok(());
                        }
                        if self.status == Status::NextMod {
                            debug!("New characteristic: {P}");
                            self.modulus = P;
                            let shift = self.shift.map(Z64::<P>::new);
                            let rec = [<RatRecMod $n>]::<P>::with_shift(self.extra_pts(), shift);
                            // Safety: RatRecMod<P> should have the same binary layout for all P
                            self.rec = unsafe {
                                std::mem::transmute(rec)
                            };
                        }
                        if P != self.modulus {
                            return Err(RecError::Mod{expected: self.modulus, found: P});
                        }
                        let rec: &'a mut [<RatRecMod $n>]<P> = unsafe {
                            std::mem::transmute(&mut self.rec)
                        };
                        let res = rec.add_pt(z, q_z);
                        match res {
                            ControlFlow::Break(()) => self.next_mod::<P>(),
                            ControlFlow::Continue(needed) => {
                                let status: Status<P, $n> = needed.into();
                                // Safety:
                                assert_eq!(P, self.modulus);
                                self.status = unsafe {
                                    std::mem::transmute(status)
                                };
                            },
                        }
                        Ok(())
                    }

                    pub fn add_pts<'a, const P: u64>(
                        &mut self,
                        pts: &[([Z64<P>; $n], Z64<P>)],
                    ) -> Result<(), RecError> {
                        if self.status == Status::Done
                            || pts.iter().map(|(z, q_z)| self.sample_done(*z, *q_z)).last().unwrap_or(false) {
                                return Ok(());
                            }
                        debug_assert_ne!(self.status, Status::NextMod);
                        if P != self.modulus {
                            return Err(RecError::Mod{expected: self.modulus, found: P});
                        }
                        let rec: &'a mut [<RatRecMod $n>]<P> = unsafe {
                            std::mem::transmute(&mut self.rec)
                        };
                        let res = rec.add_pts(pts);
                        match res {
                            ControlFlow::Break(()) => self.next_mod::<P>(),
                            ControlFlow::Continue(needed) => {
                                let status: Status<P, $n> = needed.into();
                                // Safety:
                                assert_eq!(P, self.modulus);
                                self.status = unsafe {
                                    std::mem::transmute(status)
                                };
                            },
                        }
                        Ok(())
                    }

                    fn sample_done<const P: u64>(
                        &mut self,
                        z: [Z64<P>; $n],
                        q_z: Z64<P>,
                    ) -> bool {
                        if let Ok(res) = self.res.as_ref() {
                            let our_q_z = res.try_eval(&z);
                            if our_q_z == Some(q_z) {
                                self.sample_agree += 1;
                                if self.sample_agree >= self.extra_pts {
                                    self.modulus = 0;
                                    debug!("done");
                                    self.status = Status::Done;
                                    return true;
                                }
                            } else {
                                self.sample_agree = 0;
                            }
                        }
                        false
                    }

                    /// The current reconstruction status with suggestions for the next step
                    ///
                    /// `P` has to match the current `modulus()`
                    pub fn status<'a, const P: u64>(&'a self) -> Result<&'a Status<P, $n>, RecError> {
                        if !matches!(self.status, Status::Needed(_)) || P == self.modulus {
                            // Safety:
                            // - We just checked that the modulus is either irrelevant or matches
                            // - The lifetime is conserved
                            let status: &'a Status<P, $n> = unsafe{
                                std::mem::transmute(&self.status)
                            };
                            Ok(status)
                        } else {
                            Err(RecError::Mod { expected: self.modulus, found: P })
                        }
                    }

                    /// The number of extra points used to validate the reconstruction
                    pub fn extra_pts(&self) -> usize {
                        self.extra_pts
                    }

                    /// Extract the reconstructed rational function
                    ///
                    /// Returns `None` if `status()` is not `Done`
                    pub fn into_rat(self) -> Option<Rat<FlatPoly<Integer, $n>>> {
                        self.res.ok()
                    }

                    /// The characteristic of the prime field over which we are
                    /// currently reconstructing
                    ///
                    /// Returns zero if the reconstruction has either finished,
                    /// i.e. `status()` is `Done`, or not started yet.
                    pub fn modulus(&self) -> u64 {
                        self.modulus
                    }

                    fn next_mod<const P: u64>(&mut self) {
                        let next_mod_rec = std::mem::replace(
                            &mut self.rec,
                            [<RatRecMod $n>]::new(0)
                        );
                        // SAFETY: the actual modulus is self.modulus
                        assert_eq!(P, self.modulus());
                        let next_mod_rec: [<RatRecMod $n>]<P> = unsafe{
                            std::mem::transmute(next_mod_rec)
                        };
                        let next_mod_rec = next_mod_rec.into_rat();
                        debug!("Finished reconstruction modulo {}: {next_mod_rec}", self.modulus);
                        if self.rat.modulus.is_zero() {
                            self.rat = next_mod_rec.into();
                        } else {
                            self.rat.merge_crt(next_mod_rec);
                            self.res = (&self.rat).try_into();
                        }
                        self.status = Status::NextMod;
                    }
                }
            }
        )*
    };
}

impl_rat_rec! {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}

#[derive(Copy, Clone, Debug, Error, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum RecError {
    #[error("Wrong characteristic: expected {expected}, got {found}")]
    Mod { expected: u64, found: u64 },
}

#[cfg(test)]
mod tests {
    // TODO: code duplication with thiele_linear
    use super::*;
    use ::rand::{Rng, SeedableRng};

    use crate::rec::primes::LARGE_PRIMES;
    use crate::sparse_poly::FlatMono;
    use crate::traits::Zero;
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
        while den.is_zero() {
            den = rand_poly(&mut rng);
        }
        let num = rand_poly(rng);
        Rat::from_num_den_unchecked(num, den)
    }

    seq!(NVARS in 2..=3 {
        paste! {
            #[test]
            fn [<rec_rat_explicit_ NVARS>]() {
                log_init();

                let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
                for _ in 0..NTESTS {
                    let orig = rand_rat::<NVARS>(&mut rng);
                    let shift = [(); NVARS].map(|_| rng.gen());
                    eprintln!("trying to reconstruct {orig}");
                    let mut rec = [<Rec NVARS>]::with_shift(1, shift);
                    'rec: {
                        seq!{ N in 0..20 {{
                            const P: u64 = LARGE_PRIMES[N];
                            debug!("mod {P}");
                            let z: [Z64<P>; NVARS] = [(); NVARS].map(|_| rng.gen());
                            let q_z = orig.try_eval(&z).unwrap();
                            rec.add_pt(z, q_z).unwrap();
                            loop {
                                match rec.status().unwrap() {
                                    Status::Done => break 'rec,
                                    Status::NextMod => break,
                                    Status::Needed(Needed::Pts(pts)) => {
                                        let pts = Vec::from_iter(
                                            pts.into_iter()
                                                .filter_map(
                                                    |z| orig.try_eval(z).map(|q_z| (*z, q_z))
                                                )
                                        );
                                        rec.add_pts(&pts).unwrap();
                                    },
                                    Status::Needed(Needed::Pt(z)) => {
                                        let q_z: Z64<P> = orig.try_eval(z).unwrap();
                                        rec.add_pt(*z, q_z).unwrap();
                                    }
                                }
                            }
                        }}
                        }
                    }
                    assert_eq!(rec.status, Status::Done);
                    let rec = rec.into_rat().unwrap();
                    eprintln!("reconstructed {rec}");

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
