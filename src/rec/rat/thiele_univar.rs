use std::ops::ControlFlow;

use ffnt::Z64;
use log::debug;
use rug::Integer;
use seq_macro::seq;
use thiserror::Error;

use crate::{
    algebra::{
        poly::flat::FlatPoly,
        rat::{NoneError, Rat},
    },
    rec::{
        primes::LARGE_PRIMES,
        rat::{
            ffrat::FFRat,
            finite::thiele::ThieleRec,
            util::{find_largest_missing_mod, ModPts, RecError},
        },
    },
    traits::TryEval,
};

const P0: u64 = LARGE_PRIMES[0];

/// Reconstruction status
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Status {
    /// At least one more point over the present characteristic is needed
    NeedNextPt,
    /// Points over a new characteristic are needed
    NeedNextMod,
    /// Reconstruction has finished .
    ///
    /// The result can now be extracted with `into_rat()`
    Done,
}

/// Univariate rational reconstruction using Thiele reconstruction
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Rec {
    extra_pts: usize,
    rec: ThieleRec<P0>,
    rat: FFRat<1>,
    res: Result<Rat<FlatPoly<Integer, 1>>, NoneError>,
    modulus: u64,
    sample_agree: usize,
    status: Status,
}

impl Rec {
    pub fn new(extra_pts: usize) -> Self {
        Self {
            extra_pts,
            rec: ThieleRec::new(extra_pts),
            rat: Default::default(),
            res: Err(NoneError {}),
            modulus: Default::default(),
            sample_agree: 0,
            status: Status::NeedNextMod,
        }
    }

    pub fn add_pt<'a, const P: u64>(
        &'a mut self,
        z: [Z64<P>; 1],
        q_z: Z64<P>,
    ) -> Result<(), RecError> {
        let z = z[0];
        self.sample(z, q_z);
        if self.status() == Status::Done {
            return Ok(());
        }
        if !self.rec.rec_started() {
            debug!("New characteristic: {P}");
            self.modulus = P;
        } else if P != self.modulus {
            return Err(RecError::Mod {
                expected: self.modulus,
                found: P,
            });
        }
        // SAFETY:
        // `ThieleRec` should have the same layout independently of `P`
        let rec: &'a mut ThieleRec<P> =
            unsafe { std::mem::transmute(&mut self.rec) };
        if rec.add_pt(z, q_z) == ControlFlow::Break(()) {
            self.status = Status::NeedNextMod;
            let next_mod_rec =
                std::mem::replace(rec, ThieleRec::new(self.extra_pts))
                    .into_rat();
            let next_mod_rec: Rat<FlatPoly<_, 1>> = next_mod_rec.into();
            debug!(
                "Finished reconstruction modulo {}: {next_mod_rec}",
                self.modulus
            );
            if self.rat.modulus.is_zero() {
                self.rat = next_mod_rec.into();
            } else {
                self.rat.merge_crt(next_mod_rec);
            }
            self.res = (&self.rat).try_into();
        } else {
            self.status = Status::NeedNextPt;
        }
        Ok(())
    }

    fn sample<const P: u64>(&mut self, z: Z64<P>, q_z: Z64<P>) {
        if let Ok(res) = self.res.as_ref() {
            let our_q_z = res.try_eval(&[z]);
            if our_q_z == Some(q_z) {
                self.sample_agree += 1;
                if self.sample_agree >= self.extra_pts {
                    self.status = Status::Done;
                }
            } else {
                self.sample_agree = 0;
            }
        }
    }

    /// The current reconstruction status with suggestions for the next step
    pub fn status(&self) -> Status {
        self.status
    }

    /// The number of extra points used to validate the reconstruction
    pub fn extra_pts(&self) -> usize {
        self.extra_pts
    }

    /// Extract the reconstructed rational function
    ///
    /// Returns `None` if `status()` is not `Done`
    pub fn into_rat(self) -> Option<Rat<FlatPoly<Integer, 1>>> {
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
}

pub fn rec_from_pts(
    pts: &mut [ModPts<1>],
    extra_pts: usize,
) -> Result<Rat<FlatPoly<Integer, 1>>, FailedRec> {
    if pts.is_empty() {
        return Err(FailedRec::Empty);
    }
    let mut rec = Rec::new(extra_pts);
    pts.sort_by(|pt1, pt2| {
        (pt2.pts.len(), pt2.modulus).cmp(&(pt1.pts.len(), pt1.modulus))
    });
    for ModPts { modulus, pts } in pts.iter_mut() {
        seq! { PP in 0..114 {{
            if *modulus == LARGE_PRIMES[PP] {
                const P: u64 = LARGE_PRIMES[PP];
                for pt in pts {
                    let z = pt.0.map(|z| z.into());
                    let q_z = Z64::<P>::from(pt.1);
                    rec.add_pt(z, q_z).unwrap();
                    match rec.status() {
                        Status::NeedNextPt => {},
                        Status::NeedNextMod => continue,
                        Status::Done => return Ok(rec.into_rat().unwrap()),
                    }
                }
                if rec.status() == Status::NeedNextPt {
                    return Err(FailedRec::MorePts(P));
                }
                continue;
            }
        }}}
        return Err(FailedRec::UnknownMod(*modulus));
    }
    let suggested_next_mod =
        find_largest_missing_mod(pts.iter().map(|pt| pt.modulus));
    if let Some(next_mod) = suggested_next_mod {
        Err(FailedRec::MoreMods(next_mod))
    } else {
        Err(FailedRec::NoModsLeft)
    }
}

#[derive(Error, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FailedRec {
    #[error("Need points for reconstruction")]
    Empty,
    #[error("Need more points in characteristic {0}")]
    MorePts(u64),
    #[error(
        "Need points in new characteristic. Suggested next characteristic: {0}"
    )]
    MoreMods(u64),
    #[error("Need points in new characteristic, but no supported characteristics are left.")]
    NoModsLeft,
    #[error("Unknown modulus {0}: is not in `LARGE_PRIMES`")]
    UnknownMod(u64),
}

#[cfg(test)]
mod tests {
    use super::*;
    use ::rand::{Rng, SeedableRng};

    use crate::algebra::poly::flat::FlatMono;
    use crate::rec::primes::LARGE_PRIMES;
    use crate::traits::Zero;
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

    fn rand_term(mut rng: impl Rng) -> FlatMono<Integer, 1> {
        let pow = [(); 1].map(|_| rng.gen_range(0..=MAX_POW));
        FlatMono::new(rand_int(rng), pow)
    }

    fn rand_poly(mut rng: impl Rng) -> FlatPoly<Integer, 1> {
        let nterms = rng.gen_range(0..=MAX_TERMS);
        FlatPoly::from_terms((0..nterms).map(|_| rand_term(&mut rng)).collect())
    }

    fn rand_rat(mut rng: impl Rng) -> Rat<FlatPoly<Integer, 1>> {
        let mut den = FlatPoly::zero();
        while den.is_zero() {
            den = rand_poly(&mut rng);
        }
        let num = rand_poly(rng);
        Rat::from_num_den_unchecked(num, den)
    }

    #[test]
    fn rec_rat_explicit() {
        log_init();

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        for _ in 0..NTESTS {
            let orig = rand_rat(&mut rng);
            eprintln!("trying to reconstruct {orig}");
            let mut rec = Rec::new(1);
            'rec: {
                seq! { N in 0..20 {{
                    const P: u64 = LARGE_PRIMES[N];
                    loop {
                        let (z, q_z) = std::iter::repeat_with(|| {
                            let z: [Z64<P>; 1] = [(); 1].map(|_| rng.gen());
                            orig.try_eval(&z).map(|q_z| (z, q_z))
                        })
                            .flatten()
                            .next()
                            .unwrap();
                        rec.add_pt(z, q_z).unwrap();
                        match rec.status {
                            Status::NeedNextPt => {},
                            Status::NeedNextMod => break,
                            Status::Done => break 'rec,
                        }
                    }
                }}}
            }
            assert_eq!(rec.status, Status::Done);
            let rec = rec.into_rat().unwrap();
            eprintln!("reconstructed {rec}");

            for _ in 0..EXTRA_SAMPLES {
                let n = [(); 1].map(|_| rand_int(&mut rng));
                debug!("check at x = {n:?}");
                let orig_val = orig.try_eval(&n);
                let rec_val = rec.try_eval(&n);
                debug!("{orig_val:?} == {rec_val:?}");
                assert!((orig_val == rec_val) || orig_val.is_none())
            }
        }
    }
}
