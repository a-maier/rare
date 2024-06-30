use std::{any::Any, fmt::{self, Display}, ops::ControlFlow};

use ffnt::Z64;
use log::{debug, trace};
use paste::paste;
use seq_macro::seq;
use thiserror::Error;

use crate::{
    algebra::{
        poly::flat::FlatPoly,
        rat::{NoneError, Rat},
    },
    rec::{
        primes::LARGE_PRIMES, rat::{ffrat::FFRat, finite::cuyt_lee::{Needed as ModNeeded, z_to_x}}
    },
    traits::TryEval,
    Integer,
};

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Status<const P: u64, const N: usize> {
    #[default]
    NextMod,
    Needed(ModNeeded<P, N>),
    Done,
}

impl<const P: u64, const N: usize> From<ModNeeded<P, N>> for Status<P, N> {
    fn from(value: ModNeeded<P, N>) -> Self {
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
                    sampler: Sampler,
                    shift: [u64; $n],
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
                        shift: [u64; $n],
                    ) -> Self {
                        Self {
                            extra_pts,
                            rec: [<RatRecMod $n>]::new(0),
                            rat: Default::default(),
                            res: Err(NoneError {}),
                            status: Default::default(),
                            modulus: Default::default(),
                            sampler: Sampler::new(extra_pts),
                            shift,
                        }
                    }

                    pub fn add_pt<'a, const P: u64>(
                        &'a mut self,
                        z: [Z64<P>; $n],
                        q_z: Z64<P>,
                    ) -> Result<(), RecError> {
                        if self.status == Status::Done || self.sample(z, q_z) == SampleStatus::Complete {
                            return Ok(());
                        }
                        if self.status == Status::NextMod {
                            self.update_mod::<P>();
                        }
                        if P != self.modulus {
                            return Err(RecError::Mod{expected: self.modulus, found: P});
                        }
                        let rec: &'a mut [<RatRecMod $n>]<P> = unsafe {
                            std::mem::transmute(&mut self.rec)
                        };
                        let res = rec.add_pt(z, q_z);
                        self.update_status(res);
                        Ok(())
                    }

                    pub fn add_pts<'a, const P: u64>(
                        &mut self,
                        mut pts: &[([Z64<P>; $n], Z64<P>)],
                    ) -> Result<(), RecError> {
                        loop {
                            match &self.status {
                                Status::Done => return Ok(()),
                                Status::NextMod => {
                                    if P == self.modulus {
                                        return Ok(());
                                    }
                                    let Some(((z, q_z), rest)) = pts.split_first() else {
                                        return Ok(())
                                    };
                                    pts = rest;
                                    self.add_pt(*z, *q_z)?;
                                },
                                Status::Needed(needed) => {
                                    if P != self.modulus {
                                        return Err(RecError::Mod{expected: self.modulus, found: P});
                                    }
                                    // Safety: we just checked that the actual modulus is P
                                    let needed: &'a ModNeeded<P, $n> = unsafe {
                                        std::mem::transmute(needed)
                                    };
                                    match needed {
                                        ModNeeded::Pt(next_wanted) => {
                                            trace!("next wanted point: {next_wanted:?}");
                                            trace!(
                                                "next wanted point (in x/t): {:?}",
                                                z_to_x(*next_wanted, &self.shift.map(Z64::from))
                                            );
                                            let next_pos = pts.iter()
                                                .position(|pt| &pt.0 == next_wanted);
                                            let Some(next_pos) = next_pos else {
                                                return Ok(())
                                            };
                                            let (z, q_z) = pts[next_pos];
                                            pts = &pts[(next_pos + 1)..];
                                            trace!("adding point at {next_pos}");
                                            self.add_pt(z, q_z)?;
                                        }
                                        ModNeeded::Pts(next_wanted) => {
                                            let sample = pts.iter()
                                                .map(|(q, q_z)| self.sample(*q, *q_z))
                                                .skip_while(|&res| matches!(res, SampleStatus::Success(_)))
                                                .next();
                                            if sample == Some(SampleStatus::Complete) {
                                                return Ok(());
                                            }
                                            trace!("next wanted points: {next_wanted:#?}");
                                            trace!("next wanted points (in x/t):");
                                            for z in next_wanted {
                                                trace!("{:?}", z_to_x(*z, &self.shift.map(Z64::from)));
                                            }
                                            let next_pos = pts.iter().position(|pt| pt.0 == next_wanted[0]);
                                            let Some(next_pos) = next_pos else {
                                                return Ok(())
                                            };
                                            pts = &pts[next_pos..];
                                            let has_wanted_pts = pts.iter()
                                                .take(next_wanted.len())
                                                .map(|pt| pt.0)
                                                .eq(next_wanted.iter().copied());
                                            if !has_wanted_pts {
                                                return Ok(())
                                            }
                                            trace!("adding points from {next_pos} on");
                                            let (wanted, rest) = pts.split_at(next_wanted.len());
                                            pts = rest;
                                            // Safety:
                                            assert_eq!(self.modulus, P);
                                            let rec: &'a mut [<RatRecMod $n>]<P> = unsafe {
                                                std::mem::transmute(&mut self.rec)
                                            };
                                            let res = rec.add_pts(wanted);
                                            self.update_status(res);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    pub fn add_pts_unsorted<const P: u64>(
                        &mut self,
                        pts: &mut [([Z64<P>; $n], Z64<P>)],
                    ) -> Result<(), RecError> {
                        let shift = self.shift.map(Z64::<P>::from);
                        pts.sort_by_cached_key(|pt| z_to_x(pt.0, &shift));
                        self.add_pts(pts)
                    }

                    fn update_status<const P: u64>(
                        &mut self,
                        s: ControlFlow<(), ModNeeded<P, $n>>
                    ) {
                        match s {
                            ControlFlow::Break(()) => self.finish_mod::<P>(),
                            ControlFlow::Continue(needed) => {
                                let status: Status<P, $n> = needed.into();
                                // Safety:
                                assert_eq!(P, self.modulus);
                                self.status = unsafe {
                                    std::mem::transmute(status)
                                };
                            },
                        }
                    }

                    fn sample<const P: u64>(
                        &mut self,
                        z: [Z64<P>; $n],
                        q_z: Z64<P>,
                    ) -> SampleStatus {
                        let Ok(res) = self.res.as_ref() else {
                            return SampleStatus::Failed;
                        };
                        let status = self.sampler.add_pt(&z, q_z, res);
                        if status == SampleStatus::Complete {
                            self.modulus = 0;
                            self.status = Status::Done;
                            debug!("Reconstructed function: {res}");
                        }
                        status
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

                    fn finish_mod<const P: u64>(&mut self) {
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

                    fn update_mod<const P: u64>(&mut self) {
                        debug!("New characteristic: {P}");
                        self.modulus = P;
                        let shift = self.shift.map(Z64::<P>::from);
                        let rec = [<RatRecMod $n>]::<P>::with_shift(self.extra_pts(), shift);
                        // Safety: RatRecMod<P> should have the same binary layout for all P
                        self.rec = unsafe {
                            std::mem::transmute(rec)
                        };
                        if self.sampler.status() == SampleStatus::Failed {
                            self.sampler.reset()
                        }
                    }
                }

                fn [<rec_from_pts $n>](
                    pts: &mut [ModPts<$n>],
                    shift: [u64; $n],
                    extra_pts: usize,
                ) -> Result<Rat<FlatPoly<Integer, $n>>, FailedRec<$n>> {
                    let mut rec = [<Rec $n>]::with_shift(extra_pts, shift);
                    pts.sort_by(
                        |pt1, pt2| (pt2.pts.len(), pt2.modulus)
                            .cmp(&(pt1.pts.len(), pt1.modulus))
                    );
                    for ModPts{modulus, pts} in pts.iter_mut() {
                        seq!{ PP in 0..30 {{
                            if *modulus == LARGE_PRIMES[PP] {
                                const P: u64 = LARGE_PRIMES[PP];
                                let pts = unsafe{
                                    std::mem::transmute(pts.as_mut_slice())
                                };
                                rec.add_pts_unsorted::<P>(pts).unwrap();
                                match rec.status {
                                    Status::NextMod => {},
                                    Status::Needed(needed) => {
                                        let needed = match needed {
                                            ModNeeded::Pt(z) => Needed::Pt(z.map(u64::from)),
                                            ModNeeded::Pts(z) => Needed::Pts(
                                                z.into_iter().map(|z| z.map(u64::from)).collect()
                                            )
                                        };
                                        return Err(FailedRec::MorePts{
                                            modulus: *modulus,
                                            needed
                                        });
                                    },
                                    Status::Done => return Ok(rec.into_rat().unwrap()),
                                }
                                continue;
                            }
                        }}}
                        return Err(FailedRec::UnknownMod(*modulus));
                    }
                    let suggested_next_mod = find_largest_missing_mod(
                        pts.iter().map(|pt| pt.modulus)
                    );
                    if let Some(next_mod) = suggested_next_mod {
                        Err(FailedRec::MoreMods(next_mod))
                    } else {
                        Err(FailedRec::NoModsLeft)
                    }
                }

            }
        )*
    };
}

pub(crate) fn find_largest_missing_mod(
    mods: impl Iterator<Item = u64>
) -> Option<u64> {
    use std::collections::BTreeSet;
    let mut all_mods = BTreeSet::from_iter(LARGE_PRIMES);
    for mmod in mods {
        all_mods.remove(&mmod);
    }
    all_mods.pop_last()
}

// impl_rat_rec! {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
impl_rat_rec! {2, 3, 4, 5, 6, 7, 8}

#[derive(Copy, Clone, Debug, Error, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum RecError {
    #[error("Wrong characteristic: expected {expected}, got {found}")]
    Mod { expected: u64, found: u64 },
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ModPts<const N: usize> {
    pub modulus: u64,
    pub pts: Vec<([u64; N], u64)>,
}

pub fn rec_from_pts<const N: usize>(
    pts: &mut [ModPts<N>],
    shift: [u64; N],
    extra_pts: usize,
) -> Result<Rat<FlatPoly<Integer, N>>, FailedRec<N>> {
    fn cast_res<const M: usize, const N: usize>(
        res: Result<Rat<FlatPoly<Integer, M>>, FailedRec<M>>
    ) -> Result<Rat<FlatPoly<Integer, N>>, FailedRec<N>> {
        let mut res = Some(res);
        let res: &mut dyn Any = &mut res;
        res.downcast_mut()
            .and_then(Option::take)
            .unwrap()
    }

    use std::mem::{transmute, transmute_copy};
    unsafe {
        // Safety: we only ever transmute a type to itself
        match N {
            0 => todo!(),
            1 => {
                use crate::rec::rat::thiele::FailedRec as RecErr;
                match crate::rec::rat::thiele::rec_from_pts(transmute(pts), extra_pts) {
                    Ok(res) => Ok(transmute(res)),
                    Err(RecErr::Empty) => Err(FailedRec::<N>::Empty),
                    Err(RecErr::NoModsLeft) => Err(FailedRec::<N>::NoModsLeft),
                    Err(RecErr::UnknownMod(m)) => Err(FailedRec::<N>::UnknownMod(m)),
                    Err(RecErr::MoreMods(m)) => Err(FailedRec::<N>::MoreMods(m)),
                    Err(RecErr::MorePts(modulus)) => Err(FailedRec::<N>::MorePts {
                        modulus, needed: Needed::Pts(Vec::new())
                    }),
                }
            }
            2 => cast_res(rec_from_pts2(transmute(pts), transmute_copy(&shift), extra_pts)),
            3 => cast_res(rec_from_pts3(transmute(pts), transmute_copy(&shift), extra_pts)),
            4 => cast_res(rec_from_pts4(transmute(pts), transmute_copy(&shift), extra_pts)),
            5 => cast_res(rec_from_pts5(transmute(pts), transmute_copy(&shift), extra_pts)),
            6 => cast_res(rec_from_pts6(transmute(pts), transmute_copy(&shift), extra_pts)),
            7 => cast_res(rec_from_pts7(transmute(pts), transmute_copy(&shift), extra_pts)),
            8 => cast_res(rec_from_pts8(transmute(pts), transmute_copy(&shift), extra_pts)),
            _ => unimplemented!("Multivariate reconstruction with more than 8 variables")
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Needed<const N: usize> {
    Pt([u64; N]),
    Pts(Vec<[u64; N]>),
}

impl<const N: usize> Display for Needed<N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Needed::Pt(z) => write!(f, "{z:?}"),
            Needed::Pts(z) => write!(f, "{z:?}"),
        }
    }
}

#[derive(Error, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FailedRec<const N: usize> {
    #[error("Need points for reconstruction")]
    Empty,
    #[error("Need more points in characteristic {modulus}. Next points: {needed}")]
    MorePts {
        modulus: u64,
        needed: Needed<N>,
    },
    #[error("Need points in new characteristic. Suggested next characteristic: {0}")]
    MoreMods(u64),
    #[error("Need points in new characteristic, but no supported characteristics are left.")]
    NoModsLeft,
    #[error("Unknown modulus {0}: is not in `LARGE_PRIMES`")]
    UnknownMod(u64),
}

#[derive(Debug)]
struct UnknownMod(u64);

impl<const N: usize> From<UnknownMod> for FailedRec<N> {
    fn from(source: UnknownMod) -> Self {
        Self::UnknownMod(source.0)
    }
}

#[cfg(test)]
mod tests {
    // TODO: code duplication with thiele_linear
    use super::*;
    use ::rand::{Rng, SeedableRng};

    use crate::algebra::poly::flat::FlatMono;
    use crate::rec::primes::LARGE_PRIMES;
    use crate::traits::{TryEval, Zero};
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
                                    Status::Needed(ModNeeded::Pts(pts)) => {
                                        let pts = Vec::from_iter(
                                            pts.into_iter()
                                                .filter_map(
                                                    |z| orig.try_eval(z).map(|q_z| (*z, q_z))
                                                )
                                        );
                                        rec.add_pts(&pts).unwrap();
                                    },
                                    Status::Needed(ModNeeded::Pt(z)) => {
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
