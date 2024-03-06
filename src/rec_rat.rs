use std::ops::ControlFlow;

use ffnt::Z64;
use log::{debug, trace, warn};
use paste::paste;
use rand::Rng;
use rug::Integer;
use seq_macro::seq;

use crate::{
    algebra::{
        poly::flat::FlatPoly,
        rat::{Rat, NoneError}
    },
    rec::{
        rat::{
            finite::{cuyt_lee::{self, RatRecMod}, linear::RecLinear},
            ffrat::{FFRat, combine_crt_rat}
        },
        primes::LARGE_PRIMES
    },
    traits::{Rec, TryEval, Zero},
};

seq!(N in 0..114 {
    paste! { const [<P N>]: u64 = LARGE_PRIMES[N]; }
});

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct RatRec {
    extra_pts: usize,
}

impl RatRec {
    pub fn new(extra_pts: usize) -> Self {
        Self { extra_pts }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum ReconstructionStatus {
    FirstRat,
    FirstNumPoly,
    FirstDenPoly,
    Rat,
    Done,
}

impl From<cuyt_lee::ReconstructionStatus> for ReconstructionStatus {
    fn from(source: cuyt_lee::ReconstructionStatus) -> Self {
        use ReconstructionStatus::*;
        match source {
            cuyt_lee::ReconstructionStatus::Rat => FirstRat,
            cuyt_lee::ReconstructionStatus::NumPoly => FirstNumPoly,
            cuyt_lee::ReconstructionStatus::DenPoly => FirstDenPoly,
            cuyt_lee::ReconstructionStatus::Done => Done,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Needed<const N: usize> {
    Pt([Z64<P0>; N]),
    Pts(Vec<[Z64<P0>; N]>),
    Any(usize),
}

impl<const N: usize> From<cuyt_lee::Needed<P0, N>> for Needed<N> {
    fn from(source: cuyt_lee::Needed<P0, N>) -> Self {
        use Needed::*;
        match source {
            cuyt_lee::Needed::Pt(pt) => Pt(pt),
            cuyt_lee::Needed::Pts(pts) => Pts(pts),
        }
    }
}

seq! {N in 2..=16 {
    paste! {
        use crate::rec::rat::finite::cuyt_lee::[<RatRecMod N>];

        #[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
        pub struct [<RatRec N>] {
            rec: [<RatRecMod N>]<P0>,
            rat: FFRat<N>,
            res: Result<Rat<FlatPoly<Integer, N>>, NoneError>,
            extra_pts: usize,
            done: bool,
        }

        impl [<RatRec N>] {
            pub fn new(extra_pts: usize) -> Self {
                Self::with_shift(extra_pts, [Z64::zero(); N])
            }

            pub fn with_shift(extra_pts: usize, shift: [Z64<P0>; N]) -> Self {
                Self {
                    rec: [<RatRecMod N>]::with_shift(extra_pts, shift),
                    rat: Default::default(),
                    res: Err(NoneError {  }),
                    extra_pts,
                    done: false,
                }
            }

            pub fn add_pt(&mut self, z: [Z64<P0>; N], q_z: Z64<P0>) -> ControlFlow<(), Needed<N>> {
                use ControlFlow::*;
                use ReconstructionStatus::*;
                match self.status() {
                    FirstRat | FirstNumPoly | FirstDenPoly => match self.rec.add_pt(z, q_z) {
                        Continue(needed) => Continue(needed.into()),
                        Break(()) => {
                            self.finish_first_mod_rec();
                            self.ask_for_new_mod()
                        }
                    },
                    Rat => {
                        warn!("Passed single point when a slice of points was requested");
                        self.ask_for_new_mod()
                    },
                    Done => Break(()),
                }
            }

            pub fn add_pts<'a, const P: u64>(
                &mut self,
                pts: &'a [([Z64<P>; N], Z64<P>)],
            ) -> ControlFlow<(), Needed<N>> {
                use ControlFlow::*;
                use ReconstructionStatus::*;
                match self.status() {
                    FirstRat | FirstNumPoly | FirstDenPoly => {
                        if P != P0 {
                            todo!("request points in correct characteristic")
                        }
                        // cast to the correct type
                        // SAFETY: we just checked that the type is actually &'a [([Z64<P0>; N], Z64<P0>)]
                        //         the lifetime is explicitly set to match the original one
                        let pts: &'a [_] = unsafe {
                            std::slice::from_raw_parts(pts.as_ptr() as _, pts.len())
                        };

                        match self.rec.add_pts(pts) {
                            Continue(needed) => Continue(needed.into()),
                            Break(()) => {
                                self.finish_first_mod_rec();
                                self.ask_for_new_mod()
                            }
                        }
                    },
                    Rat => self.rec_rat_mod_from_pts(pts),
                    Done => Break(()),
                }
            }

            pub fn rec_rat_mod_from_pts<const P: u64>(
                &mut self,
                pts: &[([Z64<P>; N], Z64<P>)],
            ) -> ControlFlow<(), Needed<N>> {
                // if we can already reproduce the points we are done
                if let Ok(res) = self.res.as_ref() {
                    if pts.len() >= self.extra_pts {
                        let sample_same = pts.iter()
                            .all(|(pt, val)| res.try_eval(pt) == Some(*val));

                        if sample_same {
                            self.done = true;
                            return ControlFlow::Break(());
                        }
                    }
                }

                let Some(next_mod_rec) = self.rat.rat.rec_linear(pts.iter().copied()) else {
                    return self.ask_for_new_mod();
                };
                let next_mod_rec = normalise_coeff(next_mod_rec);
                debug!("Finished reconstruction modulo {P}: {next_mod_rec}");
                let mod_rec = std::mem::take(&mut self.rat);
                self.rat = combine_crt_rat(mod_rec, next_mod_rec);
                self.res  = (&self.rat).try_into();
                self.ask_for_new_mod()
            }

            pub fn status(&self) -> ReconstructionStatus {
                use ReconstructionStatus::{Rat, Done};
                if self.rat.modulus.is_zero() {
                    // still reconstructing in first characteristic
                    let status = ReconstructionStatus::from(self.rec.status());
                    assert_ne!(status, Done);
                    status
                } else if self.done {
                    Done
                } else {
                    Rat
                }
            }

            fn finish_first_mod_rec(&mut self) {
                debug_assert_eq!(self.rec.status(), cuyt_lee::ReconstructionStatus::Done);
                let rec = std::mem::replace(&mut self.rec, [<RatRecMod N>]::new(0));
                let rec = normalise_coeff(rec.into_rat());
                self.rat = FFRat::from(rec);
                debug!("Finished reconstruction modulo {P0}: {}", self.rat.rat);
                self.res = (&self.rat).try_into();
                if let Ok(rat) = self.res.as_ref() {
                    if rat.is_zero() {
                        self.done = true;
                    }
                }
            }

            fn ask_for_new_mod(&self) -> ControlFlow<(), Needed<N>> {
                let ncoeff = self.rat.rat.num().len() + self.rat.rat.den().len() - 1;
                if ncoeff == 0 {
                    ControlFlow::Break(())
                } else {
                    ControlFlow::Continue(Needed::Any(ncoeff))
                }
            }

            pub fn request_next_arg(
                &self,
                bad_z: [Z64<P0>; N]
            ) -> [Z64<P0>; N] {
                self.rec.request_next_arg(bad_z)
            }

            pub fn into_rat(self) -> Option<Rat<FlatPoly<Integer, N>>> {
                self.res.ok()
            }
        }
    }
}}

impl<F, const N: usize> Rec<RatRec, [Integer; N]> for F
where
    F: Rec<
        RatRecMod,
        [[Z64<P0>; N]; 1],
        Output = Option<Rat<FlatPoly<Z64<P0>, N>>>,
    >,
    F: TryEval<[Z64<P1>; N], Output = Z64<P1>>,
    F: TryEval<[Z64<P2>; N], Output = Z64<P2>>,
    F: TryEval<[Z64<P3>; N], Output = Z64<P3>>,
    F: TryEval<[Z64<P4>; N], Output = Z64<P4>>,
    F: TryEval<[Z64<P5>; N], Output = Z64<P5>>,
    F: TryEval<[Z64<P6>; N], Output = Z64<P6>>,
    F: TryEval<[Z64<P7>; N], Output = Z64<P7>>,
    F: TryEval<[Z64<P8>; N], Output = Z64<P8>>,
    F: TryEval<[Z64<P9>; N], Output = Z64<P9>>,
    F: TryEval<[Z64<P10>; N], Output = Z64<P10>>,
    F: TryEval<[Z64<P11>; N], Output = Z64<P11>>,
    F: TryEval<[Z64<P12>; N], Output = Z64<P12>>,
    F: TryEval<[Z64<P13>; N], Output = Z64<P13>>,
    F: TryEval<[Z64<P14>; N], Output = Z64<P14>>,
    F: TryEval<[Z64<P15>; N], Output = Z64<P15>>,
    F: TryEval<[Z64<P16>; N], Output = Z64<P16>>,
    F: TryEval<[Z64<P17>; N], Output = Z64<P17>>,
    F: TryEval<[Z64<P18>; N], Output = Z64<P18>>,
    F: TryEval<[Z64<P19>; N], Output = Z64<P19>>,
    F: TryEval<[Z64<P20>; N], Output = Z64<P20>>,
    F: TryEval<[Z64<P21>; N], Output = Z64<P21>>,
    F: TryEval<[Z64<P22>; N], Output = Z64<P22>>,
    F: TryEval<[Z64<P23>; N], Output = Z64<P23>>,
    F: TryEval<[Z64<P24>; N], Output = Z64<P24>>,
    F: TryEval<[Z64<P25>; N], Output = Z64<P25>>,
    F: TryEval<[Z64<P26>; N], Output = Z64<P26>>,
    F: TryEval<[Z64<P27>; N], Output = Z64<P27>>,
    F: TryEval<[Z64<P28>; N], Output = Z64<P28>>,
    F: TryEval<[Z64<P29>; N], Output = Z64<P29>>,
    F: TryEval<[Z64<P30>; N], Output = Z64<P30>>,
    F: TryEval<[Z64<P31>; N], Output = Z64<P31>>,
    F: TryEval<[Z64<P32>; N], Output = Z64<P32>>,
    F: TryEval<[Z64<P33>; N], Output = Z64<P33>>,
    F: TryEval<[Z64<P34>; N], Output = Z64<P34>>,
    F: TryEval<[Z64<P35>; N], Output = Z64<P35>>,
    F: TryEval<[Z64<P36>; N], Output = Z64<P36>>,
    F: TryEval<[Z64<P37>; N], Output = Z64<P37>>,
    F: TryEval<[Z64<P38>; N], Output = Z64<P38>>,
    F: TryEval<[Z64<P39>; N], Output = Z64<P39>>,
    F: TryEval<[Z64<P40>; N], Output = Z64<P40>>,
    F: TryEval<[Z64<P41>; N], Output = Z64<P41>>,
    F: TryEval<[Z64<P42>; N], Output = Z64<P42>>,
    F: TryEval<[Z64<P43>; N], Output = Z64<P43>>,
    F: TryEval<[Z64<P44>; N], Output = Z64<P44>>,
    F: TryEval<[Z64<P45>; N], Output = Z64<P45>>,
    F: TryEval<[Z64<P46>; N], Output = Z64<P46>>,
    F: TryEval<[Z64<P47>; N], Output = Z64<P47>>,
    F: TryEval<[Z64<P48>; N], Output = Z64<P48>>,
    F: TryEval<[Z64<P49>; N], Output = Z64<P49>>,
    F: TryEval<[Z64<P50>; N], Output = Z64<P50>>,
    F: TryEval<[Z64<P51>; N], Output = Z64<P51>>,
    F: TryEval<[Z64<P52>; N], Output = Z64<P52>>,
    F: TryEval<[Z64<P53>; N], Output = Z64<P53>>,
    F: TryEval<[Z64<P54>; N], Output = Z64<P54>>,
    F: TryEval<[Z64<P55>; N], Output = Z64<P55>>,
    F: TryEval<[Z64<P56>; N], Output = Z64<P56>>,
    F: TryEval<[Z64<P57>; N], Output = Z64<P57>>,
    F: TryEval<[Z64<P58>; N], Output = Z64<P58>>,
    F: TryEval<[Z64<P59>; N], Output = Z64<P59>>,
    F: TryEval<[Z64<P60>; N], Output = Z64<P60>>,
    F: TryEval<[Z64<P61>; N], Output = Z64<P61>>,
    F: TryEval<[Z64<P62>; N], Output = Z64<P62>>,
    F: TryEval<[Z64<P63>; N], Output = Z64<P63>>,
    F: TryEval<[Z64<P64>; N], Output = Z64<P64>>,
    F: TryEval<[Z64<P65>; N], Output = Z64<P65>>,
    F: TryEval<[Z64<P66>; N], Output = Z64<P66>>,
    F: TryEval<[Z64<P67>; N], Output = Z64<P67>>,
    F: TryEval<[Z64<P68>; N], Output = Z64<P68>>,
    F: TryEval<[Z64<P69>; N], Output = Z64<P69>>,
    F: TryEval<[Z64<P70>; N], Output = Z64<P70>>,
    F: TryEval<[Z64<P71>; N], Output = Z64<P71>>,
    F: TryEval<[Z64<P72>; N], Output = Z64<P72>>,
    F: TryEval<[Z64<P73>; N], Output = Z64<P73>>,
    F: TryEval<[Z64<P74>; N], Output = Z64<P74>>,
    F: TryEval<[Z64<P75>; N], Output = Z64<P75>>,
    F: TryEval<[Z64<P76>; N], Output = Z64<P76>>,
    F: TryEval<[Z64<P77>; N], Output = Z64<P77>>,
    F: TryEval<[Z64<P78>; N], Output = Z64<P78>>,
    F: TryEval<[Z64<P79>; N], Output = Z64<P79>>,
    F: TryEval<[Z64<P80>; N], Output = Z64<P80>>,
    F: TryEval<[Z64<P81>; N], Output = Z64<P81>>,
    F: TryEval<[Z64<P82>; N], Output = Z64<P82>>,
    F: TryEval<[Z64<P83>; N], Output = Z64<P83>>,
    F: TryEval<[Z64<P84>; N], Output = Z64<P84>>,
    F: TryEval<[Z64<P85>; N], Output = Z64<P85>>,
    F: TryEval<[Z64<P86>; N], Output = Z64<P86>>,
    F: TryEval<[Z64<P87>; N], Output = Z64<P87>>,
    F: TryEval<[Z64<P88>; N], Output = Z64<P88>>,
    F: TryEval<[Z64<P89>; N], Output = Z64<P89>>,
    F: TryEval<[Z64<P90>; N], Output = Z64<P90>>,
    F: TryEval<[Z64<P91>; N], Output = Z64<P91>>,
    F: TryEval<[Z64<P92>; N], Output = Z64<P92>>,
    F: TryEval<[Z64<P93>; N], Output = Z64<P93>>,
    F: TryEval<[Z64<P94>; N], Output = Z64<P94>>,
    F: TryEval<[Z64<P95>; N], Output = Z64<P95>>,
    F: TryEval<[Z64<P96>; N], Output = Z64<P96>>,
    F: TryEval<[Z64<P97>; N], Output = Z64<P97>>,
    F: TryEval<[Z64<P98>; N], Output = Z64<P98>>,
    F: TryEval<[Z64<P99>; N], Output = Z64<P99>>,
    F: TryEval<[Z64<P100>; N], Output = Z64<P100>>,
    F: TryEval<[Z64<P101>; N], Output = Z64<P101>>,
    F: TryEval<[Z64<P102>; N], Output = Z64<P102>>,
    F: TryEval<[Z64<P103>; N], Output = Z64<P103>>,
    F: TryEval<[Z64<P104>; N], Output = Z64<P104>>,
    F: TryEval<[Z64<P105>; N], Output = Z64<P105>>,
    F: TryEval<[Z64<P106>; N], Output = Z64<P106>>,
    F: TryEval<[Z64<P107>; N], Output = Z64<P107>>,
    F: TryEval<[Z64<P108>; N], Output = Z64<P108>>,
    F: TryEval<[Z64<P109>; N], Output = Z64<P109>>,
    F: TryEval<[Z64<P110>; N], Output = Z64<P110>>,
    F: TryEval<[Z64<P111>; N], Output = Z64<P111>>,
    F: TryEval<[Z64<P112>; N], Output = Z64<P112>>,
    F: TryEval<[Z64<P113>; N], Output = Z64<P113>>,
{
    type Output = Option<Rat<FlatPoly<Integer, N>>>;

    fn rec_with_ran(&mut self, rec: RatRec, mut rng: impl Rng) -> Self::Output {
        let rec_mod = RatRecMod::new(rec.extra_pts);

        debug!("Trying rational reconstruction over characteristic {P0}");
        let mod_rec = Rec::<RatRecMod, [[Z64<P0>; N]; 1]>::rec_with_ran(
            self, rec_mod, &mut rng,
        )?;
        let ncoeff = mod_rec.num().len() + mod_rec.den().len() - 1;
        let mod_rec = normalise_coeff(mod_rec);
        debug!("Reconstructed {mod_rec}");
        let mut mod_rec = FFRat::from(mod_rec);
        let mut res: Result<Rat<FlatPoly<Integer, N>>, _> =
            (&mod_rec).try_into();

        seq!( M in 1..114 {{
            const P: u64 = paste!{ [<P M>] };
            debug!("Trying rational reconstruction over characteristic {P}");
            // next batch of sampling points
            let mut pts = Vec::with_capacity(ncoeff);
            while pts.len() < ncoeff {
                let pt: [Z64<P>; N] = [(); N].map(|_| rng.gen());
                if let Some(val) = self.try_eval(&pt) {
                    pts.push((pt, val))
                }
            }
            // if we can already reproduce the points we are done
            if let Ok(res_ref) = res.as_ref() {
                let sample_same = pts.iter()
                    .all(|(pt, val)| res_ref.try_eval(pt) == Some(*val));

                if sample_same {
                    return Some(res.unwrap());
                }
            }

            let next_mod_rec = mod_rec.rat.rec_linear(pts)?;
            let next_mod_rec = normalise_coeff(next_mod_rec);
            debug!("Reconstructed {next_mod_rec}");
            mod_rec = combine_crt_rat(mod_rec, next_mod_rec);
            res  = (&mod_rec).try_into();
        }});

        debug!("Rational reconstruction failed");
        trace!("Final value: {res:?}");
        None
    }
}

fn normalise_coeff<const P: u64, const N: usize>(
    rat: Rat<FlatPoly<Z64<P>, N>>,
) -> Rat<FlatPoly<Z64<P>, N>> {
    let Some(term) = rat.num().terms().last() else {
        return Zero::zero();
    };
    let norm = term.coeff.inv();
    let (mut num, mut den) = rat.into_num_den();
    num *= norm;
    den *= norm;
    Rat::from_num_den_unchecked(num, den)
}
