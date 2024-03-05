use std::ops::ControlFlow;

use ffnt::Z64;
use log::{debug, trace, warn};
use paste::paste;
use rand::Rng;
use rug::{integer::IntegerExt64, ops::RemRounding, Integer, Rational};
use seq_macro::seq;

use crate::{
    algebra::{
        poly::flat::{FlatMono, FlatPoly},
        rat::{Rat, NoneError}
    },
    arr::Arr,
    rec_linear::RecLinear,
    rec_rat_mod::{self, RatRecMod},
    traits::{One, Rec, TryEval, Zero},
};

pub const LARGE_PRIMES: [u64; 114] = [
    1152921504606846883,
    1152921504606846869,
    1152921504606846803,
    1152921504606846797,
    1152921504606846719,
    1152921504606846697,
    1152921504606846607,
    1152921504606846581,
    1152921504606846577,
    1152921504606846523,
    1152921504606846419,
    1152921504606846397,
    1152921504606846347,
    1152921504606846307,
    1152921504606846281,
    1152921504606846269,
    1152921504606846259,
    1152921504606846251,
    1152921504606846223,
    1152921504606846199,
    1152921504606846179,
    1152921504606846097,
    1152921504606846043,
    1152921504606845993,
    1152921504606845977,
    1152921504606845849,
    1152921504606845839,
    1152921504606845789,
    1152921504606845777,
    1152921504606845683,
    1152921504606845657,
    1152921504606845647,
    1152921504606845527,
    1152921504606845503,
    1152921504606845473,
    1152921504606845471,
    1152921504606845399,
    1152921504606845387,
    1152921504606845321,
    1152921504606845317,
    1152921504606845269,
    1152921504606845213,
    1152921504606845167,
    1152921504606845161,
    1152921504606845147,
    1152921504606845101,
    1152921504606845027,
    1152921504606844957,
    1152921504606844933,
    1152921504606844913,
    1152921504606844849,
    1152921504606844829,
    1152921504606844811,
    1152921504606844787,
    1152921504606844741,
    1152921504606844717,
    1152921504606844691,
    1152921504606844591,
    1152921504606844519,
    1152921504606844513,
    1152921504606844447,
    1152921504606844417,
    1152921504606844411,
    1152921504606844289,
    1152921504606844279,
    1152921504606844267,
    1152921504606844243,
    1152921504606844177,
    1152921504606844127,
    1152921504606843871,
    1152921504606843863,
    1152921504606843779,
    1152921504606843733,
    1152921504606843731,
    1152921504606843707,
    1152921504606843703,
    1152921504606843607,
    1152921504606843547,
    1152921504606843527,
    1152921504606843493,
    1152921504606843479,
    1152921504606843431,
    1152921504606843421,
    1152921504606843371,
    1152921504606843353,
    1152921504606843347,
    1152921504606843341,
    1152921504606843313,
    1152921504606843299,
    1152921504606843259,
    1152921504606843233,
    1152921504606843227,
    1152921504606843221,
    1152921504606843199,
    1152921504606843173,
    1152921504606843163,
    1152921504606843107,
    1152921504606843073,
    1152921504606843031,
    1152921504606842989,
    1152921504606842951,
    1152921504606842911,
    1152921504606842909,
    1152921504606842833,
    1152921504606842809,
    1152921504606842791,
    1152921504606842753,
    1152921504606842731,
    1152921504606842669,
    1152921504606842651,
    1152921504606842569,
    1152921504606842527,
    1152921504606842513,
    1152921504606842491,
];

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

impl From<rec_rat_mod::ReconstructionStatus> for ReconstructionStatus {
    fn from(source: rec_rat_mod::ReconstructionStatus) -> Self {
        use ReconstructionStatus::*;
        match source {
            rec_rat_mod::ReconstructionStatus::Rat => FirstRat,
            rec_rat_mod::ReconstructionStatus::NumPoly => FirstNumPoly,
            rec_rat_mod::ReconstructionStatus::DenPoly => FirstDenPoly,
            rec_rat_mod::ReconstructionStatus::Done => Done,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Needed<const N: usize> {
    Pt([Z64<P0>; N]),
    Pts(Vec<[Z64<P0>; N]>),
    Any(usize),
}

impl<const N: usize> From<rec_rat_mod::Needed<P0, N>> for Needed<N> {
    fn from(source: rec_rat_mod::Needed<P0, N>) -> Self {
        use Needed::*;
        match source {
            rec_rat_mod::Needed::Pt(pt) => Pt(pt),
            rec_rat_mod::Needed::Pts(pts) => Pts(pts),
        }
    }
}

seq! {N in 2..=16 {
    paste! {
        use crate::rec_rat_mod::[<RatRecMod N>];

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
                debug_assert_eq!(self.rec.status(), rec_rat_mod::ReconstructionStatus::Done);
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

pub(crate) fn combine_crt_rat<const P: u64, const N: usize>(
    rat: FFRat<N>,
    new_rat: Rat<FlatPoly<Z64<P>, N>>,
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
    let new_terms = new_num.into_iter().chain(new_den);
    for (term, new_term) in terms.zip(new_terms) {
        debug_assert_eq!(term.powers, new_term.powers);
        merge_crt(&mut term.coeff, new_term.coeff, &modulus);
    }
    let num = FlatPoly::from_raw_terms(num);
    let den = FlatPoly::from_raw_terms(den);
    let rat = Rat::from_num_den_unchecked(num, den);
    FFRat {
        rat,
        modulus: modulus * P,
    }
}

fn merge_crt<const P: u64>(c: &mut Integer, d: Z64<P>, modulus: &Integer) {
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

// rational function over finite characteristic that does not necessarily fit in a `Z64`
#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub(crate) struct FFRat<const N: usize> {
    pub(crate) rat: Rat<FlatPoly<Integer, N>>,
    pub(crate) modulus: Integer,
}

impl<const P: u64, const N: usize> From<Rat<FlatPoly<Z64<P>, N>>>
    for FFRat<N>
{
    fn from(source: Rat<FlatPoly<Z64<P>, N>>) -> Self {
        let rat = source.into();
        Self {
            rat,
            modulus: P.into(),
        }
    }
}

impl<'a, const N: usize> TryFrom<&'a FFRat<N>>
    for Rat<FlatPoly<Rational, N>>
{
    type Error = NoneError;

    fn try_from(source: &'a FFRat<N>) -> Result<Self, Self::Error> {
        let mut num = Vec::with_capacity(source.rat.num().len());
        for term in source.rat.num().terms() {
            let rat_coeff = rat_reconstruct(&term.coeff, &source.modulus)
                .ok_or(NoneError {})?;
            num.push(FlatMono::new(rat_coeff, term.powers));
        }
        let num = FlatPoly::from_raw_terms(num);

        let mut den = Vec::with_capacity(source.rat.den().len());
        for term in source.rat.den().terms().iter() {
            let rat_coeff = rat_reconstruct(&term.coeff, &source.modulus)
                .ok_or(NoneError {})?;
            den.push(FlatMono::new(rat_coeff, term.powers));
        }
        let den = FlatPoly::from_raw_terms(den);

        Ok(Rat::from_num_den_unchecked(num, den))
    }
}

impl<'a, const N: usize> TryFrom<&'a FFRat<N>> for Rat<FlatPoly<Integer, N>> {
    type Error = NoneError;

    fn try_from(source: &'a FFRat<N>) -> Result<Self, Self::Error> {
        Rat::<FlatPoly<Rational, N>>::try_from(source).map(|r| r.into())
    }
}

fn rat_reconstruct(coeff: &Integer, modulus: &Integer) -> Option<Rational> {
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

#[cfg(test)]
mod tests {
    use ::rand::SeedableRng;
    use rug::integer::Order;

    use crate::rec_rat_mod::find_shift;

    use super::*;

    const NTESTS: usize = 100;
    const EXTRA_SAMPLES: usize = 10;
    const MAX_TERMS: usize = 5;
    const MAX_POW: u32 = 5;
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
        FlatPoly::from_terms(
            (0..nterms).map(|_| rand_term(&mut rng)).collect(),
        )
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

    seq!(NVARS in 1..=3 {
        paste! {
            #[test]
            fn [<rec_rat_ NVARS>]() {
                log_init();

                let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
                let mut rat_rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(128420185);
                for _ in 0..NTESTS {

                    let orig = rand_rat::<NVARS>(&mut rng);
                    eprintln!("trying to reconstruct {orig}");
                    let rec = orig.clone().rec_with_ran(
                        RatRec::new(1),
                        &mut rat_rng
                    ).unwrap();
                    eprintln!("reconstructed {rec}");

                    for _ in 0..EXTRA_SAMPLES {
                        let n = [(); NVARS].map(|_| rand_int(&mut rng));
                        let orig_val = orig.try_eval(&n);
                        let rec_val = rec.try_eval(&n);
                        assert!((orig_val == rec_val) || orig_val.is_none())
                    }
                }
            }
        }
    });

    seq!(NVARS in 2..=3 {
        paste! {
            #[test]
            fn [<rec_rat_explicit_ NVARS>]() {
                log_init();

                let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
                for _ in 0..NTESTS {

                    let orig = rand_rat::<NVARS>(&mut rng);
                    eprintln!("trying to reconstruct {orig}");
                    let shift = find_shift(|z| orig.try_eval(&z), &mut rng);
                    let mut rec = [<RatRec NVARS>]::with_shift(1, shift);
                    let (z, q_z) = loop {
                        let z: [Z64<P0>; NVARS] = [(); NVARS].map(|_| rng.gen());
                        if let Some(q_z) = orig.try_eval(&z) {
                            break (z, q_z);
                        }
                    };
                    let mut next = rec.add_pt(z, q_z);
                    loop {
                        use ControlFlow::*;
                        use Needed::*;
                        match next {
                            Continue(Pt(mut z)) => {
                                let q_z = loop {
                                    if let Some(q_z) = orig.try_eval(&z) {
                                        break q_z;
                                    }
                                    z = rec.request_next_arg(z);
                                };
                                next = rec.add_pt(z, q_z);
                            },
                            Continue(Pts(zs)) => {
                                let mut pts = Vec::from_iter(
                                    zs.iter().copied()
                                        .filter_map(|z| orig.try_eval(&z).map(|q_z| (z, q_z)))
                                );
                                let mut z = *zs.last().unwrap();
                                while pts.len() < zs.len() {
                                    z = rec.request_next_arg(z);
                                    if let Some(q_z) = orig.try_eval(&z) {
                                        pts.push((z, q_z))
                                    }
                                }
                                next = rec.add_pts(&pts);
                            },
                            Continue(Any(n)) => {
                                seq!{ N in 1..20 {{
                                    const P: u64 = LARGE_PRIMES[N];
                                    let pts = Vec::from_iter(
                                        std::iter::repeat_with(|| {
                                            let z: [Z64<P>; NVARS] = [(); NVARS].map(|_| rng.gen());
                                            orig.try_eval(&z).map(|q_z| (z, q_z))
                                        })
                                            .flatten()
                                            .take(n)
                                    );
                                    let next = rec.add_pts(&pts);
                                    if next == Break(()) {
                                        break;
                                    }
                                    assert_eq!(next, Continue(Any(n)))
                                }}}
                                panic!("Need more than 20 characteristics!");
                            }
                            Break(()) => break,
                        }
                    }
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
