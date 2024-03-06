use std::ops::{ControlFlow, Range};

use ffnt::Z64;
use log::{debug, trace};
use paste::paste;
use rand::Rng;
use rug::Integer;
use seq_macro::seq;
use thiserror::Error;

use crate::{
    algebra::{
        poly::{dense::DensePoly, flat::{FlatPoly, FlatMono}},
        rat::{Rat, NoneError}
    },
    rec::rat::finite::thiele::ThieleRec,
    rec_linear::{RecLinear, UNIT, Unit},
    rec_rat::{combine_crt_rat, FFRat, LARGE_PRIMES},
    traits::{TryEval, Zero, Rec},
};

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Needed<const P: u64, const N: usize> {
    Pt([Z64<P>; N]),
    Any(usize),
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum ReconstructionStatus {
    Degrees,
    Rat,
    Done,
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct RatRecLinear {
    extra_pts: usize,
}

impl RatRecLinear {
    pub fn new(extra_pts: usize) -> Self {
        Self {extra_pts}
    }

    pub fn extra_pts(&self) -> usize {
        self.extra_pts
    }
}

fn nterms_with_max_pows<const N: usize>(
    max_pows: [u32; N]
) -> usize {
    max_pows.iter().map(|&n| n as usize).product()
}

fn gen_poly_with_max_pows<const N: usize>(max_pows: [u32; N]) -> FlatPoly<Unit, N> {
    let num_terms = nterms_with_max_pows(max_pows);
    let mut terms = Vec::with_capacity(num_terms);
    for mut i in 0..(num_terms as u32) {
        let mut pows = [0u32; N];
        for (pow, max) in pows.iter_mut().zip(&max_pows) {
            *pow = i % *max;
            i /= max;
        }
        terms.push(FlatMono::new(UNIT, pows));
    }
    // TODO: would be better to generate them in the
    // correct order instead of sorting
    terms.sort_unstable();
    FlatPoly::from_raw_terms(terms)
}

seq!(N in 0..114 {
    paste! { const [<P N>]: u64 = LARGE_PRIMES[N]; }
});

impl<F, const N: usize> Rec<RatRecLinear, [Integer; N]> for F
where
    F: TryEval<[Z64<P0>; N], Output = Z64<P0>>,
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

    fn rec_with_ran(&mut self, rec: RatRecLinear, mut rng: impl Rng) -> Self::Output {
        let z: [Z64<P0>; N] = [(); N].map(|_| rng.gen());
        let mut num_pows = [0; N];
        let mut den_pows = [0; N];
        for i in 0..N {
            let mut q = |y: Z64<P0>| {
                let mut z = z;
                z[i] += y;
                self.try_eval(&z)
            };
            let res = q.rec_with_ran(
                ThieleRec::new(rec.extra_pts()),
                &mut rng
            )?;
            let res: Rat<DensePoly<_>> = res.into();
            num_pows[i] = res.num().len() as u32;
            den_pows[i] = res.den().len() as u32;
            if num_pows[i] == 0 {
                return Some(Zero::zero())
            }
        }
        debug!("Numerator powers: {num_pows:?}");
        debug!("Denominator powers: {den_pows:?}");
        // TODO: reuse points
        let num = gen_poly_with_max_pows(num_pows);
        let den = gen_poly_with_max_pows(den_pows);
        let pts = Vec::from_iter(
            std::iter::repeat_with(|| {
                let z: [Z64<P0>; N] = [(); N].map(|_| rng.gen());
                self.try_eval(&z).map(|q_z| (z, q_z))
            })
                .flatten()
                .take(num.len() + den.len() + rec.extra_pts() - 1)
        );
        let ansatz = Rat::from_num_den_unchecked(num, den);
        debug!("Trying rational reconstruction over characteristic {P0}");
        let mod_rec = ansatz.rec_linear(pts)?;
        let ncoeff = mod_rec.num().len() + mod_rec.den().len() - 1;

        // TODO: code duplication with `Rec<RatRec, [Integer; N]> for F`
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
            debug!("Reconstructed {next_mod_rec}");
            mod_rec = combine_crt_rat(mod_rec, next_mod_rec);
            res  = (&mod_rec).try_into();
        }});

        debug!("Rational reconstruction failed");
        trace!("Final value: {res:?}");
        None
    }
}

#[derive(Error, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FailedRec<const N: usize> {
    #[error("Need points for reconstruction")]
    Empty,
    #[error("Need more points for arguments {:?} in characteristic {modulus}", format_args(*args, *ncoord))]
    MoreAt{
        modulus: u64,
        args: [u64; N],
        ncoord: usize,
    },
    #[error("Need more points in characteristic {modulus}. Estimated total number: {}", nexpected_to_string(*.nexpected))]
    MorePts{
        modulus: u64,
        nexpected: Option<usize>,
    },
    #[error("Need approximately {0} points in new characteristic")]
    MoreMods(usize),
    #[error("Unknown modulus {0}: is not in `LARGE_PRIMES`")]
    UnknownMod(u64)
}

#[derive(Debug)]
struct UnknownMod(u64);

impl<const N: usize> From<UnknownMod> for FailedRec<N> {
    fn from(source: UnknownMod) -> Self {
        Self::UnknownMod(source.0)
    }
}

fn format_args<const N: usize>(
    args: [u64; N],
    ncoord: usize
) -> [String; N] {
    let mut res = args.map(|a| a.to_string());
    res[ncoord] = "x".to_string();
    res
}

fn nexpected_to_string(
    n: Option<usize>
) -> String {
    if let Some(n) = n {
        n.to_string()
    } else {
        "unknown".to_string()
    }
}

pub fn rec_linear_from_pts<const N: usize>(
    pts: &mut [(u64, Vec<([u64; N], u64)>)],
    extra_pts: usize,
) -> Result<Rat<FlatPoly<Integer, N>>, FailedRec<N>>
{
    use FailedRec::*;
    if pts.is_empty() {
        return Err(Empty)
    }

    let mut num_pows = [0; N];
    let mut den_pows = [0; N];
    for ncoord in 0..N {
        let (modulus, indices) = get_longest_streak(pts, ncoord);
        let streak = &pts.iter().find(|(m, _)| *m == modulus).unwrap().1[indices];
        let pts = streak.iter().map(|(z, q_z)| (z[ncoord], *q_z));
        let rec_pow = try_rec_max_num_den_pow(modulus, pts, extra_pts);
        let Some([max_num_pow, max_den_pow]) = rec_pow? else {
            let args = streak[0].0;
            return Err(MoreAt{modulus, args, ncoord})
        };
        num_pows[ncoord] = max_num_pow;
        den_pows[ncoord] = max_den_pow;
    }
    debug!("Numerator powers: {num_pows:?}");
    debug!("Denominator powers: {den_pows:?}");
    let nexpected = nterms_with_max_pows(num_pows) + nterms_with_max_pows(den_pows) + extra_pts - 1;
    let (modulus, most_pts) = pts.iter().max_by_key(|(_, pts)| pts.len()).unwrap();
    let modulus = *modulus;
    let Some(mut mod_rec) = rec_with_pows(most_pts, modulus, num_pows, den_pows)? else {
        debug!("Reconstruction over mod {modulus} failed");
        return Err(MorePts{modulus, nexpected: Some(nexpected)})
    };
    if max_pows_reached(&mod_rec.rat) != [num_pows, den_pows] {
        debug!("Reconstruction over mod {modulus} did not reach powers {num_pows:?}, {den_pows:?}");
        return Err(MorePts{modulus, nexpected: None})
    }
    let mut res: Result<Rat<FlatPoly<Integer, N>>, _> = (&mod_rec).try_into();
    let first_modulus = modulus;

    // TODO: might be better to iterate over LARGE_PRIMES
    for (m, pts) in pts {
        use ControlFlow::*;
        if *m == first_modulus {
            continue;
        }
        if let Break(res) = add_rec(&mut res, &mut mod_rec, pts, *m, extra_pts)? {
            return res;
        }
    }
    let nexpected = mod_rec.rat.num().len() + mod_rec.rat.den().len() - 1;
    Err(MoreMods(nexpected))
}

fn add_rec<'a, const N: usize>(
    res: &mut Result<Rat<FlatPoly<Integer, N>>, NoneError>,
    rat: &mut FFRat<N>,
    pts: &'a [([u64; N], u64)],
    modulus: u64,
    extra_pts: usize
) -> Result<ControlFlow<Result<Rat<FlatPoly<Integer, N>>, FailedRec<N>>>, UnknownMod> {
    seq!{ Q in 0..114 {
          if modulus == LARGE_PRIMES[Q] {
              const P: u64 = LARGE_PRIMES[Q];
              let pts: &'a [([Z64<P>; N], Z64<P>)] =  unsafe {
                  std::slice::from_raw_parts(pts.as_ptr() as _, pts.len())
              };
              return Ok(add_rec_mod(res, rat, pts, extra_pts));
          }
    } }
    Err(UnknownMod(modulus))
}

fn add_rec_mod<const P: u64, const N: usize>(
    res: &mut Result<Rat<FlatPoly<Integer, N>>, NoneError>,
    rat: &mut FFRat<N>,
    pts: &[([Z64<P>; N], Z64<P>)],
    extra_pts: usize
) -> ControlFlow<Result<Rat<FlatPoly<Integer, N>>, FailedRec<N>>> {
    use ControlFlow::*;
    // if we can already reproduce the points we are done
    if let Ok(res_ref) = res.as_ref() {
        if pts.len() > extra_pts {
            let sample_same = pts.iter()
                .all(|(pt, val)| res_ref.try_eval(pt) == Some(*val));

            if sample_same {
                let res = std::mem::take(res.as_mut().unwrap());
                return Break(Ok(res));
            }
        }
    }

    debug!("Trying rational reconstruction over characteristic {P}");
    let Some(next_mod_rec) = rat.rat.rec_linear(pts.iter().copied()) else {
        let nexpected = rat.rat.num().len() + rat.rat.den().len() - 1;
        return Break(Err(FailedRec::MorePts { modulus: P, nexpected: Some(nexpected) }))
    };
    debug!("Reconstructed {next_mod_rec}");
    let [num_pows, den_pows] = max_pows_reached(&rat.rat);
    if max_pows_reached(&next_mod_rec) != [num_pows, den_pows] {
        debug!("Reconstruction over mod {P} did not reach powers {num_pows:?}, {den_pows:?}");
        return Break(Err(FailedRec::MorePts{modulus: P, nexpected: None}))
    }
    *rat = combine_crt_rat(std::mem::take(rat), next_mod_rec);
    *res  = (&*rat).try_into();
    Continue(())
}

fn max_pows_reached<T, const N: usize>(
    rat: &Rat<FlatPoly<T, N>>,
) -> [[u32; N]; 2] {
    use std::cmp::max;
    let mut num_pows = [0; N];
    for term in rat.num().terms().iter() {
        for i in 0..N {
            num_pows[i] = max(term.powers[i] + 1, num_pows[i]);
        }
    }
    let mut den_pows = [0; N];
    for term in rat.den().terms().iter() {
        for i in 0..N {
            den_pows[i] = max(term.powers[i] + 1, den_pows[i]);
        }
    }
    [num_pows, den_pows]
}

fn rec_with_pows<'a, const N: usize>(
    pts: &'a [([u64; N], u64)],
    modulus: u64,
    num_pows: [u32; N],
    den_pows: [u32; N],
) -> Result<Option<FFRat<N>>, UnknownMod> {
    seq!{ Q in 0..114 {
          if modulus == LARGE_PRIMES[Q] {
              const P: u64 = LARGE_PRIMES[Q];
              let pts: &'a [([Z64<P>; N], Z64<P>)] =  unsafe {
                  std::slice::from_raw_parts(pts.as_ptr() as _, pts.len())
              };
              return Ok(rec_with_pows_mod(pts, num_pows, den_pows));
          }
    } }
    Err(UnknownMod(modulus))
}

fn rec_with_pows_mod<const P: u64, const N: usize>(
    pts: &[([Z64<P>; N], Z64<P>)],
    num_pows: [u32; N],
    den_pows: [u32; N],
) -> Option<FFRat<N>> {
    let num = gen_poly_with_max_pows(num_pows);
    let den = gen_poly_with_max_pows(den_pows);
    let ansatz = Rat::from_num_den_unchecked(num, den);
    debug!("Trying rational reconstruction over characteristic {P}");
    let mod_rec = ansatz.rec_linear(pts.iter().copied())?;
    debug!("Reconstructed {mod_rec}");
    Some(FFRat::from(mod_rec))
}

fn get_longest_streak<const N: usize>(
    pts: &mut [(u64, Vec<([u64; N], u64)>)],
    ncoord: usize
) -> (u64, Range<usize>) {
    let mut longest_streak = (0, 0..0);
    for (modulus, pts) in pts {
        pts.sort_unstable_by_key(|(coord, _)| {
            let mut coord = *coord;
            coord[ncoord] = 0;
            coord
        });
        let streak = longest_streak_pos_by(
            &pts,
            |a, b| {
                let mut a_coord = a.0;
                a_coord[ncoord] = 0;
                let mut b_coord = b.0;
                b_coord[ncoord] = 0;
                a_coord == b_coord
            }
        );
        if streak.len() > longest_streak.1.len() {
            longest_streak = (*modulus, streak);
        }
    }
    longest_streak
}

fn try_rec_max_num_den_pow(
    modulus: u64,
    pts: impl IntoIterator<Item = (u64, u64)>,
    extra_pts: usize
) -> Result<Option<[u32; 2]>, UnknownMod> {
    seq!{ Q in 0..114 {
          if modulus == LARGE_PRIMES[Q] {
              const P: u64 = LARGE_PRIMES[Q];
              let pts = pts.into_iter().map(|(z, q_z)| unsafe{
                  (Z64::<P>::new_unchecked(z), Z64::<P>::new_unchecked(q_z))
              });
              return Ok(try_rec_max_num_den_pow_mod(pts, extra_pts));
          }
    } }
    Err(UnknownMod(modulus))
}

fn try_rec_max_num_den_pow_mod<const P: u64>(
    pts: impl IntoIterator<Item = (Z64<P>, Z64<P>)>,
    extra_pts: usize
) -> Option<[u32; 2]> {
    let rec = ThieleRec::new(extra_pts).rec_from_seq(pts)?;
    let rec: Rat<DensePoly<_>> = rec.into();
    Some([rec.num().len() as u32, rec.den().len() as u32])
}

fn longest_streak_pos_by<T, F>(
    slice: &[T],
    mut is_eq: F
) -> Range<usize>
where
    F: FnMut(&T, &T) -> bool
{
    if slice.is_empty() {
        return 0..0;
    }
    let mut longest_streak = 0..0;
    let mut streak_start = 0;
    for (n, e) in slice.iter().enumerate() {
        if !is_eq(e, &slice[streak_start]) {
            if n - streak_start > longest_streak.len() {
                longest_streak = streak_start..n;
            }
            streak_start = n;
        }
    }
    if slice.len() - streak_start > longest_streak.len() {
        longest_streak = streak_start..slice.len();
    }
    longest_streak
}

#[cfg(test)]
mod tests {
    use super::*;
    use ::rand::{Rng, SeedableRng};

    use crate::rec_rat::LARGE_PRIMES;
    use crate::sparse_poly::FlatMono;
    use crate::traits::Zero;
    use rug::integer::Order;

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

    seq!(NVARS in 2..=3 {
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
                        RatRecLinear::new(1),
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

}
