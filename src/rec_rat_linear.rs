use std::ops::{ControlFlow, Range};

use ffnt::Z64;
use log::debug;
use rug::Integer;
use seq_macro::seq;
use thiserror::Error;

use crate::{
    algebra::{
        poly::{
            dense::DensePoly,
            flat::{FlatMono, FlatPoly},
        },
        rat::{NoneError, Rat},
    },
    rec::{
        primes::LARGE_PRIMES,
        rat::{
            ffrat::FFRat,
            finite::{
                linear::{RecLinear, Unit, UNIT},
                thiele::ThieleRec,
            },
        },
    },
    traits::TryEval,
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
        Self { extra_pts }
    }

    pub fn extra_pts(&self) -> usize {
        self.extra_pts
    }
}

pub(crate) fn nterms_with_max_pows<const N: usize>(
    max_pows: [u32; N],
) -> usize {
    max_pows.iter().map(|&n| n as usize).product()
}

pub(crate) fn gen_poly_with_max_pows<const N: usize>(
    max_pows: [u32; N],
) -> FlatPoly<Unit, N> {
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

#[derive(Error, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FailedRec<const N: usize> {
    #[error("Need points for reconstruction")]
    Empty,
    #[error("Need more points for arguments {:?} in characteristic {modulus}", format_args(*args, *ncoord))]
    MoreAt {
        modulus: u64,
        args: [u64; N],
        ncoord: usize,
    },
    #[error("Need more points in characteristic {modulus}. Estimated total number: {}", nexpected_to_string(*.nexpected))]
    MorePts {
        modulus: u64,
        nexpected: Option<usize>,
    },
    #[error("Need approximately {0} points in new characteristic")]
    MoreMods(usize),
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

fn format_args<const N: usize>(args: [u64; N], ncoord: usize) -> [String; N] {
    let mut res = args.map(|a| a.to_string());
    res[ncoord] = "x".to_string();
    res
}

fn nexpected_to_string(n: Option<usize>) -> String {
    if let Some(n) = n {
        n.to_string()
    } else {
        "unknown".to_string()
    }
}

pub fn rec_linear_from_pts<const N: usize>(
    pts: &mut [(u64, Vec<([u64; N], u64)>)],
    extra_pts: usize,
) -> Result<Rat<FlatPoly<Integer, N>>, FailedRec<N>> {
    use FailedRec::*;
    if pts.is_empty() {
        return Err(Empty);
    }

    let mut num_pows = [0; N];
    let mut den_pows = [0; N];
    for ncoord in 0..N {
        let (modulus, indices) = get_longest_streak(pts, ncoord);
        let streak =
            &pts.iter().find(|(m, _)| *m == modulus).unwrap().1[indices];
        let pts = streak.iter().map(|(z, q_z)| (z[ncoord], *q_z));
        let rec_pow = try_rec_max_num_den_pow(modulus, pts, extra_pts);
        let Some([max_num_pow, max_den_pow]) = rec_pow? else {
            let args = streak[0].0;
            return Err(MoreAt {
                modulus,
                args,
                ncoord,
            });
        };
        num_pows[ncoord] = max_num_pow;
        den_pows[ncoord] = max_den_pow;
    }
    debug!("Numerator powers: {num_pows:?}");
    debug!("Denominator powers: {den_pows:?}");
    let nexpected = nterms_with_max_pows(num_pows)
        + nterms_with_max_pows(den_pows)
        + extra_pts
        - 1;
    let (modulus, most_pts) =
        pts.iter().max_by_key(|(_, pts)| pts.len()).unwrap();
    let modulus = *modulus;
    let Some(mut mod_rec) =
        rec_with_pows(most_pts, modulus, num_pows, den_pows)?
    else {
        debug!("Reconstruction over mod {modulus} failed");
        return Err(MorePts {
            modulus,
            nexpected: Some(nexpected),
        });
    };
    if max_pows_reached(&mod_rec.rat) != [num_pows, den_pows] {
        debug!("Reconstruction over mod {modulus} did not reach powers {num_pows:?}, {den_pows:?}");
        return Err(MorePts {
            modulus,
            nexpected: None,
        });
    }
    let mut res: Result<Rat<FlatPoly<Integer, N>>, _> = (&mod_rec).try_into();
    let first_modulus = modulus;

    // TODO: might be better to iterate over LARGE_PRIMES
    for (m, pts) in pts {
        use ControlFlow::*;
        if *m == first_modulus {
            continue;
        }
        if let Break(res) = add_rec(&mut res, &mut mod_rec, pts, *m, extra_pts)?
        {
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
    extra_pts: usize,
) -> Result<
    ControlFlow<Result<Rat<FlatPoly<Integer, N>>, FailedRec<N>>>,
    UnknownMod,
> {
    seq! { Q in 0..114 {
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
    extra_pts: usize,
) -> ControlFlow<Result<Rat<FlatPoly<Integer, N>>, FailedRec<N>>> {
    use ControlFlow::*;
    // if we can already reproduce the points we are done
    if let Ok(res_ref) = res.as_ref() {
        if pts.len() > extra_pts {
            let sample_same = pts
                .iter()
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
        return Break(Err(FailedRec::MorePts {
            modulus: P,
            nexpected: Some(nexpected),
        }));
    };
    debug!("Reconstructed {next_mod_rec}");
    let [num_pows, den_pows] = max_pows_reached(&rat.rat);
    if max_pows_reached(&next_mod_rec) != [num_pows, den_pows] {
        debug!("Reconstruction over mod {P} did not reach powers {num_pows:?}, {den_pows:?}");
        return Break(Err(FailedRec::MorePts {
            modulus: P,
            nexpected: None,
        }));
    }
    rat.merge_crt(next_mod_rec);
    *res = (&*rat).try_into();
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
    seq! { Q in 0..114 {
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
    ncoord: usize,
) -> (u64, Range<usize>) {
    let mut longest_streak = (0, 0..0);
    for (modulus, pts) in pts {
        pts.sort_unstable_by_key(|(coord, _)| {
            let mut coord = *coord;
            coord[ncoord] = 0;
            coord
        });
        let streak = longest_streak_pos_by(pts, |a, b| {
            let mut a_coord = a.0;
            a_coord[ncoord] = 0;
            let mut b_coord = b.0;
            b_coord[ncoord] = 0;
            a_coord == b_coord
        });
        if streak.len() > longest_streak.1.len() {
            longest_streak = (*modulus, streak);
        }
    }
    longest_streak
}

fn try_rec_max_num_den_pow(
    modulus: u64,
    pts: impl IntoIterator<Item = (u64, u64)>,
    extra_pts: usize,
) -> Result<Option<[u32; 2]>, UnknownMod> {
    seq! { Q in 0..114 {
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
    extra_pts: usize,
) -> Option<[u32; 2]> {
    let rec = ThieleRec::new(extra_pts).rec_from_seq(pts)?;
    let rec: Rat<DensePoly<_>> = rec.into();
    Some([rec.num().len() as u32, rec.den().len() as u32])
}

fn longest_streak_pos_by<T, F>(slice: &[T], mut is_eq: F) -> Range<usize>
where
    F: FnMut(&T, &T) -> bool,
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
