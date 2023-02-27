use std::{ops::{AddAssign, ControlFlow}, num::NonZeroUsize};

use galois_fields::Z64;
use log::{debug, trace};
use rand::Rng;

use crate::{sparse_poly::{SparsePoly, SparseMono}, traits::{Eval, Rec, WithVars, Zero, One, Shift}, rat::Rat, rec_thiele::ThieleRec, dense_poly::{DensePoly, DensePoly1, DensePoly2}, rec_newton::NewtonPolyRec, rec_linear::LinearRec};


/// Rational function reconstruction
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct RatRec {
    extra_pts: usize
}

impl RatRec {
    fn new(extra_pts: usize) -> Self {
        Self { extra_pts }
    }
}

const N: usize = 2;
const M: usize = 1;

// transform from [x0, x1, ..., t] to
// [z0, z1, z2, ...] = [t*x0 + s0, t*x1 + s1, t*x2 + s2, ..., t + sM]
// with a `shift` = [s0, s1, s2, ..., sM]
fn x_to_z<const N: usize, const P: u64>(
    mut coord: [Z64<P>; N],
    shift: &[Z64<P>; N]
) -> [Z64<P>; N] {
    let (t, xs) = coord.split_last_mut().unwrap();
    for x in xs {
        *x *= *t
    }
    for (z, s) in coord.iter_mut().zip(shift.iter()) {
        *z += s;
    }
    coord
}

// iterate over coordinates [x0, x1, ..., t] with varying values of t
// transforming to z
// with a shift s = [s0, s1, ...]
// returns an iterator over tuples (t, z)
fn transformed_coord_iter<const N: usize, const P: u64>(
    mut xcoord: [Z64<P>; N],
    shift: [Z64<P>; N]
) -> impl Iterator<Item = (Z64<P>, [Z64<P>; N])> {
    std::iter::repeat_with(move || {
        xcoord[N - 1] += Z64::one();
        (xcoord[N - 1], x_to_z(xcoord, &shift))
    }).take(P as usize)
}

// find a shift such that the shifted function fs(z) = f(z + s)
// has a constant leading term in the denominator
// In other words, find s such that f(s) is well-defined
fn find_shift<F, const P: u64, const N: usize>(
    mut f: F,
    mut rng: impl Rng,
) -> [Z64<P>; N]
where F: FnMut([Z64<P>; N]) -> Option<Z64<P>>
{
    let mut shift = [Z64::zero(); N];
    for i in (0..N).cycle() {
        if f(shift).is_some() {
            return shift;
        }
        shift[i] = rng.gen()
    }
    unreachable!()
}

impl<F, const P: u64> Rec<RatRec, [Z64<P>; N]> for F
where F: FnMut([Z64<P>; N]) -> Option<Z64<P>> {
    type Output = Option<Rat<SparsePoly<Z64<P>, N>>>;

    fn rec_with_ran(
        &mut self,
        rec: RatRec,
        mut rng: impl Rng
    ) -> Self::Output {
        let shift = find_shift(|z| (self)(z), &mut rng);
        debug!("Shift: {shift:?}");

        // initial [x0, x1, ..., t]
        let base_coord: [Z64<P>; N] = rng.gen();
        debug!("Starting multivariate rational reconstruction with base point {base_coord:?}");

        // varying t while keeping x0, ... fixed
        let coord = transformed_coord_iter(base_coord, shift);

        // reconstruct in t
        let thiele_rec = ThieleRec::new(rec.extra_pts);
        let rat = thiele_rec.rec_from_seq(
            coord.filter_map(|(t, z)| (self)(z).map(|v| (t, v)))
        )?;
        let rat: Rat<_> = rat.into();
        debug!("From Thiele reconstruction: {}", rat.with_vars(&["t"]));
        debug_assert!(!rat.den().coeff(0).is_zero());

        // first point to be used for starting the reconstruction
        let pt = RatPt::<P, M>::new([Z64::zero(); M], rat);

        let mut rec = RecHelper::new(pt, rec.extra_pts, base_coord, shift);
        use NumOrDen::*;
        let num_rec = rec.rec(|x| (self)(x), Num)?;
        debug!("Reconstructed numerator: {num_rec}");
        let den_rec = rec.rec(|x| (self)(x), Den)?;
        debug!("Reconstructed denominator: {den_rec}");
        Some(Rat::from_num_den_unchecked(num_rec, den_rec))
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
enum NumOrDen {
    Num,
    Den
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct RecHelper<const P: u64> {
    known_pts: Vec<RatPt<P, M>>,
    extra_pts: usize,
    base_coord: [Z64<P>; N],
    shift: [Z64<P>; N],
}

impl<const P: u64> RecHelper<P> {
    fn new(
        first_pt: RatPt<P, M>,
        extra_pts: usize,
        base_coord: [Z64<P>; N],
        shift: [Z64<P>; N],
    ) -> Self {
        Self {
            known_pts: vec![first_pt],
            extra_pts,
            base_coord,
            shift,
        }
    }

    fn rec<F>(
        &mut self,
        mut f: F,
        wot: NumOrDen
    ) -> Option<SparsePoly<Z64<P>, N>>
    where F: FnMut([Z64<P>; N]) -> Option<Z64<P>>
    {
        debug!("Reconstruct {wot:?} coefficients");
        use NumOrDen::*;
        let base_x_coord: [Z64<P>; M] = self.base_coord[..M].try_into().unwrap();

        let nnum = self.known_pts[0].num_coeff.len();
        let nden = self.known_pts[0].den_coeff.len();
        let den_terms = NonZeroUsize::new(1 + nden).unwrap();
        // on to the meat:
        // iteratively reconstruct the highest-degree polynomial in the numerator / denominator
        // storing also the results for all other coefficients for later use
        let mut rec: SparsePoly<_, N> = Zero::zero();
        let max_idx = match wot {
            Num => nnum,
            Den => usize::from(nden)
        };
        for coeff_idx in (0..max_idx).rev() {
            debug!("Reconstruct coefficient of t^{coeff_idx}");
            let mut poly = NewtonPolyRec::<P>::new(self.extra_pts);
            // what do we know at this point?
            // TODO: maybe use a BTreeMap so we don't have to sort?
            self.known_pts.sort_unstable();
            let known = self.known_pts.iter();
            let known: Vec<_> = match wot {
                    Num => known.map(|pt| (pt.offset, pt.num_coeff[coeff_idx])).collect(),
                    Den => known.map(|pt| (pt.offset, pt.den_coeff[coeff_idx])).collect(),
            };
            let mut known = known.iter().peekable();

            let first_pt = known.next().unwrap();
            let mut next_offset = first_pt.0;
            let mut next_val = first_pt.1;
            loop {
                let coord = add(base_x_coord, next_offset);
                match poly.add_pt(&coord, next_val) {
                    ControlFlow::Continue(n) => {
                        next_offset[n] += Z64::one();
                        next_offset[(n + 1)..].fill(Z64::zero());
                        // skip known points with a lesser offset
                        // they are not needed for this reconstruction
                        // TODO: some way to do this with `skip_while`?
                        while let Some((offset, _)) = known.peek() {
                            if offset < &next_offset {
                                known.next();
                            } else {
                                break;
                            }
                        }
                        // TODO: accept points with larger offset in the same coordinate
                        match known.peek() {
                            Some((offset, val)) if offset == &next_offset => {
                                // TODO: apply shift corrections
                                trace!("Using known value at offset {next_offset:?}: {val}");
                                next_val = *val;
                            },
                            _ => {
                                trace!("Need new value at offset {next_offset:?}");
                                let mut coord = self.base_coord;
                                for (x, s) in coord[..M].iter_mut().zip(next_offset.iter()) {
                                    *x += s;
                                }
                                let eval_pts = transformed_coord_iter(coord, self.shift);

                                // TODO: don't reconstruct the already known coefficients
                                let linear_rec = LinearRec::new(nnum, den_terms);
                                let rat = linear_rec.rec_from_seq(
                                    eval_pts.filter_map(|(t, z)| f(z).map(|v| (t, v)))
                                )?;
                                trace!("From rational reconstruction in t: {}", rat.with_vars(&["t"]));
                                next_val = match wot {
                                    Num => *rat.num().coeff(coeff_idx),
                                    Den => *rat.den().coeff(coeff_idx),
                                };
                                let pt = RatPt::<P, M>::new(next_offset, rat);

                                self.known_pts.push(pt);
                            },
                        }

                    },
                    ControlFlow::Break(()) => {
                        let rec_poly: DensePoly1<Z64<P>> = poly.into_poly().into();
                        // turn it back from a polynomial in [x0, x1, ..., t] to [z0, z1, ...]
                        // first expand it out
                        let rec_poly: SparsePoly<Z64<P>, M> = rec_poly.into();
                        // rescale, adding back the power in t
                        let terms = rec_poly.into_terms().into_iter()
                            .map(|c| {
                                let mut powers = [0; N];
                                powers[0..M].copy_from_slice(&c.powers);
                                powers[M] = coeff_idx as u32 - c.powers.iter().sum::<u32>();
                                SparseMono::new(c.coeff, powers)
                            });
                        // // subtract effect of shift in this polynomial from
                        // // previously computed function values
                        // let corr = rec_poly.clone().shift(shift) - &rec_poly;
                        // debug!("New correction term: {corr}");

                        let rec_poly = SparsePoly::from_terms(terms.collect());
                        debug!("Reconstructed numerator contribution {rec_poly}");
                        rec += rec_poly;
                        break;
                    },
                }
            }
        }
        Some(rec)
    }
}



fn add<T: AddAssign, const N: usize>(mut a: [T; N], b: [T; N]) -> [T; N] {
    for (a, b) in a.iter_mut().zip(b.into_iter()) {
        *a += b;
    }
    a
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct RatPt<const P: u64, const N: usize> {
    offset: [Z64<P>; N],
    num_coeff: Vec<Z64<P>>,
    den_coeff: Vec<Z64<P>>,
}

impl<const P: u64, const N: usize> RatPt<P, N> {
    fn new(pt: [Z64<P>; N], rat: Rat<DensePoly<Z64<P>>>) -> Self {
        let (num, den) = rat.into_num_den();
        Self {
            offset: pt,
            num_coeff: num.into_coeff(),
            den_coeff: den.into_coeff(),
        }
    }
}


#[cfg(test)]
mod tests {
    use std::iter::repeat_with;

    use rand::Rng;
    use rand_xoshiro::rand_core::SeedableRng;

    use crate::{traits::{TryEval, One, Zero}, sparse_poly::SparseMono};

    use super::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn rec_rat2_small_noshift() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_POW: u32 = 2;
        const MAX_COEFF: u64 = 29;
        const P: u64 = 1152921504606846883;
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = RatRec::new(1);

        for _ in 0..NTESTS {
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let coeff = repeat_with(
                || SparseMono::new(
                    Z64::<P>::new_unchecked(rng.gen_range(0..MAX_COEFF)),
                    [(); 2].map(|_| rng.gen_range(0..=MAX_POW))
                )
            ).take(nterms).collect();
            let num = SparsePoly::from_terms(coeff);

            let den = if num.is_zero() {
                One::one()
            } else {
                let mut den: SparsePoly<_, 2> = Zero::zero();
                while den.is_zero() {
                    let coeff = repeat_with(
                        || SparseMono::new(
                            Z64::new_unchecked(rng.gen_range(0..MAX_COEFF)),
                            [(); 2].map(|_| rng.gen_range(0..=MAX_POW))
                        )
                    ).take(nterms).collect();
                    den = SparsePoly::from_terms(coeff);
                }
                // ensure we don't need a shift
                let mut den = den.into_terms();
                den[0].powers = Zero::zero();
                den[0].coeff = One::one();
                SparsePoly::from_terms(den)
            };

            let rat = Rat::from_num_den_unchecked(num, den);
            eprintln!("trying to reconstruct {rat}");

            let res = (|z| rat.try_eval(&z))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            eprintln!("reconstructed {res}");
            assert_eq!(rat, res);
        }
    }
}
