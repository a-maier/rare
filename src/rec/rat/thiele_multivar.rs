use std::{any::Any, ops::ControlFlow::{self, *}};

use log::{debug, trace};
use rand::Rng;
use thiserror::Error;

use crate::{
    algebra::{poly::{dense::DensePoly, flat::{FlatMono, FlatPoly}}, rat::{NoneError, Rat}},
    rec::rat::{
        ffrat::FFRat, finite::{degree_rec, thiele::{ThieleRat, ThieleRec}}, sampler::{self, Sampler}, util
    }, traits::{One, Zero}, Integer, Z64
};

pub type IntRat<const N: usize> = Rat<FlatPoly<Integer, N>>;

#[derive(Debug)]
pub struct Rec<const P: u64, const N: usize> {
    extra_pts: usize,
    rec: DegreeOrScaledRec<P, N>,
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct DegreeRec<const P: u64, const N: usize> {
    rec: degree_rec::DegreeRec<P, N>,
    shift: [u64; N],
}

#[derive(Debug)]
struct ScaledRec<const N: usize> {
    rec: Box<dyn Any>,
    modulus: u64,
    rat: FFRat<N>,
    res: Result<Rat<FlatPoly<Integer, N>>, NoneError>,
    transform: Transform<N>,
    sampler: Sampler,
}

impl<const N: usize> Default for ScaledRec<N> {
    fn default() -> Self {
        Self {
            rec: Box::new(()),
            modulus: 0,
            rat: Default::default(),
            res: Err(NoneError {}),
            transform: Default::default(),
            sampler: Default::default(),
        }
    }
}

#[derive(Debug)]
enum DegreeOrScaledRec<const P: u64, const N: usize> {
    DegreeRec(DegreeRec<P, N>),
    ScaledRec(ScaledRec<N>),
}

impl<const P: u64, const N: usize> Rec<P, N> {
    pub fn with_shift(extra_pts: usize, shift: [u64; N]) -> Self {
        let rec = if N == 1 {
            DegreeOrScaledRec::ScaledRec(ScaledRec::default())
        } else {
            DegreeOrScaledRec::DegreeRec(DegreeRec{
                rec: degree_rec::DegreeRec::new(extra_pts),
                shift,
            })
        };
        Self {
            extra_pts,
            rec,
        }
    }

    pub fn with_random_shift(extra_pts: usize, mut rng: impl Rng) -> Self {
        Self::with_shift(extra_pts, [(); N].map(|_| rng.gen()))
    }

    pub fn add_pt<const Q: u64>(
        &mut self,
        z: [Z64<Q>; N],
        q_z: Z64<Q>
    ) -> Result<ControlFlow<IntRat<N>, Status<N>>, Error<P, N>> {
        use Status::*;
        trace!("Adding: f({z:?}) = {q_z}");
        match self.rec {
            DegreeOrScaledRec::DegreeRec(ref mut rec) => {
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
                // we don't actually need the power of the last
                // variable for anything, so we stop before
                // reconstructing it. We set the power to `u32::MAX`,
                // since that simplifies `Transform::undo`
                let res = match rec.rec.add_pt(z, q_z)?  {
                    Continue(n) if n < (N - 1) => Varying(n),
                    _ => {
                        let mut powers = rec.rec.cur_powers();
                        powers[N - 1] = [u32::MAX; 2];
                        let mut scaled_rec = ScaledRec::default();
                        scaled_rec.transform = Transform {
                            scale: powers.map(|p| p.into_iter().max().unwrap()),
                            shift: rec.shift,
                        };
                        self.rec = DegreeOrScaledRec::ScaledRec(scaled_rec);
                        Scaling
                    }
                };
                Ok(Continue(res))
            },
            DegreeOrScaledRec::ScaledRec(ref mut rec) => {
                if rec.modulus == 0 {
                    debug!("New modulus: {Q}");
                    rec.rec = Box::new(ThieleRec::<Q>::new(self.extra_pts));
                    if rec.sampler.status() == sampler::Status::Failed {
                        rec.sampler = Sampler::new(self.extra_pts);
                    }
                    rec.modulus = Q;
                } else if rec.modulus != Q {
                    return Err(util::RecError::Mod{
                        expected: rec.modulus,
                        found: Q,
                    }.into());
                }
                rec.transform.check_pt(z)?;
                if let Ok(res) = rec.res.as_mut() {
                    if rec.sampler.add_pt(&z, q_z, res) == sampler::Status::Complete {
                        return Ok(Break(std::mem::take(res)));
                    }
                }
                let scaled_rec = rec.rec.downcast_mut::<ThieleRec<Q>>().unwrap();
                let res = match scaled_rec.add_pt(z[0], q_z) {
                    Continue(()) => Scaling,
                    Break(()) => {
                        let rat = std::mem::take(scaled_rec).into_rat();
                        rec.modulus = 0;
                        let rat = rec.transform.undo(rat);
                        debug!("Multivariate reconstruction result: {rat}");
                        if rec.rat.modulus == 0 {
                            rec.rat = rat.into();
                        } else {
                            rec.rat.merge_crt(rat);
                        }
                        rec.res = (&rec.rat).try_into();
                        NextMod
                    },
                };
                Ok(Continue(res))
            }
        }
    }

    pub fn to_args<const Q: u64>(&self, z: Z64<Q>) -> Option<[Z64<Q>; N]> {
        if let DegreeOrScaledRec::ScaledRec(ref rec) = self.rec {
            Some(rec.transform.transform(z))
        } else {
            None
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct Transform<const N: usize> {
    scale: [u32; N],
    shift: [u64; N],
}

impl<const N: usize> Transform<N> {
    fn transform<const P: u64>(&self, z: Z64<P>) -> [Z64<P>; N] {
        let mut res = [z; N];
        for n in 1..N {
            res[n] = res[n - 1].powu(self.scale[n - 1] as u64);
        }
        for (lhs, shift) in res.iter_mut().skip(1).zip(self.shift) {
            *lhs += Z64::from(shift);
        }
        res
    }

    fn check_pt<const P: u64, const Q: u64>(
        &self,
        z: [Z64<P>; N],
    ) -> Result<(), Error<Q, N>> {
        let expected = self.transform(z[0]);
        if z != expected {
            Err(Error::Scaling {
                expected: expected.map(|z| u64::from(z)),
                got: z.map(|z| u64::from(z))
            })
        } else {
            Ok(())
        }
    }

    fn undo<const P: u64>(
        &self,
        rat: ThieleRat<Z64<P>>
    ) -> Rat<FlatPoly<Z64<P>, N>> {
        debug!("Undo scaling");
        let rat: Rat<DensePoly<_>> = rat.into();
        debug!("Expanded: {rat}");
        if rat.is_zero() {
            return Rat::zero();
        }
        let (num, den) = rat.into_num_den();
        let mut num = self.undo_poly(num);
        let mut den = self.undo_poly(den);
        let norm = num.terms().last().unwrap().coeff.inv();
        num *= norm;
        den *= norm;
        Rat::from_num_den_unchecked(num, den)
    }

    fn undo_poly<const P: u64>(
        &self,
        poly: DensePoly<Z64<P>>
    ) -> FlatPoly<Z64<P>, N> {
        if N == 1 {
            return FlatPoly::from_raw_terms(
                poly.into_coeff()
                    .into_iter()
                    .enumerate()
                    .filter_map(|(pow, coeff)| if coeff.is_zero() {
                        None
                    } else {
                        Some(FlatMono{coeff, powers: [pow as u32; N]})
                    })
                    .collect()
            );
        }
        assert_eq!(self.scale[N - 1], u32::MAX);
        let shift = self.shift.map(|s| s.into());
        let mut res = FlatPoly::new();
        for (pow, coeff) in poly.into_coeff().into_iter().enumerate() {
            let mut pow = pow as u32;
            let mut expanded = FlatPoly::new();
            for (n, s) in self.scale.iter().enumerate() {
                let mut res_powers = [0; N];
                let zn_pow = pow % s;
                pow /= s;
                if n == 0 {
                    res_powers[n] = zn_pow;
                    expanded = FlatMono {
                        coeff: Z64::one(),
                        powers: res_powers,
                    }.into();
                } else {
                    res_powers[n] = 1;
                    let base = FlatMono {
                        coeff: Z64::one(),
                        powers: res_powers,
                    } - FlatMono {
                        coeff: shift[n - 1],
                        powers: [0; N],
                    };
                    let tmp = &expanded * &base.powu(zn_pow);
                    expanded = tmp;
                }
            }
            let expanded: FlatPoly<_, N> = expanded.into();
            res += expanded * &coeff;
        }
        res
    }
}

impl<const N: usize> Default for Transform<N> {
    fn default() -> Self {
        Self{
            scale: [0; N],
            shift: [0; N],
        }
    }
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
    use log::debug;
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
        while den.is_zero() {
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
                    let mut rec: Rec<P, NVARS> = Rec::with_random_shift(1, &mut rng);
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
                                    let mut z = rec.to_args(z0).unwrap();
                                    let q_z = loop {
                                        if let Some(val) = orig.try_eval(&z) {
                                            break val;
                                        }
                                        z0 += Z64::one();
                                        z = rec.to_args(z0).unwrap();
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
