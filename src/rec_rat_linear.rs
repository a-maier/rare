use std::ops::ControlFlow;

use ffnt::Z64;
use log::{debug, warn, trace};
use paste::paste;
use rand::Rng;
use rug::Integer;
use seq_macro::seq;

use crate::{
    rat::{NoneError, Rat},
    rec_linear::{RecLinear, UNIT, Unit},
    rec_rat::{combine_crt_rat, FFRat, LARGE_PRIMES},
    sparse_poly::{SparsePoly, SparseMono},
    traits::{TryEval, One, Zero, Rec}, rec_thiele::ThieleRec, dense_poly::DensePoly,
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

seq! {N in 2..=16 {
    paste! {
        #[derive(Clone, Debug, PartialEq, PartialOrd)]
        pub struct [<RatRecLinear N>]<const P: u64> {
            extra_pts: usize,
            degree_rec: ThieleRec<P>,
            last_z: Option<[Z64<P>; N]>,
            powers: [[u32; 2]; N],
            next_needed_pow: usize,
            rat: FFRat<N>,
            res: Result<Rat<SparsePoly<Integer, N>>, NoneError>,
            done: bool,
        }

        impl<const P: u64> [<RatRecLinear N>]<P> {

            pub fn new(extra_pts: usize) -> Self {
                Self {
                    extra_pts,
                    powers: [[0; 2]; N],
                    degree_rec: ThieleRec::new(extra_pts),
                    last_z: None,
                    next_needed_pow: 0,
                    rat: Default::default(),
                    res: Err(NoneError {  }),
                    done: false,
                }
            }

            pub fn add_pt(
                &mut self,
                z: [Z64<P>; N],
                q_z: Z64<P>,
            ) -> ControlFlow<(), Needed<P, N>> {
                use ControlFlow::*;
                use Needed::Pt;
                if self.status() != ReconstructionStatus::Degrees {
                    return self.ask_for_new_mod();
                }
                if let Some(last_z) = self.last_z.as_ref() {
                    let shift_ok = z.iter().zip(last_z)
                        .enumerate()
                        .all(|(n, (z, zp))| n == self.next_needed_pow || z == zp);
                    if !shift_ok {
                        warn!("Ignoring point {z:?}: should differ from {last_z:?} only in coordinate {}", self.next_needed_pow);
                        return Continue(Pt(self.request_next_pt()));
                    }
                }
                self.last_z = Some(z);
                match self.degree_rec.add_pt(z[self.next_needed_pow], q_z) {
                    Continue(()) => Continue(Pt(self.request_next_pt())),
                    Break(()) => {
                        let rat = std::mem::replace(
                            &mut self.degree_rec,
                            ThieleRec::new(self.extra_pts)
                        ).into_rat();
                        let rat: Rat<DensePoly<_>> = rat.into();
                        let num_pow = rat.num().len().try_into().unwrap();
                        let den_pow = rat.den().len().try_into().unwrap();
                        self.powers[self.next_needed_pow] = [num_pow, den_pow];
                        debug!(
                            "Powers in variable {}: {:?}",
                            self.next_needed_pow,
                            self.powers[self.next_needed_pow]
                        );
                        if num_pow == 0 {
                            self.res = Ok(Zero::zero());
                            self.done = true;
                            return Break(());
                        }
                        self.next_needed_pow += 1;
                        if self.next_needed_pow < self.powers.len() {
                            Continue(Pt(self.request_next_pt()))
                        } else {
                            self.ask_for_new_mod()
                        }
                    }
                }
            }

            pub fn add_pts<'a, const Q: u64>(
                &mut self,
                pts: &[([Z64<Q>; N], Z64<Q>)]
            ) -> ControlFlow<(), Needed<P, N>> {
                use ControlFlow::*;
                use Needed::Pt;
                use ReconstructionStatus::*;
                match self.status() {
                    Degrees => {
                        if Q == P {
                            // cast to the correct type
                            // SAFETY: we just checked that the type is actually &'a [([Z64<P>; N], Z64<P>)]
                            //         the lifetime is explicitly set to match the original one
                            let pts: &'a [_] = unsafe {
                                std::slice::from_raw_parts(pts.as_ptr() as _, pts.len())
                            };
                            let mut res = Break(());
                            for (z, q_z) in pts.iter().copied() {
                                res = self.add_pt(z, q_z);
                            }
                            res
                        } else {
                            warn!("Ignoring points in characteristic {Q}, need characteristic {P}");
                            Continue(Pt(self.request_next_pt()))
                        }
                    },
                    Rat => if self.rat.modulus.is_zero() {
                        self.rec_first_rat_mod(pts)
                    } else {
                        self.rec_rat_mod(pts)
                    },
                    Done => Break(()),
                }
            }

            fn rec_first_rat_mod<const Q: u64>(
                &mut self,
                pts: &[([Z64<Q>; N], Z64<Q>)],
            ) -> ControlFlow<(), Needed<P, N>> {
                let num = gen_poly_with_max_pows(self.num_pows());
                let den = gen_poly_with_max_pows(self.den_pows());
                let ansatz = Rat::from_num_den_unchecked(num, den);
                trace!("ansatz: {ansatz:#?}");
                if let Some(next_mod_rec) = ansatz.rec_linear(pts.iter().copied()) {
                    debug!("Finished reconstruction modulo {Q}: {next_mod_rec}");
                    self.rat = FFRat::from(next_mod_rec);
                };
                self.ask_for_new_mod()
            }

            fn num_pows(&self) -> [u32; N] {
                self.powers.map(|d| d[0])
            }

            fn den_pows(&self) -> [u32; N] {
                self.powers.map(|d| d[1])
            }

            fn rec_rat_mod<const Q: u64>(
                &mut self,
                pts: &[([Z64<Q>; N], Z64<Q>)],
            ) -> ControlFlow<(), Needed<P, N>> {
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
                debug!("Finished reconstruction modulo {Q}: {next_mod_rec}");
                let mod_rec = std::mem::take(&mut self.rat);
                self.rat = combine_crt_rat(mod_rec, next_mod_rec);
                self.res  = (&self.rat).try_into();
                self.ask_for_new_mod()
            }

            pub fn status(&self) -> ReconstructionStatus {
                use ReconstructionStatus::*;
                if self.next_needed_pow < self.powers.len() {
                    Degrees
                } else if self.done {
                    Done
                } else {
                    Rat
                }
            }

            fn ask_for_new_mod(&self) -> ControlFlow<(), Needed<P, N>> {
                use ControlFlow::*;
                if self.rat.modulus.is_zero() {
                    let nnum = nterms_with_max_pows(self.num_pows());
                    let nden = nterms_with_max_pows(self.den_pows());
                    return Continue(Needed::Any(nnum + nden + self.extra_pts - 1));
                }
                let ncoeff = self.rat.rat.num().len() + self.rat.rat.den().len() - 1;
                if ncoeff == 0 {
                    Break(())
                } else {
                    Continue(Needed::Any(ncoeff))
                }
            }

            pub fn request_next_pt(&mut self) -> [Z64<P>; N] {
                let z = self.last_z.as_mut().unwrap();
                z[self.next_needed_pow] += Z64::one();
                *z
            }

            pub fn into_rat(self) -> Option<Rat<SparsePoly<Integer, N>>> {
                self.res.ok()
            }

            pub fn extra_pts(&self) -> usize {
                self.extra_pts
            }
        }
    }
}}

fn nterms_with_max_pows<const N: usize>(
    max_pows: [u32; N]
) -> usize {
    max_pows.iter().map(|&n| n as usize).product()
}

fn gen_poly_with_max_pows<const N: usize>(max_pows: [u32; N]) -> SparsePoly<Unit, N> {
    let num_terms = nterms_with_max_pows(max_pows);
    let mut terms = Vec::with_capacity(num_terms);
    for mut i in 0..(num_terms as u32) {
        let mut pows = [0u32; N];
        for (pow, max) in pows.iter_mut().zip(&max_pows) {
            *pow = i % *max;
            i /= max;
        }
        terms.push(SparseMono::new(UNIT, pows));
    }
    // TODO: would be better to generate them in the
    // correct order instead of sorting
    terms.sort_unstable();
    SparsePoly::from_raw_terms(terms)
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
    type Output = Option<Rat<SparsePoly<Integer, N>>>;

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
        let mut res: Result<Rat<SparsePoly<Integer, N>>, _> =
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

#[cfg(test)]
mod tests {
    use super::*;
    use ::rand::{Rng, SeedableRng};

    use crate::rec_rat::LARGE_PRIMES;
    use crate::sparse_poly::SparseMono;
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

    fn rand_term<const N: usize>(mut rng: impl Rng) -> SparseMono<Integer, N> {
        let pow = [(); N].map(|_| rng.gen_range(0..=MAX_POW));
        SparseMono::new(rand_int(rng), pow)
    }

    fn rand_poly<const N: usize>(mut rng: impl Rng) -> SparsePoly<Integer, N> {
        let nterms = rng.gen_range(0..=MAX_TERMS);
        SparsePoly::from_terms(
            (0..nterms).map(|_| rand_term(&mut rng)).collect(),
        )
    }

    fn rand_rat<const N: usize>(
        mut rng: impl Rng,
    ) -> Rat<SparsePoly<Integer, N>> {
        let mut den = SparsePoly::zero();
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

    seq!(NVARS in 2..=3 {
        paste! {
            #[test]
            fn [<rec_rat_explicit_ NVARS>]() {
                log_init();

                let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
                for _ in 0..NTESTS {
                    const P: u64 = LARGE_PRIMES[0];

                    let orig = rand_rat::<NVARS>(&mut rng);
                    eprintln!("trying to reconstruct {orig}");
                    let mut rec = [<RatRecLinear NVARS>]::new(1);
                    let (z, q_z) = std::iter::repeat_with(|| {
                        let z: [Z64<P>; NVARS] = [(); NVARS].map(|_| rng.gen());
                        orig.try_eval(&z).map(|q_z| (z, q_z))
                    })
                        .flatten()
                        .next()
                        .unwrap();
                    let mut pts = vec![(z, q_z)];
                    loop {
                        use ControlFlow::*;
                        use Needed::*;
                        let (z, q_z) = pts[pts.len() - 1];
                        match rec.add_pt(z, q_z) {
                            Continue(Pt(mut z)) => {
                                let q_z = loop {
                                    if let Some(q_z) = orig.try_eval(&z) {
                                        break q_z;
                                    }
                                    z = rec.request_next_pt();
                                };
                                pts.push((z, q_z));
                            },
                            Continue(Any(mut n)) => {
                                // TODO: reuse points
                                seq!{ N in 0..20 {{
                                    const P: u64 = LARGE_PRIMES[N];
                                    let pts = Vec::from_iter(
                                        std::iter::repeat_with(|| {
                                            let z: [Z64<P>; NVARS] = [(); NVARS].map(|_| rng.gen());
                                            orig.try_eval(&z).map(|q_z| (z, q_z))
                                        })
                                            .flatten()
                                            .take(n)
                                    );
                                    match rec.add_pts(&pts) {
                                        Continue(Any(nnew)) => n = nnew,
                                        Break(()) => break,
                                        _ => panic!("Unexpected reconstruction return value")
                                    }
                                }}}
                                let _n = n;
                                panic!("Need more than 20 characteristics!");
                            },
                            Break(()) => break
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
