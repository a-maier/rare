use ffnt::Z64;
use log::{debug, trace, warn};
use rug::Integer;

use crate::{
    algebra::{
        poly::flat::{FlatMono, FlatPoly},
        rat::{NoneError, Rat},
    },
    rec::rat::ffrat::FFRat,
    rec::{
        probe::Probe,
        rat::finite::linear::{RecLinear, Unit, UNIT},
    },
    traits::{TryEval, Zero},
};

use super::finite::degree_rec::DegreeRec;

/// Reconstruction status
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Status<const P: u64, const N: usize> {
    /// Reconstructing numerator and denominator degrees
    ///
    /// We are currently reconstructing the highest powers in the
    /// variable with the given number. The next point should
    /// only differ in this variable from the previous one.
    Degree(usize),
    /// Reconstructing the rational function.
    ///
    /// The argument is the estimated number of probes that will be needed
    /// for reconstruction over a new characteristic. Characteristics
    /// over which reconstruction has already succeeded add no new
    /// information and should be avoided.
    Rat(usize),
    /// Reconstruction has finished .
    ///
    /// The result can now be extracted with `into_rat()`
    Done,
}

/// Rational reconstruction using a combination of Thiele reconstruction and linear systems of equations
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Rec<const P: u64, const N: usize> {
    extra_pts: usize,
    degree_rec: DegreeRec<P, N>,
    powers: [[u32; 2]; N],
    rat: FFRat<N>,
    res: Result<Rat<FlatPoly<Integer, N>>, NoneError>,
    status: Status<P, N>,
}

impl<const P: u64, const N: usize> Rec<P, N> {
    pub fn new(extra_pts: usize) -> Self {
        Self {
            extra_pts,
            powers: [[0; 2]; N],
            degree_rec: DegreeRec::new(extra_pts),
            rat: Default::default(),
            res: Err(NoneError {}),
            status: Status::Degree(0),
        }
    }

    pub fn add_pt(&mut self, z: [Z64<P>; N], q_z: Z64<P>) -> Status<P, N> {
        use std::ops::ControlFlow::{Break, Continue};
        use Status::Degree;
        if !matches!(self.status, Degree(_)) {
            warn!("Not reconstructing degrees, ignoring single point");
            return self.status();
        }
        match self.degree_rec.add_pt(z, q_z) {
            Ok(Break(res)) => {
                debug!("Finished reconstructing variable powers: {res:?}");
                self.powers = res;
                if self.powers.is_zero() {
                    self.res = Ok(Zero::zero());
                    self.status = Status::Done;
                } else {
                    let nnum = nterms_with_max_pows(self.num_pows());
                    let nden = nterms_with_max_pows(self.den_pows());
                    self.status = Status::Rat(nnum + nden + self.extra_pts - 1);
                }
            }
            Ok(Continue(n)) => {
                self.status = Degree(n);
            }
            Err(err) => {
                // TODO: return error
                warn!("{err}");
            }
        }
        self.status()
    }

    pub fn add_pts<'a, const Q: u64>(
        &mut self,
        pts: &'a [Probe<Q, N>],
    ) -> Status<P, N> {
        match self.status() {
            Status::Degree { .. } => {
                if Q == P {
                    // cast to the correct type
                    // SAFETY: we just checked that the type is actually &'a [Probe<P, N>]
                    //         the lifetime is explicitly set to match the original one
                    let pts: &'a [_] = unsafe {
                        std::slice::from_raw_parts(pts.as_ptr() as _, pts.len())
                    };
                    for Probe { arg, val } in pts.iter().copied() {
                        let status = self.add_pt(arg, val);
                        if !matches!(status, Status::Degree { .. }) {
                            break;
                        }
                    }
                } else {
                    warn!("Ignoring points in characteristic {Q}, need characteristic {P}");
                }
            }
            Status::Rat(_) => {
                if self.rat.modulus.is_zero() {
                    self.rec_first_rat_mod(pts)
                } else {
                    self.rec_rat_mod(pts)
                }
            }
            Status::Done => {}
        }
        self.status()
    }

    fn rec_first_rat_mod<const Q: u64>(&mut self, pts: &[Probe<Q, N>]) {
        let num = gen_poly_with_max_pows(self.num_pows());
        let den = gen_poly_with_max_pows(self.den_pows());
        let ansatz = Rat::from_num_den_unchecked(num, den);
        trace!("ansatz: {ansatz:#?}");
        let pts = pts.iter().map(|Probe { arg, val }| (*arg, *val));
        if let Some(next_mod_rec) = ansatz.rec_linear(pts) {
            debug!("Finished reconstruction modulo {Q}: {next_mod_rec}");
            self.rat = FFRat::from(next_mod_rec);
            let ncoeff =
                self.rat.rat.num().len() + self.rat.rat.den().len() - 1;
            self.status = Status::Rat(ncoeff);
        };
    }

    fn rec_rat_mod<const Q: u64>(&mut self, pts: &[Probe<Q, N>]) {
        // if we can already reproduce the points we are done
        if let Ok(res) = self.res.as_ref() {
            if pts.len() >= self.extra_pts {
                let sample_same = pts
                    .iter()
                    .all(|Probe { arg, val }| res.try_eval(arg) == Some(*val));

                if sample_same {
                    self.status = Status::Done;
                    return;
                }
            }
        }

        let pts = pts.iter().map(|Probe { arg, val }| (*arg, *val));
        let Some(next_mod_rec) = self.rat.rat.rec_linear(pts) else {
            // TODO: better return value
            warn!("Failed reconstruction modulo {Q}");
            return;
        };
        debug!("Finished reconstruction modulo {Q}: {next_mod_rec}");
        self.rat.merge_crt(next_mod_rec);
        self.res = (&self.rat).try_into();
    }

    fn num_pows(&self) -> [u32; N] {
        self.powers.map(|d| d[0])
    }

    fn den_pows(&self) -> [u32; N] {
        self.powers.map(|d| d[1])
    }

    /// The current reconstruction status with suggestions for the next step
    pub fn status(&self) -> Status<P, N> {
        self.status
    }

    /// Extract the reconstructed rational function
    ///
    /// Returns `None` if `status()` is not `Done`
    pub fn into_rat(self) -> Option<Rat<FlatPoly<Integer, N>>> {
        self.res.ok()
    }

    /// The number of extra points used to validate the reconstruction
    pub fn extra_pts(&self) -> usize {
        self.extra_pts
    }

    /// The current modulus
    pub fn modulus(&self) -> &Integer {
        &self.rat.modulus
    }
}

fn nterms_with_max_pows<const N: usize>(max_pows: [u32; N]) -> usize {
    max_pows.iter().map(|&n| n as usize).product()
}

fn gen_poly_with_max_pows<const N: usize>(
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

#[cfg(test)]
mod tests {
    use super::*;
    use ::rand::{Rng, SeedableRng};

    use crate::algebra::poly::flat::FlatMono;
    use crate::rec::primes::LARGE_PRIMES;
    use crate::traits::{One, Zero};
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
                    const P: u64 = LARGE_PRIMES[0];

                    let orig = rand_rat::<NVARS>(&mut rng);
                    eprintln!("trying to reconstruct {orig}");
                    let mut rec = Rec::new(1);
                    let (mut z, mut q_z) = std::iter::repeat_with(|| {
                        let z: [Z64<P>; NVARS] = [(); NVARS].map(|_| rng.gen());
                        orig.try_eval(&z).map(|q_z| (z, q_z))
                    })
                        .flatten()
                        .next()
                        .unwrap();
                    while let Status::Degree(n) = rec.add_pt(z, q_z) {
                        z[n] += Z64::one();
                        q_z = loop {
                            if let Some(val) = orig.try_eval(&z) {
                                break val;
                            }
                            z[n] += Z64::one() ;
                        };
                    }
                    // TODO: re-use points
                    if let Status::Rat(mut npts) = rec.status() {
                        'rec: {
                            seq!{ N in 0..20 {{
                                const P: u64 = LARGE_PRIMES[N];
                                let pts = Vec::from_iter(
                                    std::iter::repeat_with(|| {
                                        let arg: [Z64<P>; NVARS] = [(); NVARS].map(|_| rng.gen());
                                        orig.try_eval(&arg).map(|val| Probe{arg, val})
                                    })
                                        .flatten()
                                        .take(npts)
                                );
                                match rec.add_pts(&pts) {
                                    Status::Rat(n) => npts = n,
                                    Status::Done => break 'rec,
                                    _ => unreachable!("Unexpected reconstruction return value")
                                }
                            }}}
                            let _ = npts;
                            panic!("Need more than 20 characteristics!");
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
