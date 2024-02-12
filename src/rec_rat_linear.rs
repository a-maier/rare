use std::ops::ControlFlow;

use ffnt::Z64;
use log::debug;
use paste::paste;
use rug::Integer;
use seq_macro::seq;

use crate::{
    rat::{NoneError, Rat},
    rec_linear::{RecLinear, UnknownDegreeRec},
    rec_rat::{combine_crt_rat, FFRat},
    sparse_poly::SparsePoly,
    traits::{TryEval, Zero},
};

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Needed {
    MoreThan(usize),
    Exact(usize),
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum ReconstructionStatus {
    FirstRat,
    Rat,
    Done,
}

seq! {N in 2..=16 {
    paste! {
        #[derive(Clone, Debug, PartialEq, PartialOrd)]
        pub struct [<RatRecLinear N>] {
            rec: UnknownDegreeRec,
            rat: FFRat<N>,
            res: Result<Rat<SparsePoly<Integer, N>>, NoneError>,
            done: bool,
        }

        impl [<RatRecLinear N>] {
            pub fn new(extra_pts: usize) -> Self {
                Self::with_degree_ratio(extra_pts, 1.)
            }

            pub fn with_degree_ratio(
                extra_pts: usize,
                degree_ratio: f64
            ) -> Self {
                Self {
                    rec: UnknownDegreeRec{extra_pts, degree_ratio},
                    rat: Default::default(),
                    res: Err(NoneError {  }),
                    done: false,
                }
            }

            pub fn add_pts<'a, const P: u64>(
                &mut self,
                pts: &[([Z64<P>; N], Z64<P>)]
            ) -> ControlFlow<(), Needed> {
                use ControlFlow::*;
                use ReconstructionStatus::*;
                let npts = pts.len();
                match self.status() {
                    FirstRat => if let Some(rat) = self.rec.rec(pts.iter().copied()) {
                        self.rat = FFRat::from(rat);
                        debug!("Finished reconstruction modulo {P}: {}", self.rat.rat);
                        self.res = (&self.rat).try_into();
                        if let Ok(rat) = self.res.as_ref() {
                            if rat.is_zero() {
                                self.done = true;
                            }
                        }
                        self.ask_for_new_mod()
                    } else {
                        Continue(Needed::MoreThan(npts))
                    },
                    Rat => self.rec_rat_mod_from_pts(pts),
                    Done => Break(()),
                }
            }

            pub fn rec_rat_mod_from_pts<const P: u64>(
                &mut self,
                pts: &[([Z64<P>; N], Z64<P>)],
            ) -> ControlFlow<(), Needed> {
                // if we can already reproduce the points we are done
                if let Ok(res) = self.res.as_ref() {
                    if pts.len() >= self.rec.extra_pts {
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
                debug!("Finished reconstruction modulo {P}: {next_mod_rec}");
                let mod_rec = std::mem::take(&mut self.rat);
                self.rat = combine_crt_rat(mod_rec, next_mod_rec);
                self.res  = (&self.rat).try_into();
                self.ask_for_new_mod()
            }

            pub fn status(&self) -> ReconstructionStatus {
                use ReconstructionStatus::*;
                if self.rat.modulus.is_zero() {
                    FirstRat
                } else if self.done {
                    Done
                } else {
                    Rat
                }
            }

            fn ask_for_new_mod(&self) -> ControlFlow<(), Needed> {
                let ncoeff = self.rat.rat.num().len() + self.rat.rat.den().len() - 1;
                if ncoeff == 0 {
                    ControlFlow::Break(())
                } else {
                    ControlFlow::Continue(Needed::Exact(ncoeff))
                }
            }

            pub fn into_rat(self) -> Option<Rat<SparsePoly<Integer, N>>> {
                self.res.ok()
            }
        }
    }
}}

#[cfg(test)]
mod tests {
    use super::*;
    use ::rand::{Rng, SeedableRng};

    use crate::rec_rat::LARGE_PRIMES;
    use crate::sparse_poly::SparseMono;
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
            fn [<rec_rat_explicit_ NVARS>]() {
                log_init();

                let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
                for _ in 0..NTESTS {
                    const P: u64 = LARGE_PRIMES[0];

                    let orig = rand_rat::<NVARS>(&mut rng);
                    eprintln!("trying to reconstruct {orig}");
                    let mut rec = [<RatRecLinear NVARS>]::new(1);
                    let mut pt_iter = std::iter::repeat_with(|| {
                        let z: [Z64<P>; NVARS] = [(); NVARS].map(|_| rng.gen());
                        orig.try_eval(&z).map(|q_z| (z, q_z))
                    })
                        .flatten();
                    let mut pts = Vec::from_iter(pt_iter.by_ref().take(3));
                    loop {
                        use ControlFlow::*;
                        use Needed::*;
                        match rec.add_pts(&pts) {
                            Continue(MoreThan(_)) => {
                                pts.push(pt_iter.next().unwrap())
                            },
                            Continue(Exact(n)) => {
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
                                    assert_eq!(next, Continue(Exact(n)))
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
