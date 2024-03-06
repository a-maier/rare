/// Implements the algorithm of
///
/// A. Cuyt, W.-s. Lee,
/// Sparse interpolation of multivariate rational functions
/// Theor. Comp. Sci. 412 (2011) 1445–1456.
///
/// as summarised in the FireFly paper:
///
/// J. Klappert, F. Lange
/// Reconstructing Rational Functions with FireFly
/// Comput.Phys.Commun. 264 (2021) 107968
///
use std::{num::NonZeroUsize, ops::ControlFlow};

use ffnt::Z64;
use log::{debug, error, trace, warn};
use paste::paste;
use rand::Rng;

use crate::{
    algebra::{
        poly::{dense::DensePoly, flat::{FlatPoly, FlatMono}},
        rat::Rat
    },
    arr::Arr,
    rec::rat::finite::{linear::{LinearRec, RecLinear}, thiele::ThieleRec},
    traits::{Eval, One, Rec, Shift, TryEval, WithVars, Zero},
    util::{slice_start, ALL_VARS_Z},
};

/// Multivariate rational function reconstruction over a finite field
// Note that this is mostly a marker struct!
// The actual reconstruction is implemented in RatRecMod2, RatRecMod3, etc.
#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct RatRecMod {
    extra_pts: usize,
}

impl RatRecMod {
    pub fn new(extra_pts: usize) -> Self {
        Self { extra_pts }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum ReconstructionStatus {
    Rat,
    NumPoly,
    DenPoly,
    Done,
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Needed<const P: u64, const N: usize> {
    Pt([Z64<P>; N]),
    Pts(Vec<[Z64<P>; N]>),
}

impl<const P: u64, const N: usize> From<[Z64<P>; N]> for Needed<P, N> {
    fn from(source: [Z64<P>; N]) -> Self {
        Self::Pt(source)
    }
}

impl<const P: u64, const N: usize> From<Vec<[Z64<P>; N]>> for Needed<P, N> {
    fn from(source: Vec<[Z64<P>; N]>) -> Self {
        Self::Pts(source)
    }
}

macro_rules! impl_rat_rec_mod_recursive {
    ( $($n:literal, $n_minus_one:literal), * ) => {
        $(
            paste! {
                use crate::rec::poly::finite::newton::[<NewtonPolyRec $n_minus_one>];
                use crate::algebra::poly::dense::[<DensePoly $n_minus_one>];

                /// Rational function reconstruction of $n variables over finite field with characteristic `P'
                // Note that the _internal_ normalisation convention
                // is different from all other parts of the code: the
                // coefficient of the *lowest-degree* term in the
                // *denominator* is set to one.
                // The rationale is that we can only fix the normalisation
                // for coefficients that don't depend on any of the variables
                // x0, x1, etc. With the current implementation of
                // shifts we can only guarantee that for the
                // lowest-degree term in the denominator.
                // The returned rational function follows the usual
                // convention that the highest-degree term in the
                // numerator is normalised to one.
                #[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
                pub struct [<RatRecMod $n>]<const P: u64> {
                    extra_pts: usize,
                    rat_in_t_rec: ThieleRec<P>,
                    poly_rec: [<NewtonPolyRec $n_minus_one>]<P>,
                    t_pow: usize,
                    pts_for_coeff: Vec<([Z64<P>; $n_minus_one], Z64<P>)>,
                    base_x_pt: [Z64<P>; $n],
                    x_offset: [Z64<P>; $n_minus_one],
                    shift: [Z64<P>; $n],
                    known_pts: Vec<RatPt<P, $n_minus_one>>,
                    num: FlatPoly<Z64<P>, $n>,
                    den: FlatPoly<Z64<P>, $n>,
                    status: ReconstructionStatus,
                    subtr: Vec<[<DensePoly $n_minus_one>]<Z64<P>>>
                }

                impl<const P: u64> [<RatRecMod $n>]<P> {
                    pub fn new(extra_pts: usize) -> Self {
                        Self::with_shift(extra_pts, [Z64::zero(); $n])
                    }

                    pub fn with_shift(extra_pts: usize, shift: [Z64<P>; $n]) -> Self {
                        Self {
                            extra_pts,
                            rat_in_t_rec: ThieleRec::new(extra_pts),
                            base_x_pt: [Z64::zero(); $n],
                            x_offset: [Z64::zero(); $n_minus_one],
                            shift,
                            known_pts: vec![],
                            poly_rec: [<NewtonPolyRec $n_minus_one>]::new(extra_pts),
                            t_pow: Default::default(),
                            pts_for_coeff: vec![],
                            num: Default::default(),
                            den: Default::default(),
                            status: ReconstructionStatus::Rat,
                            subtr: Default::default(),
                        }
                    }

                    /// Add a single point `(z, q(z))` to the reconstruction
                    ///
                    /// This should only by used if `status()` returns
                    /// `Rat`, preferably with a value of `z` that was
                    /// returned from a previous call
                    pub fn add_pt(
                        &mut self,
                        z: [Z64<P>; $n],
                        q_z: Z64<P>
                    ) -> ControlFlow<(), Needed<P, $n>> {
                        use ReconstructionStatus::*;
                        match self.status() {
                            Rat => self.add_pt_for_rat_in_t(z, q_z),
                            NumPoly | DenPoly => self.ask_for_new_rat_pt(),
                            Done => ControlFlow::Break(())
                        }
                    }

                    /// Add a slice of points `[(z, q(z))]` to the reconstruction
                    ///
                    /// This should only by used if `status()` returns
                    /// `NumPoly` or `DenPoly`. The length of the
                    /// slice and the values of `z` should correspond
                    /// to the previous return from `add_pt` or
                    /// `add_pts`.
                    pub fn add_pts(
                        &mut self,
                        pts: &[([Z64<P>; $n], Z64<P>)],
                    ) -> ControlFlow<(), Needed<P, $n>> {
                        use ReconstructionStatus::*;
                        match self.status() {
                            Rat => {
                                let mut res = ControlFlow::Break(());
                                for (z, q_z) in pts {
                                    res = self.add_pt(*z, *q_z);
                                    if self.status() != Rat {
                                        break;
                                    }
                                }
                                res
                            },
                            NumPoly | DenPoly => {
                                self.add_rec_rat_pt(pts)
                            },
                            Done => ControlFlow::Break(())
                        }
                    }

                    // Add a point for the Thiele reconstruction of
                    // q(x, t) = (a0(x)*t^0 + a1(x)*t^1 + ...) / (1 + b1(x)*t^1 + ...)
                    // in t
                    //
                    // If this completes the reconstruction in `t`, we
                    // immediately proceed to the reconstruction of
                    // the coefficients `a` (or finish with 0 if all a vanish identically)
                    fn add_pt_for_rat_in_t(
                        &mut self,
                        z: [Z64<P>; $n],
                        q_z: Z64<P>
                    ) -> ControlFlow<(), Needed<P, $n>> {
                        use ControlFlow::*;
                        let x_pt = z_to_x(z, &self.shift);
                        if self.rat_in_t_rec.rec_started() {
                            assert_eq!(x_pt[..$n_minus_one], self.base_x_pt[..$n_minus_one]);
                        } else {
                            debug!("Starting multivariate rational reconstruction with base point {x_pt:?}");
                            self.base_x_pt = x_pt;
                        }
                        let t = *x_pt.last().unwrap();
                        match self.rat_in_t_rec.add_pt(t, q_z) {
                            Continue(_) => {
                                let mut next_x = x_pt;
                                next_x[$n - 1] += Z64::one();
                                Continue(x_to_z(next_x, &self.shift).into())
                            },
                            Break(_) => {
                                let rat_in_t: Rat<_> = std::mem::take(&mut self.rat_in_t_rec)
                                    .into_rat()
                                    .into();
                                let rat_in_t = normalise_lowest_den_coeff_to_one(rat_in_t);

                                debug!("From Thiele reconstruction: {}", rat_in_t.with_vars(&["t"]));
                                debug_assert!(!rat_in_t.den().coeff(0).is_zero());
                                if rat_in_t.is_zero() {
                                    self.status = ReconstructionStatus::Done;
                                    self.den = One::one();
                                    return Break(());
                                }

                                // first point to be used for starting the polynomial
                                // reconstruction (numerator and denominator)
                                let nnum = rat_in_t.num().len();
                                assert!(nnum > 0);
                                self.known_pts.push(RatPt::new([Z64::zero(); $n_minus_one], rat_in_t));

                                self.subtr = vec![Zero::zero(); nnum];
                                self.t_pow = nnum - 1;
                                self.status = ReconstructionStatus::NumPoly;
                                self.init_poly_rec()
                            },
                        }
                    }

                    // If the reconstruction returned an `z` where `q(z)` is not defined
                    // you can call this function to get a new value of `z`
                    pub fn request_next_arg(
                        &self,
                        bad_z: [Z64<P>; $n]
                    ) -> [Z64<P>; $n] {
                        let mut xt = z_to_x(bad_z, &self.shift);
                        xt[$n - 1] += Z64::one();
                        x_to_z(xt, &self.shift)
                    }

                    // The value of `x` for the first point used for reconstruction
                    fn x_base(&self) -> [Z64<P>; $n_minus_one] {
                        self.base_x_pt[..$n_minus_one].try_into().unwrap()
                    }

                    // Initialise the reconstruction of a new coefficient `a` or `b` in
                    // q(x, t) = (a0(x)*t^0 + a1(x)*t^1 + ...) / (1 + b1(x)*t^1 + ...)
                    fn init_poly_rec(&mut self) -> ControlFlow<(), Needed<P, $n>> {
                        use ReconstructionStatus::*;

                        debug!("Reconstruct numerator coefficient of t^{}", self.t_pow);

                        let known_pts = self.known_pts.iter();
                        self.pts_for_coeff = if self.status() == NumPoly {
                            known_pts.map(|pt| (pt.offset, pt.num_coeff[self.t_pow])).collect()
                        } else {
                            debug_assert_eq!(self.status(), DenPoly);
                            known_pts.map(|pt| (pt.offset, pt.den_coeff[self.t_pow])).collect()
                        };
                        self.pts_for_coeff.sort_unstable_by(|a, b| b.cmp(a));
                        trace!("Reusing points: {:#?}", self.pts_for_coeff.len());
                        self.poly_rec = [<NewtonPolyRec $n_minus_one>]::<P>::new(self.extra_pts);
                        self.x_offset = [Z64::zero(); $n_minus_one];

                        self.rec_coeff_from_pts()
                    }

                    // Reconstruction of a coefficient `a` or `b` in
                    // q(x, t) = (a0(x)*t^0 + a1(x)*t^1 + ...) / (1 + b1(x)*t^1 + ...)
                    // possibly continuing a reconstruction that was
                    // interrupted due to lack of points.
                    fn rec_coeff_from_pts(&mut self) -> ControlFlow<(), Needed<P, $n>> {
                        while let Some((offset, mut val)) = self.pts_for_coeff.pop() {
                            use std::cmp::Ordering;
                            match offset.cmp(&self.x_offset) {
                                Ordering::Less => continue, // not needed for reconstruction
                                Ordering::Equal => {
                                    trace!("Using known value at offset {offset:?}: {val}");
                                    let x = (Arr(self.x_base()) + Arr(offset)).0;
                                    val -= self.subtr[self.t_pow].eval(&x);
                                    match self.poly_rec.add_pt(&x, val) {
                                        ControlFlow::Continue(n) => {
                                            self.x_offset[n] += Z64::one();
                                            self.x_offset[(n + 1)..].fill(Z64::zero());
                                        },
                                        ControlFlow::Break(()) => {
                                            return self.add_reconstructed_poly();
                                        }
                                    }
                                },
                                Ordering::Greater => {
                                    // TODO: allow offset > 1 in the
                                    // variable that should be
                                    // incremented
                                    self.pts_for_coeff.push((offset, val));
                                    break;
                                }
                            }
                        }
                        self.ask_for_new_rat_pt()
                    }

                    // Ask for a new set of coordinates z for a fixed value of x = c
                    // and varying t so that we can reconstruct
                    // q(x = c, t) = (a0(c)*t^0 + a1(c)*t^1 + ...) / (1 + b1(c)*t^1 + ...)
                    // where we know the degrees of the numerator and the denominator
                    // from earlier Thiele reconstruction
                    fn ask_for_new_rat_pt(&self) -> ControlFlow<(), Needed<P, $n>> {
                        // we need to reconstruct the rational function in t again
                        // for a new value of x calculated from the current offset
                        trace!("Need new value at offset {:?}", self.x_offset);
                        let nnum = self.known_pts[0].num_coeff.len();
                        let nden = self.known_pts[0].den_coeff.len();
                        let num_needed = (nnum + nden - 1) as u64;
                        assert!(num_needed <= P);
                        let mut xt = self.base_x_pt;
                        for (x, s) in xt.iter_mut().zip(self.x_offset.iter()) {
                            *x += s;
                        }
                        let needed = (0..num_needed).map(|t_offset| {
                            let mut xt = xt;
                            xt[$n - 1] += unsafe{ Z64::new_unchecked(t_offset) };
                            x_to_z(xt, &self.shift)
                        }).collect();
                        ControlFlow::Continue(Needed::Pts(needed))
                    }

                    // Add
                    // q(x = c, t) = (a0(c)*t^0 + a1(c)*t^1 + ...) / (1 + b1(c)*t^1 + ...)
                    // for a new fixed value of x = c and continue
                    // reconstructing the next coefficient `a` or `b`.
                    fn add_rec_rat_pt(
                        &mut self,
                        pts: &[([Z64<P>; $n], Z64<P>)],
                    ) -> ControlFlow<(), Needed<P, $n>> {
                        use ReconstructionStatus::*;

                        if pts.is_empty() {
                            error!("Received empty slice of points for reconstruction");
                            return self.ask_for_new_rat_pt();
                        }
                        let xt0 = z_to_x(pts[0].0, &self.shift);
                        let x0: [_; $n - 1] = xt0[..($n - 1)].try_into().unwrap();
                        let x_offset = (Arr(x0) - Arr(self.x_base())).0;
                        if let Some(next_pt) = self.pts_for_coeff.last() {
                            if next_pt.0 <= x_offset {
                                error!("Received function values at wrong values of argument");
                                return self.ask_for_new_rat_pt();
                            }
                        }

                        let mut t_pts = Vec::with_capacity(pts.len());
                        for (z, q_z) in pts {
                            let xt = z_to_x(*z, &self.shift);
                            if &xt[..($n - 1)] != &x0 {
                                error!("Received function values at wrong values of argument");
                                return self.ask_for_new_rat_pt();
                            }
                            t_pts.push((xt[$n - 1], *q_z));
                        }

                        // TODO: subtract the coefficients we have already reconstructed
                        let nnum = self.known_pts[0].num_coeff.len();
                        let nden = self.known_pts[0].den_coeff.len();
                        let linear_rec = LinearRec::new(
                            nnum,
                            NonZeroUsize::new(nden).unwrap()
                        );
                        let Some(rat) = linear_rec.rec_linear(t_pts) else {
                            warn!("Reconstruction of rational function failed! Trying again.");
                            return self.ask_for_new_rat_pt()
                        };
                        let rat = normalise_lowest_den_coeff_to_one(rat);
                        trace!("Reconstructed {}", rat.with_vars(&["t"]));

                        let coeff_val = if self.status() == NumPoly {
                            *rat.num().coeff(self.t_pow)
                        } else {
                            debug_assert_eq!(self.status(), DenPoly);
                            *rat.den().coeff(self.t_pow)
                        };
                        self.pts_for_coeff.push((x_offset, coeff_val));
                        let pt = RatPt::new(x_offset, rat);
                        self.known_pts.push(pt);
                        self.rec_coeff_from_pts()
                    }

                    // We finished the reconstruction of a coefficient `a` or `b` in
                    // q(x, t) = (a0(x)*t^0 + a1(x)*t^1 + ...) / (1 + b1(x)*t^1 + ...)
                    // Here, we turn it back from a(x)*t^? ot b(x)*t^?
                    // to a polynomial in z and add it to the
                    // previously reconstructed terms in the numerator
                    // / denominator. We then proceed with the next
                    // coefficient or declare victory if there are
                    // none left to reconstruct.
                    fn add_reconstructed_poly(&mut self) -> ControlFlow<(), Needed<P, $n>> {
                        use ReconstructionStatus::*;
                        let rec_poly: [<DensePoly $n_minus_one>]<Z64<P>> =
                            std::mem::take(&mut self.poly_rec).into_poly().into();
                        // turn it back from a polynomial in [x0, x1, ..., t] to [z0, z1, ...]
                        // first expand it out
                        let rec_poly: FlatPoly<Z64<P>, $n_minus_one> = rec_poly.into();
                        // rescale, adding back the power in t
                        let terms = rec_poly.into_terms().into_iter()
                            .map(|c| {
                                let mut powers = [0; $n];
                                powers[0..$n_minus_one].copy_from_slice(&c.powers);
                                powers[$n_minus_one] = self.t_pow as u32 - c.powers.iter().sum::<u32>();
                                FlatMono::new(c.coeff, powers)
                            });
                        let rec_poly = FlatPoly::from_terms(terms.collect());
                        debug!("Reconstructed contribution {}", rec_poly.with_vars(Self::zvars()));
                        if !self.shift.is_zero() {
                            self.add_shift_subtr_from(rec_poly.clone());
                        }

                        match self.status() {
                            NumPoly => self.num += rec_poly,
                            DenPoly => self.den += rec_poly,
                            _ => unreachable!()
                        }

                        if self.t_pow > 0 {
                            self.t_pow -= 1;
                            self.init_poly_rec()
                        } else if self.status() == NumPoly {
                            let nden = self.known_pts[0].den_coeff.len();
                            self.subtr = vec![Zero::zero(); nden];
                            self.t_pow = nden - 1;
                            self.status = DenPoly;
                            self.init_poly_rec()
                        } else {
                            self.status = Done;
                            ControlFlow::Break(())
                        }
                    }

                    fn add_shift_subtr_from(
                        &mut self,
                        poly: FlatPoly<Z64<P>, $n>
                    ) {
                        use crate::algebra::poly::dense::[<DensePoly $n>];

                        // calculate P(z + s) - P(Z)
                        let collected = [<DensePoly $n>]::from(poly);
                        let shifted = collected.clone().shift(self.shift);
                        let subtr = shifted - &collected;
                        trace!("shift subtraction in z: {}", subtr.with_vars(Self::zvars()));
                        // turn it into a polynomial in [t, x0, x1, ...]
                        let subtr: FlatPoly<Z64<P>, $n> = subtr.into();
                        let terms = subtr.into_terms().into_iter()
                            .map(|c| {
                                let mut powers = [0; $n];
                                powers[0] = c.powers.iter().sum::<u32>();
                                powers[1..].copy_from_slice(&c.powers[0..$n_minus_one]);
                                FlatMono::new(c.coeff, powers)
                            });
                        let subtr = FlatPoly::from_terms(terms.collect());
                        let subtr = [<DensePoly $n>]::from(subtr);
                        trace!("shift subtraction in [t, x1, ...]: {subtr}");
                        let new_subtr = subtr.into_coeff().into_iter();
                        for (subtr, new_subtr) in self.subtr.iter_mut().zip(new_subtr) {
                            *subtr += new_subtr;
                        }
                    }

                    pub fn status(&self) -> ReconstructionStatus {
                        self.status
                    }

                    pub fn into_rat(self) -> Rat<FlatPoly<Z64<P>, $n>> {
                        if self.status() != ReconstructionStatus::Done {
                            warn!("Taking out rational function before reconstruction is finished");
                        }
                        let Self{mut num, mut den, ..} = self;
                        let Some(last) = num.terms().last() else {
                            return Zero::zero();
                        };
                        let norm = last.coeff.inv();
                        num *= norm;
                        den *= norm;
                        Rat::from_num_den_unchecked(num, den)
                    }

                    fn zvars() -> &'static [&'static str; $n] {
                        slice_start(&ALL_VARS_Z)
                    }
                }
            }

        )*
    };
}

impl_rat_rec_mod_recursive! {16, 15, 15, 14, 14, 13, 13, 12, 12, 11, 11, 10, 10, 9, 9, 8, 8, 7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1}

fn normalise_lowest_den_coeff_to_one<const P: u64>(
    rat: Rat<DensePoly<Z64<P>>>,
) -> Rat<DensePoly<Z64<P>>> {
    let (mut num, mut den) = rat.into_num_den();
    let norm = den.coeff(0).inv();
    num *= &norm;
    den *= &norm;
    Rat::from_num_den_unchecked(num, den)
}

// transform from [x0, x1, ..., t] to
// [z0, z1, z2, ...] = [t*x0 + s0, t*x1 + s1, t*x2 + s2, ..., t + sM]
// with a `shift` = [s0, s1, s2, ..., sM]
fn x_to_z<const N: usize, const P: u64>(
    mut coord: [Z64<P>; N],
    shift: &[Z64<P>; N],
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

// transform from [z0, z1, ...] to [x0, x1, ..., t]
// with a `shift` = [s0, s1, s2, ..., sM]
fn z_to_x<const N: usize, const P: u64>(
    mut coord: [Z64<P>; N],
    shift: &[Z64<P>; N],
) -> [Z64<P>; N] {
    for (z, s) in coord.iter_mut().zip(shift.iter()) {
        *z -= s;
    }
    let (t, zs) = coord.split_last_mut().unwrap();
    let tinv = t.inv();
    for z in zs {
        *z *= tinv
    }
    coord
}

// find a shift such that the shifted function fs(z) = f(z + s)
// has a constant leading term in the denominator
// In other words, find s such that f(s) is well-defined
pub fn find_shift<F, const P: u64, const N: usize>(
    mut f: F,
    mut rng: impl Rng,
) -> [Z64<P>; N]
where
    F: FnMut([Z64<P>; N]) -> Option<Z64<P>>,
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

macro_rules! impl_rec_with_ran {
    ( $($x:literal), *) => {
        $(
            paste! {
                impl<F, const P: u64> Rec<RatRecMod, [Z64<P>; $x]> for F
                where F: FnMut([Z64<P>; $x]) -> Option<Z64<P>> {
                    type Output = Option<Rat<FlatPoly<Z64<P>, $x>>>;

                    fn rec_with_ran(
                        &mut self,
                        rec: RatRecMod,
                        mut rng: impl Rng
                    ) -> Self::Output {
                        // it's impossible to pass &mut self to `find_shift` in rust 1.75.0
                        // this is a workaround
                        let mut tmp = |z| (self)(z);
                        let shift = find_shift(&mut tmp, &mut rng);
                        let mut rec = [<RatRecMod $x>]::with_shift(rec.extra_pts, shift);
                        let (z, q_z) = loop {
                            let z: [Z64<P>; $x] = [(); $x].map(|_| rng.gen());
                            if let Some(q_z) = (self)(z) {
                                break (z, q_z);
                            }
                        };
                        let mut next = rec.add_pt(z, q_z);
                        loop {
                            use ControlFlow::*;
                            match next {
                                Continue(Needed::Pt(mut z)) => {
                                    let q_z = loop {
                                        if let Some(q_z) = (self)(z) {
                                            break q_z;
                                        }
                                        z = rec.request_next_arg(z);
                                    };
                                    next = rec.add_pt(z, q_z);
                                },
                                Continue(Needed::Pts(zs)) => {
                                    let mut pts = Vec::from_iter(
                                        zs.iter().copied()
                                            .filter_map(|z| (self)(z).map(|q_z| (z, q_z)))
                                    );
                                    let mut z = *zs.last().unwrap();
                                    while pts.len() < zs.len() {
                                        z = rec.request_next_arg(z);
                                        if let Some(q_z) = (self)(z) {
                                            pts.push((z, q_z))
                                        }
                                    }
                                    next = rec.add_pts(&pts);
                                },
                                Break(()) => {
                                    return Some(rec.into_rat());
                                }
                            }
                        }
                    }
                }

                impl<F, const P: u64> Rec<RatRecMod, [[Z64<P>; $x]; 1]> for F
                where F: TryEval<[Z64<P>; $x], Output = Z64<P>> {
                    type Output = Option<Rat<FlatPoly<Z64<P>, $x>>>;

                    fn rec_with_ran(
                        &mut self,
                        rec: RatRecMod,
                        rng: impl Rng
                    ) -> Self::Output {
                        (|x: [Z64<P>; $x]| self.try_eval(&x)).rec_with_ran(rec, rng)
                    }
                }
            }
        )*
    };
}

impl<F, const P: u64> Rec<RatRecMod, [Z64<P>; 1]> for F
where
    F: FnMut([Z64<P>; 1]) -> Option<Z64<P>>,
{
    type Output = Option<Rat<FlatPoly<Z64<P>, 1>>>;

    fn rec_with_ran(&mut self, rec: RatRecMod, rng: impl Rng) -> Self::Output {
        let rec = ThieleRec::new(rec.extra_pts);
        rec.rec_univariate_with_ran(|x: Z64<P>| (self)([x]), rng)
            .map(|r| r.into())
    }
}

impl<F, const P: u64> Rec<RatRecMod, [[Z64<P>; 1]; 1]> for F
where
    F: TryEval<[Z64<P>; 1], Output = Z64<P>>,
{
    type Output = Option<Rat<FlatPoly<Z64<P>, 1>>>;

    fn rec_with_ran(&mut self, rec: RatRecMod, rng: impl Rng) -> Self::Output {
        (|x: [Z64<P>; 1]| self.try_eval(&x)).rec_with_ran(rec, rng)
    }
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

impl_rec_with_ran!(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);

#[cfg(test)]
mod tests {
    use std::iter::repeat_with;

    use rand::Rng;
    use rand_xoshiro::rand_core::SeedableRng;

    use crate::{
        sparse_poly::FlatMono,
        traits::{One, TryEval, Zero},
    };

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
        let rec = RatRecMod::new(1);

        fn ran_coeff(mut rng: impl Rng) -> Z64<P> {
            let z = rng.gen_range(0..MAX_COEFF);
            unsafe { Z64::new_unchecked(z) }
        }

        for _ in 0..NTESTS {
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let coeff = repeat_with(|| {
                FlatMono::new(
                    ran_coeff(&mut rng),
                    [(); 2].map(|_| rng.gen_range(0..=MAX_POW)),
                )
            })
            .take(nterms)
            .collect();
            let num = FlatPoly::from_terms(coeff);

            let den = if num.is_zero() {
                One::one()
            } else {
                let mut den: FlatPoly<_, 2> = Zero::zero();
                while den.is_zero() {
                    let coeff = repeat_with(|| {
                        FlatMono::new(
                            ran_coeff(&mut rng),
                            [(); 2].map(|_| rng.gen_range(0..=MAX_POW)),
                        )
                    })
                    .take(nterms)
                    .collect();
                    den = FlatPoly::from_terms(coeff);
                }
                // ensure we don't need a shift
                let mut den = den.into_terms();
                den[0].powers = Zero::zero();
                den[0].coeff = One::one();
                FlatPoly::from_terms(den)
            };

            let rat = Rat::from_num_den_unchecked(num, den);
            eprintln!("trying to reconstruct {rat}");

            let res =
                (|z| rat.try_eval(&z)).rec_with_ran(rec, &mut rng).unwrap();
            eprintln!("reconstructed {res}");
            assert!(sample_eq(&rat, &res, &mut rng));
        }
    }

    #[test]
    fn rec_rat2_simple() {
        log_init();

        const P: u64 = 1152921504606846883;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = RatRecMod::new(1);

        let x = FlatMono::new(Z64::one(), [1, 0]);
        let y = FlatMono::new(Z64::one(), [0, 1]);

        // 1 / y
        let num = FlatPoly::<Z64<P>, 2>::one();
        let den = FlatPoly::<Z64<P>, 2>::from_terms(vec![y]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let res = (|z: [Z64<P>; 2]| rat.try_eval(&z))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        eprintln!("reconstructed {res}");
        assert!(sample_eq(&rat, &res, &mut rng));

        // 1 / y^2
        let num = FlatPoly::<Z64<P>, 2>::one();
        let den = FlatPoly::<Z64<P>, 2>::from_terms(vec![y * y]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let res = (|z: [Z64<P>; 2]| rat.try_eval(&z))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        eprintln!("reconstructed {res}");
        assert!(sample_eq(&rat, &res, &mut rng));

        // 1 / x
        let num = FlatPoly::<Z64<P>, 2>::one();
        let den = FlatPoly::<Z64<P>, 2>::from_terms(vec![x]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let res = (|z: [Z64<P>; 2]| rat.try_eval(&z))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        eprintln!("reconstructed {res}");
        assert!(sample_eq(&rat, &res, &mut rng));
    }

    #[test]
    fn rec_rat2_small_shift() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_POW: u32 = 2;
        const MAX_COEFF: u64 = 29;
        const P: u64 = 1152921504606846883;
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = RatRecMod::new(1);

        fn ran_coeff(mut rng: impl Rng) -> Z64<P> {
            let z = rng.gen_range(0..MAX_COEFF);
            unsafe { Z64::new_unchecked(z) }
        }

        for _ in 0..NTESTS {
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let coeff = repeat_with(|| {
                FlatMono::new(
                    ran_coeff(&mut rng),
                    [(); 2].map(|_| rng.gen_range(0..=MAX_POW)),
                )
            })
            .take(nterms)
            .collect();
            let num = FlatPoly::from_terms(coeff);

            let den = if num.is_zero() {
                One::one()
            } else {
                let mut den: FlatPoly<_, 2> = Zero::zero();
                while den.is_zero() {
                    let coeff = repeat_with(|| {
                        FlatMono::new(
                            ran_coeff(&mut rng),
                            [(); 2].map(|_| rng.gen_range(0..=MAX_POW)),
                        )
                    })
                    .take(nterms)
                    .collect();
                    den = FlatPoly::from_terms(coeff);
                }
                den
            };
            let rat = Rat::from_num_den_unchecked(num, den);
            eprintln!("trying to reconstruct {rat}");

            let res =
                (|z| rat.try_eval(&z)).rec_with_ran(rec, &mut rng).unwrap();
            eprintln!("reconstructed {res}");
            assert!(sample_eq(&rat, &res, &mut rng));
        }
    }

    #[test]
    fn rec_rat3() {
        log_init();

        const NTESTS: u32 = 40;
        const MAX_POW: u32 = 2;
        const MAX_COEFF: u64 = 29;
        const P: u64 = 1152921504606846883;
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let rec = RatRecMod::new(1);

        fn ran_coeff(mut rng: impl Rng) -> Z64<P> {
            let z = rng.gen_range(0..MAX_COEFF);
            unsafe { Z64::new_unchecked(z) }
        }

        for _ in 0..NTESTS {
            let max_pow = rng.gen_range(0..=MAX_POW);
            let nterms = 2usize.pow(max_pow);
            let coeff = repeat_with(|| {
                FlatMono::new(
                    ran_coeff(&mut rng),
                    [(); 3].map(|_| rng.gen_range(0..=MAX_POW)),
                )
            })
            .take(nterms)
            .collect();
            let num = FlatPoly::from_terms(coeff);

            let den = if num.is_zero() {
                One::one()
            } else {
                let mut den: FlatPoly<_, 3> = Zero::zero();
                while den.is_zero() {
                    let coeff = repeat_with(|| {
                        FlatMono::new(
                            ran_coeff(&mut rng),
                            [(); 3].map(|_| rng.gen_range(0..=MAX_POW)),
                        )
                    })
                    .take(nterms)
                    .collect();
                    den = FlatPoly::from_terms(coeff);
                }
                den
            };

            let rat = Rat::from_num_den_unchecked(num, den);
            eprintln!("trying to reconstruct {rat}");

            let res =
                (|z| rat.try_eval(&z)).rec_with_ran(rec, &mut rng).unwrap();
            eprintln!("reconstructed {res}");
            assert!(sample_eq3(&rat, &res, &mut rng));
        }
    }

    fn sample_eq<const P: u64>(
        orig: &Rat<FlatPoly<Z64<P>, 2>>,
        rec: &Rat<FlatPoly<Z64<P>, 2>>,
        mut rng: impl Rng,
    ) -> bool {
        const TESTS: usize = 10;
        for _ in 0..TESTS {
            let pt: [Z64<P>; 2] = rng.gen();
            if orig.try_eval(&pt) != rec.try_eval(&pt) {
                return false;
            }
        }
        true
    }

    fn sample_eq3<const P: u64>(
        orig: &Rat<FlatPoly<Z64<P>, 3>>,
        rec: &Rat<FlatPoly<Z64<P>, 3>>,
        mut rng: impl Rng,
    ) -> bool {
        const TESTS: usize = 10;
        for _ in 0..TESTS {
            let pt: [Z64<P>; 3] = rng.gen();
            if orig.try_eval(&pt) != rec.try_eval(&pt) {
                return false;
            }
        }
        true
    }
}
