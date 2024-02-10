use std::num::NonZeroUsize;

use ffnt::Z64;
use log::{debug, trace};

use crate::{
    matrix::Matrix,
    rand::pt_iter,
    rat::Rat,
    sparse_poly::{SparseMono, SparsePoly},
    traits::{Eval, One, Zero, Rec}, dense_poly::DensePoly,
};

/// Reconstruction using a linear system of equations built from the given points
pub(crate) trait RecLinear<Pts> {
    type Output;

    /// Attempt the reconstruction
    fn rec_linear(&self, pts: Pts) -> Option<Self::Output>;
}

/// Reconstruct a rational function by solving linear systems of equations
///
/// The degrees of both numerator and denominator have to be known.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct LinearRec {
    num_len: usize,
    den_len: NonZeroUsize,
}

impl LinearRec {
    pub fn new(num_len: usize, den_len: NonZeroUsize) -> Self {
        Self { num_len, den_len }
    }
}

impl<Pts, const P: u64> RecLinear<Pts> for LinearRec
where
    Pts: IntoIterator<Item = (Z64<P>, Z64<P>)>,
{
    type Output = Rat<DensePoly<Z64<P>>>;

    fn rec_linear(
        &self,
        pts: Pts
    ) -> Option<Rat<DensePoly<Z64<P>>>> {
        let pts = pts.into_iter();
        debug!(
            "1d rational function reconstruction with known degrees {}/{}",
            self.num_len as i64 - 1,
            usize::from(self.den_len) - 1
        );
        if self.num_len == 0 {
            return Some(Zero::zero());
        }
        let num_coeffs = self.num_len + usize::from(self.den_len);
        let neqs = std::cmp::max(pts.size_hint().0, num_coeffs - 1);
        debug!("Reconstructing with {neqs} equations");
        let mut eqs = Vec::with_capacity(neqs * num_coeffs);
        for (x, q_x) in pts.take(neqs) {
            let mut den_coeff = -q_x;
            for _ in 0..usize::from(self.den_len) {
                eqs.push(den_coeff);
                den_coeff *= x;
            }
            let mut x_to_i = One::one();
            for _ in 0..self.num_len {
                eqs.push(x_to_i);
                x_to_i *= x;
            }
        }
        let (den_coeffs, num_coeffs) = solve_eqs(
            neqs,
            eqs,
            usize::from(self.den_len),
        )?;
        let num = DensePoly::from_coeff(num_coeffs);
        let den = DensePoly::from_coeff(den_coeffs);
        Some(Rat::from_num_den_unchecked(num, den))
    }
}

impl<Pts, T, const P: u64, const N: usize> RecLinear<Pts> for Rat<SparsePoly<T, N>>
where
    Pts: IntoIterator<Item = ([Z64<P>; N], Z64<P>)>,
{
    type Output = Rat<SparsePoly<Z64<P>, N>>;

    fn rec_linear(
        &self,
        pts: Pts
    ) -> Option<Rat<SparsePoly<Z64<P>, N>>> {
        let pts = pts.into_iter();
        assert!(!self.den().is_empty());
        if self.num().len() == 0 {
            return Some(Zero::zero());
        }
        let num_coeffs = self.num().len() + self.den().len();
        let neqs = std::cmp::max(pts.size_hint().0, num_coeffs - 1);
        debug!("Reconstructing with {neqs} equations");

        // construct equations num(x) - q(x) * den(x) from known points x, q(x)
        // by convention the coefficient of the highest-order monomial in
        // the numerator is set to one
        let mut eqs = Vec::with_capacity(neqs * num_coeffs);
        for (x, q_x) in pts.take(neqs) {
            for term in self.den().terms() {
                eqs.push(-q_x * eval_pow(term, &x));
            }
            for term in self.num().terms() {
                eqs.push(eval_pow(term, &x));
            }
        }
        let (den_coeffs, num_coeffs) = solve_eqs(
            neqs,
            eqs,
            self.den().len()
        )?;
        let num = num_coeffs
            .into_iter()
            .zip(self.num().terms().iter().map(|t| t.powers))
            .filter_map(|(coeff, powers)| if coeff.is_zero() {
                None
            } else {
                Some(SparseMono { powers, coeff })
            })
            .collect();
        let num = SparsePoly::from_raw_terms(num);
        let den = den_coeffs
            .into_iter()
            .zip(self.den().terms().iter().map(|t| t.powers))
            .filter_map(|(coeff, powers)| if coeff.is_zero() {
                None
            } else {
                Some(SparseMono { powers, coeff })
            })
            .collect();
        let den = SparsePoly::from_raw_terms(den);
        Some(Rat::from_num_den_unchecked(num, den))
    }
}

// Solve `neqs` linear equations Σ c_i x_i = 0
//
// Returns `None` iff there is no non-trivial solution. Otherwise
// returns one non-trivial solution with the first `nden` x_i in one
// vector and the remainder in a second vector. The solution is chosen
// such that as many x_i with the highest i as possible are set to
// zero and one set to unity.
fn solve_eqs<const P: u64>(
    neqs: usize,
    eqs: Vec<Z64<P>>,
    nden: usize,
) -> Option<(Vec<Z64<P>>, Vec<Z64<P>>)> {
    debug_assert!(nden <= neqs);
    // if eqs.len() < neqs * (neqs + 1) {
    //     debug!("Found {} points, need at least {}", eqs.len(), neqs * (neqs + 1));
    //     return None;
    // }
    let mut eqs = Matrix::from_vec(neqs, eqs);
    eqs.row_reduce();
    eqs.trim_end();
    // Set as many variables as possible to zero, starting with the
    // last one. We can set nth variable to zero if
    // 1. The variable is not already expressed in terms of other
    //    variables.
    // 2. Setting the variable to zero mustn't turn the whole system
    //    zero.
    // Setting a variable to zero means setting the nth column to zero.
    debug!("After row reduction:\n{eqs}");
    let mut coeffs = vec![Z64::zero(); eqs.ncols()];
    for n in (0..eqs.ncols()).rev() {
        if let Some(eq) = find_solved_for(&eqs, n) {
            // all variables on the rhs have already been set to 0 or 1
            coeffs[n] = - eq[n+1..].iter().fold(Z64::zero(), |a, b| a + b);
        } else if can_set_to_zero(&eqs, n) {
            set_var_to_zero(&mut eqs, n);
        } else {
            coeffs[n] = One::one();
        }
    }
    let num_coeffs = coeffs.split_off(nden);
    let den_coeffs = coeffs;
    if den_coeffs.iter().all(|c| c.is_zero()) {
        None
    } else {
        Some((den_coeffs, num_coeffs))
    }
}

fn set_var_to_zero<const P: u64>(eqs: &mut Matrix<Z64<P>>, n: usize) {
    for r in eqs.rows_mut() {
        r[n] = Zero::zero();
    }
}

// check which equation is the solution for the nth variable
// i.e. the nth entry is the first non-vanishing one
fn find_solved_for<const P: u64>(
    eqs: &Matrix<Z64<P>>,
    n: usize
) -> Option<&[Z64<P>]> {
    eqs.rows().find(|eq| eq.iter().position(|e| !e.is_zero()) == Some(n))
}

// Check if we can set a variable to zero without making all zeroes
// the only solution. That means there has to be at least one row with
// more than one non-vanishing entry left.
fn can_set_to_zero<const P: u64>(eqs: &Matrix<Z64<P>>, n: usize) -> bool {
    eqs.rows().any(
        |eq| eq.iter().enumerate().filter(|(c, e)| !e.is_zero() && *c != n).count() > 1
    )
}

fn eval_pow<T, const P: u64, const N: usize>(
    term: &SparseMono<T, N>,
    x: &[Z64<P>; N],
) -> Z64<P> {
    SparseMono {
        coeff: Z64::<P>::one(),
        powers: term.powers,
    }
    .eval(x)
}

impl<F, const P: u64> Rec<LinearRec, Z64<P>> for F
where
    F: FnMut(Z64<P>) -> Option<Z64<P>>,
{
    type Output = Option<Rat<DensePoly<Z64<P>>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: LinearRec,
        rng: impl ::rand::Rng,
    ) -> Self::Output {
        reconstructor.rec_linear(
            pt_iter(rng).filter_map(|pt| (self)(pt).map(|fy| (pt, fy))),
        )
    }
}

impl<F, const P: u64> Rec<LinearRec, [Z64<P>; 1]> for F
where
    F: FnMut([Z64<P>; 1]) -> Option<Z64<P>>,
{
    type Output = Option<Rat<DensePoly<Z64<P>>>>;

    fn rec_with_ran(
        &mut self,
        reconstructor: LinearRec,
        rng: impl ::rand::Rng,
    ) -> Self::Output {
        (|pt| (self)([pt])).rec_with_ran(reconstructor, rng)
    }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct Unit { }

impl Zero for Unit {
    fn zero() -> Self {
        unimplemented!("Internal logic error: `Unit::zero()` must never be called")
    }

    fn is_zero(&self) -> bool {
        false
    }
}

const UNIT: Unit = Unit{};

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct UnknownDegreeRec {
    pub degree_ratio: f64,
    pub extra_pts: usize,
}

impl UnknownDegreeRec {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn rec<Pts, const P: u64, const N: usize>(
        &self,
        pts: Pts
    ) -> Option<Rat<SparsePoly<Z64<P>, N>>>
    where
        Pts: Iterator<Item = ([Z64<P>; N], Z64<P>)> + ExactSizeIterator
    {
        if pts.len() < self.extra_pts + 1 {
            return None;
        }
        let ncoeff_total = 1 + pts.len() - self.extra_pts;
        let num_frac = self.degree_ratio / (1. + self.degree_ratio);
        let ncoeff_num = (num_frac * ncoeff_total as f64) as usize;
        let ncoeff_den = ncoeff_total - ncoeff_num;
        if ncoeff_num < 1 || ncoeff_den < 1 {
            return None;
        }
        let num_ansatz = PowerIter::new()
            .take(ncoeff_num)
            .map(|p| SparseMono::new(UNIT, p))
            .collect();
        let num_ansatz = SparsePoly::from_raw_terms(num_ansatz);
        let den_ansatz = PowerIter::new()
            .take(ncoeff_den)
            .map(|p| SparseMono::new(UNIT, p))
            .collect();
        let den_ansatz = SparsePoly::from_raw_terms(den_ansatz);
        let rat = Rat::from_num_den_unchecked(num_ansatz, den_ansatz);
        debug_assert_eq!(rat.num().len() + rat.den().len(), ncoeff_total);
        trace!("Reconstruction ansatz: {rat:?}");
        rat.rec_linear(pts)
    }
}

impl Default for UnknownDegreeRec {
    fn default() -> Self {
        Self {
            degree_ratio: 1.,
            extra_pts: 1
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct PowerIter<const N: usize> (
    [u32; N]
);

impl<const N: usize> PowerIter<N> {
    fn new() -> Self {
        Self::default()
    }
}

impl<const N: usize> Default for PowerIter<N> {
    fn default() -> Self {
        Self([0; N])
    }
}

impl<const N: usize> Iterator for PowerIter<N> {
    type Item = [u32; N];

    fn next(&mut self) -> Option<Self::Item> {
        let res = self.0;
        let pos = self.0.iter().rposition(|&e| e > 0).unwrap_or(0);
        if pos > 0 {
            let val = std::mem::take(&mut self.0[pos]);
            *self.0.last_mut().unwrap() = val - 1;
            self.0[pos - 1] += 1;
        } else {
            let val = std::mem::take(&mut self.0[0]);
            *self.0.last_mut().unwrap() = val + 1;
        }
        Some(res)
    }
}

#[cfg(test)]
mod tests {
    use std::iter::repeat_with;

    use rand::{SeedableRng, Rng};

    use crate::{_test_util::{gen_sparse_rat, sample_eq, gen_dense_rat1}, traits::TryEval};

    use super::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    fn gen_pts<const P: u64, const N: usize>(
        rat: &impl TryEval<[Z64<P>; N], Output = Z64<P>>,
        num: usize,
        mut rng: impl Rng,
    ) -> Vec<([Z64<P>; N], Z64<P>)> {
        let mut pts = Vec::with_capacity(num);
        while pts.len() < num {
            let pt = [(); N].map(|_| rng.gen::<Z64<P>>());
            if let Some(val) = rat.try_eval(&pt) {
                pts.push((pt, val))
            }
        }
        pts
    }

    #[test]
    fn rec_rat2_small() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_POW: u32 = 2;
        const P: u64 = 1152921504606846883;
        const N: usize = 2;
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        for _ in 0..NTESTS {
            let rat = gen_sparse_rat::<P, N>(MAX_POW, MAX_POW, &mut rng);
            eprintln!("trying to reconstruct {rat}");
            let nneeded = rat.num().len() + rat.den().len() - 1;
            let pts = gen_pts(&rat, nneeded, &mut rng);
            let rec = rat.rec_linear(pts).unwrap();
            eprintln!("reconstructed {rec}");
            assert!(sample_eq(&rat, &rec, &mut rng));
        }
    }

    #[test]
    fn rec_rat2_small_unknown_degree() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_POW: u32 = 2; // update also MAX_PTS_NEEDED if this is changed!
        const P: u64 = 1152921504606846883;
        const N: usize = 2;
        const EXTRA_PTS: usize = 2;
        const MAX_PTS_NEEDED: usize = 25 + EXTRA_PTS;
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);
        let mut rec = UnknownDegreeRec::new();
        rec.extra_pts = EXTRA_PTS;
        for _ in 0..NTESTS {
            let rat = gen_sparse_rat::<P, N>(MAX_POW, MAX_POW, &mut rng);
            eprintln!("trying to reconstruct {rat}");
            let min_nneeded = rat.num().len() + rat.den().len() - 1;
            let pts = gen_pts(&rat, min_nneeded, &mut rng);
            let res = rec.rec(pts.into_iter());
            assert!(res.is_none());
            let pts = gen_pts(&rat, MAX_PTS_NEEDED, &mut rng);
            eprintln!("trying with {} points", pts.len());
            let res = rec.rec(pts.into_iter()).unwrap();
            eprintln!("reconstructed {res}");
            assert!(sample_eq(&rat, &res, &mut rng));
        }
    }

    // TODO: duplication with rec_rat_mod
    #[test]
    fn rec_rat2_simple() {
        log_init();

        const P: u64 = 1152921504606846883;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

        let x = SparseMono::new(Z64::one(), [1, 0]);
        let y = SparseMono::new(Z64::one(), [0, 1]);

        // 1 / y
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![y]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let nneeded = rat.num().len() + rat.den().len() - 1;
        let pts = gen_pts(&rat, nneeded, &mut rng);
        let res = rat.rec_linear(pts).unwrap();
        eprintln!("reconstructed {res}");
        assert!(sample_eq(&rat, &res, &mut rng));

        // 1 / y^2
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![y * y]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let nneeded = rat.num().len() + rat.den().len() - 1;
        let pts = gen_pts(&rat, nneeded, &mut rng);
        let res = rat.rec_linear(pts).unwrap();
        eprintln!("reconstructed {res}");
        assert!(sample_eq(&rat, &res, &mut rng));

        // 1 / x
        let num = SparsePoly::<Z64<P>, 2>::one();
        let den = SparsePoly::<Z64<P>, 2>::from_terms(vec![x]);
        let rat = Rat::from_num_den_unchecked(num, den);
        eprintln!("trying to reconstruct {rat}");
        let nneeded = rat.num().len() + rat.den().len() - 1;
        let pts = gen_pts(&rat, nneeded, &mut rng);
        let res = rat.rec_linear(pts).unwrap();
        eprintln!("reconstructed {res}");
        assert!(sample_eq(&rat, &res, &mut rng));
    }

    #[test]
    fn rec_rat_small() {
        log_init();

        const NTESTS: u32 = 10;
        const MAX_POW: u32 = 1;
        const P: u64 = 29;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

        for _ in 0..NTESTS {
            let rat = gen_dense_rat1(&[MAX_POW], &mut rng);
            eprintln!("trying to reconstruct {rat}");
            let rec = LinearRec::new(
                rat.num().len(),
                rat.den().len().try_into().unwrap(),
            );
            let reconstructed = (|x: Z64<P>| rat.try_eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            eprintln!("{reconstructed}");
            let reconstructed: Rat<DensePoly<Z64<P>>> = reconstructed.into();
            eprintln!("{reconstructed}");
            assert!(sample_eq(&rat, &reconstructed, &mut rng))
        }
    }

    #[test]
    fn rec_rat_large() {
        log_init();

        const NTESTS: u32 = 100;
        const MAX_POW: u32 = 8;
        const P: u64 = 1152921504606846883;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

        for _ in 0..NTESTS {
            let rat = gen_dense_rat1(&[MAX_POW], &mut rng);
            let rec = LinearRec::new(
                rat.num().len(),
                rat.den().len().try_into().unwrap(),
            );
            let reconstructed = (|x: Z64<P>| rat.try_eval(&x))
                .rec_with_ran(rec, &mut rng)
                .unwrap();
            let reconstructed: Rat<DensePoly<Z64<P>>> = reconstructed.into();
            assert!(sample_eq(&rat, &reconstructed, &mut rng))
        }
    }

    #[test]
    fn rec_unknown_degree() {
        log_init();
        const P: u64 = 1152921504606846883;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

        let rat = Rat::from_num_den_unchecked(
            DensePoly::from_coeff(vec![Z64::one(), Z64::one()]),
            DensePoly::one()
        );

        let rec = LinearRec::new(3, 3.try_into().unwrap());
        let reconstructed = (|x: Z64<P>| rat.try_eval(&x))
            .rec_with_ran(rec, &mut rng)
            .unwrap();
        assert_eq!(rat, reconstructed);
    }

    #[test]
    fn failed_rec() {
        log_init();
        const P: u64 = 1152921504606846883;

        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(1);

        let rat = Rat::from_num_den_unchecked(
            DensePoly::from_coeff(vec![Z64::<P>::one(), Z64::one()]),
            DensePoly::one()
        );

        let rec = LinearRec::new(1, 1.try_into().unwrap());
        let pts = Vec::from_iter(
            repeat_with(|| {
                let x: Z64<P> = rng.gen();
                rat.try_eval(&x).map(|q_x| (x, q_x))
            })
            .filter_map(|x| x)
            .take(2)
        );
        let rec = rec.rec_linear(pts);
        assert!(rec.is_none());
    }

}
