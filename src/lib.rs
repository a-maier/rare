/// Utilities for tests and benchmarks
#[doc(hidden)]
pub mod _test_util;
mod arr;
#[cfg(feature = "linear")]
mod matrix;
pub mod rand;
#[cfg(all(feature = "thiele-linear-rec", feature = "rec-from-pts"))]
pub mod rec_rat_linear;
pub mod traits;
mod util;

/// Polynomials and rational functions
pub mod algebra;

/// Function reconstruction algorithms
pub mod rec;

// #[cfg(feature = "recfn")]
// mod rec_rat_fn;
#[cfg(all(feature = "thiele-linear-rec", feature = "recfn"))]
mod rec_rat_fn_linear;

pub use ffnt::Z64;
pub use rug::Integer;
