/// Utilities for tests and benchmarks
#[doc(hidden)]
pub mod _test_util;
mod arr;
mod matrix;
pub mod rand;
pub mod rec_linear;
pub mod rec_newton;
pub mod rec_rat;
pub mod rec_rat_linear;
pub mod rec_rat_mod;
pub mod rec_thiele;
pub mod traits;
mod util;

/// Polynomials and rational functions
pub mod algebra;

/// Function reconstruction algorithms
pub mod rec;

#[deprecated]
pub mod dense_poly;
#[deprecated]
pub mod sparse_poly;
#[deprecated]
pub mod rat;


pub use ffnt::Z64;
pub use rug::Integer;
