# rare

Implementations of a number of algorithms for the black-box
reconstruction of rational functions. The task is to reconstruct a
rational function `f(x, y, z, ...)` from a number of probes `(x, y, z,
..., f(x, y, z, ...))`.

This crate is highly experimental. The current API is not very user
friendly and likely to change drastically in the future.

# Usage

Add the following lines to your Cargo.toml:

```toml
[dependencies.rare]
git = "https://github.com/a-maier/rare.git"
features = ["thiele-multivar"]
tag = "0.9.2"
```

The `features` entry determines the available reconstruction
algorithms, see the [algorithms section](#Implemented-algorithms).

# Example

Reconstruct the function `f(x) = x/(1 + x + x^2)` from its function
values. This example uses [rand](https://crates.io/crates/rand) to
generate random probes.

```
use rand::Rng;
use rare::{
   rec::{
      primes::LARGE_PRIMES,
      rat::thiele_univar::{Rec, Status}
   },
   Z64
};

let mut rec = Rec::new(1);
let mut rng = rand::thread_rng();

'rec: {
    loop {
        const P: u64 = LARGE_PRIMES[0];
        const ONE: Z64<P> = Z64::new(1);

        let z: Z64<P> = rng.gen();
        let q_z = z / (ONE + z * (ONE + z));
        rec.add_pt([z], q_z).unwrap();
        match rec.status() {
            Status::NeedNextPt => {},
            Status::NeedNextMod => break,
            Status::Done => break 'rec,
        }
    }

    // repeat for `LARGE_PRIMES[1]` etc. until done
    loop {
        const P: u64 = LARGE_PRIMES[1];
        const ONE: Z64<P> = Z64::new(1);

        let z: Z64<P> = rng.gen();
        let q_z = z / (ONE + z * (ONE + z));
        rec.add_pt([z], q_z).unwrap();
        match rec.status() {
            Status::NeedNextPt => {},
            Status::NeedNextMod => break,
            Status::Done => break 'rec,
        }
    }
}

assert_eq!(rec.status(), Status::Done);
let rec = rec.into_rat().unwrap();
eprintln!("reconstructed {rec}");
```

# Implemented algorithms

- Thiele interpolation (only univariate rational functions).
  This algorithm is implemented in [rare::rec::rat::thiele_univar].
- The algorithm by Cuyt and Lee [1] as described in the Firefly
  paper [2] using Newton interpolation for polynomials. Requires the
  `cuyt-lee-rec` feature and is implemented in
  [rare::rec::rat::cuyt_lee].
- Thiele interpolation to determine the degrees followed by solving
  linear systems of equations. Requires the `thiele-linear-rec`
  feature and is implemented in [rare::rec::rat::thiele_linear].
- Multivariate reconstruction via rescaled univariate Thiele
  interpolation. Requires the `thiele-multivar` feature and is
  implemented in [rare::rec::rat::thiele_multivar].

[1]: A. A. M. Cuyt and W. Lee, Sparse interpolation of multivariate rational functions, Theor. Comput. Sci. 412 (2011) 1445.

[2]: J. Klappert and F. Lange, Reconstructing rational functions with FireFly, Comput. Phys. Commun. 247 (2020) 106951 [arXiv:1904.00009](https://arxiv.org/abs/1904.00009).
