use ffnt::Z64;

/// A probe, aka sampling point, of the function to be reconstructed
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Probe<const P: u64, const N: usize> {
    /// Function argument
    pub arg: [Z64<P>; N],
    /// Function value for the given argument
    pub val: Z64<P>,
}

impl<const P: u64, const N: usize> Default for Probe<P, N> {
    fn default() -> Self {
        Self {
            arg: [(); N].map(|_| Default::default()),
            val: Default::default(),
        }
    }
}
