use thiserror::Error;

use crate::rec::primes::LARGE_PRIMES;

pub(crate) fn find_largest_missing_mod(
    mods: impl Iterator<Item = u64>
) -> Option<u64> {
    use std::collections::BTreeSet;
    let mut all_mods = BTreeSet::from_iter(LARGE_PRIMES);
    for mmod in mods {
        all_mods.remove(&mmod);
    }
    all_mods.pop_last()
}

#[derive(Copy, Clone, Debug, Error, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum RecError {
    #[error("Wrong characteristic: expected {expected}, got {found}")]
    Mod { expected: u64, found: u64 },
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ModPts<const N: usize> {
    pub modulus: u64,
    pub pts: Vec<([u64; N], u64)>,
}
