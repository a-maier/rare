use std::hash::Hash;

use ahash::AHashSet;
use rand::{Rng, distributions::Standard, prelude::Distribution};

use crate::traits::Cardinality;

#[derive(Clone, Debug, Default)]
pub struct UniqueRan<T, R> {
    seen: AHashSet<T>,
    rng: R
}

impl<T: Clone + Cardinality + Hash + Eq, R: Rng> UniqueRan<T, R>
where
    Standard: Distribution<T>
{
    pub fn new(rng: R) -> Self {
        Self {rng, seen: Default::default()}
    }
}

impl<T: Clone + Cardinality + Hash + Eq, R: Rng> Iterator for UniqueRan<T, R>
where
    Standard: Distribution<T>
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.seen.len() >= T::cardinality().unwrap_or(usize::MAX) {
            return None;
        }
        let mut res = self.rng.gen();
        while !self.seen.insert(res.clone()) {
            res = self.rng.gen();
        }
        Some(res)
    }
}
