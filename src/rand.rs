use std::hash::Hash;

use ahash::AHashSet;
use rand::{Rng, distributions::Standard, prelude::Distribution};

use crate::traits::Cardinality;

#[derive(Clone, Debug)]
pub struct UniqueRand<T> {
    seen: AHashSet<T>,
}

impl<T> Default for UniqueRand<T> {
    fn default() -> Self {
        Self { seen: Default::default() }
    }
}

impl<T: Clone + Cardinality + Hash + Eq> UniqueRand<T>
where
    Standard: Distribution<T>
{
    pub fn new() -> Self {
        Self::default()
    }

    pub fn try_gen(&mut self, mut rng: impl Rng) -> Option<T> {
        if self.seen.len() >= T::cardinality().unwrap_or(usize::MAX) {
            return None;
        }
        let mut res = rng.gen();
        while !self.seen.insert(res.clone()) {
            res = rng.gen();
        }
        Some(res)
    }
}

#[derive(Clone, Debug, Default)]
pub struct UniqueRandIter<T, R> {
    gen: UniqueRand<T>,
    rng: R
}

impl<T: Clone + Cardinality + Hash + Eq, R: Rng> UniqueRandIter<T, R> {
    pub fn new(rng: R) -> Self {
        Self { rng, gen: Default::default() }
    }
}

impl<T: Clone + Cardinality + Hash + Eq, R: Rng> Iterator for UniqueRandIter<T, R>
where
    Standard: Distribution<T>
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.gen.try_gen(&mut self.rng)
    }
}
