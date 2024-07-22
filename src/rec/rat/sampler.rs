use log::{debug, trace};

use crate::traits::TryEval;

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub(crate) struct Sampler {
    successes_needed: usize,
    status: Status,
}

impl Sampler {
    pub fn new(successes_needed: usize) -> Self {
        Self {
            successes_needed,
            status: Default::default()
        }
    }

    pub fn add_pt<Arg, Val, F>(
        &mut self,
        z: &Arg,
        q_z: Val,
        rat: &F
    ) -> Status
    where
        F: TryEval<Arg, Output = Val>,
        Val: Eq,
    {
        use Status::*;
        if self.status == Failed {
            return Failed;
        }
        if rat.try_eval(z) != Some(q_z) {
            debug!("Sampling failed");
            self.status = Failed;
        } else if let Success(n) = self.status {
            self.status = if n >= self.successes_needed {
                debug!("Complete");
                Complete
            } else {
                trace!("Sampling success ({n})");
                Success(n + 1)
            };
        }
        self.status
    }

    pub(crate) fn status(&self) -> Status {
        self.status
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub(crate) enum Status {
    Complete,
    Success(usize),
    Failed,
}

impl Default for Status {
    fn default() -> Self {
        Self::Success(0)
    }
}
