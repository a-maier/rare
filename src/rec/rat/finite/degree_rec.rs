use std::ops::ControlFlow;

use log::debug;
use thiserror::Error;

use crate::{
    algebra::{poly::dense::DensePoly, rat::Rat},
    rec::rat::finite::thiele::ThieleRec,
    traits::{One, Zero},
    Z64,
};

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct DegreeRec<const P: u64, const N: usize> {
    extra_pts: usize,
    rec: ThieleRec<P>,
    powers: [[u32; 2]; N],
    cur_idx: usize,
    last_arg: Option<[Z64<P>; N]>,
}

impl<const P: u64, const N: usize> DegreeRec<P, N> {
    pub fn new(extra_pts: usize) -> Self {
        Self {
            extra_pts,
            powers: [[0; 2]; N],
            cur_idx: 0,
            last_arg: None,
            rec: ThieleRec::new(extra_pts),
        }
    }

    pub fn add_pt(
        &mut self,
        mut z: [Z64<P>; N],
        q_z: Z64<P>,
    ) -> Result<ControlFlow<[[u32; 2]; N], usize>, Error<P, N>> {
        let idx = self.cur_idx;
        if idx >= N {
            return Err(Error::Finished);
        }
        self.check_arg(z)?;
        match self.rec.add_pt(z[idx], q_z) {
            ControlFlow::Continue(()) => {
                z[self.cur_idx] += Z64::one();
                self.last_arg = Some(z);
                Ok(ControlFlow::Continue(self.cur_idx))
            }
            ControlFlow::Break(()) => {
                let res = std::mem::replace(
                    &mut self.rec,
                    ThieleRec::new(self.extra_pts),
                );
                let res: Rat<DensePoly<_>> = res.into_rat().into();
                let num_pow = res.num().len().try_into().unwrap();
                if num_pow == 0 {
                    assert!(self.powers.is_zero());
                    return Ok(ControlFlow::Break(self.powers));
                }
                let den_pow = res.den().len().try_into().unwrap();
                self.powers[idx] = [num_pow, den_pow];
                debug!("Powers in variable {idx}: {:?}", self.powers[idx]);
                self.cur_idx += 1;
                if self.cur_idx < N {
                    self.last_arg = None;
                    Ok(ControlFlow::Continue(self.cur_idx))
                } else {
                    Ok(ControlFlow::Break(self.powers))
                }
            }
        }
    }

    pub fn cur_powers(&self) -> [[u32; 2]; N] {
        self.powers
    }

    fn check_arg(&self, z: [Z64<P>; N]) -> Result<(), Error<P, N>> {
        if let Some(last) = self.last_arg {
            let idx = self.cur_idx;
            let ok = z
                .into_iter()
                .zip(last)
                .enumerate()
                .all(|(n, (z, last))| z == last || n == idx);
            if !ok {
                return Err(Error::BadPt { z, last, idx });
            }
        }
        Ok(())
    }
}

#[derive(Debug, Error)]
pub enum Error<const P: u64, const N: usize> {
    #[error("Cannot use point: Found {z:?}, expected a difference in coordinate {idx} from {last:?}")]
    BadPt {
        z: [Z64<P>; N],
        idx: usize,
        last: [Z64<P>; N],
    },
    #[error("Reconstruction finished, no further points accepted")]
    Finished,
}
