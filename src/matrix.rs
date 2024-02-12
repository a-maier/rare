// TODO: there has to be something better
//       problems with swapping rows:
//       - ndarray is beyond clumsy to use
//       - nalgebra uses column-major storage
//       - grid can't do it

use std::{
    fmt::{self, Display},
    ops::{Index, IndexMut, Mul, MulAssign, SubAssign},
    slice::{Chunks, ChunksMut, RChunks},
};

use log::trace;
use num_traits::Inv;

use crate::traits::{One, Zero};

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub(crate) struct Matrix<T> {
    nrows: usize,
    ncols: usize,
    elem: Vec<T>,
}

impl<T> Matrix<T> {
    pub(crate) fn from_vec(nrows: usize, elem: Vec<T>) -> Self {
        let ncols = if elem.is_empty() {
            0
        } else {
            elem.len() / nrows
        };
        assert_eq!(
            nrows * ncols,
            elem.len(),
            "Number of elements has to be a multiple of the number of rows"
        );
        Self { nrows, ncols, elem }
    }

    pub(crate) fn nrows(&self) -> usize {
        self.nrows
    }

    pub(crate) fn row(&self, r: usize) -> &[T] {
        let row_length = self.ncols();
        let start_idx = r * row_length;
        &self.elem[start_idx..(start_idx + row_length)]
    }

    pub(crate) fn row_mut(&mut self, r: usize) -> &mut [T] {
        let row_length = self.ncols();
        let start_idx = r * row_length;
        &mut self.elem[start_idx..(start_idx + row_length)]
    }

    pub(crate) fn rows(&self) -> Chunks<'_, T> {
        self.elem.chunks(self.ncols())
    }

    pub(crate) fn rrows(&self) -> RChunks<'_, T> {
        self.elem.rchunks(self.ncols())
    }

    pub(crate) fn rows_mut(&mut self) -> ChunksMut<'_, T> {
        let ncols = self.ncols();
        self.elem.chunks_mut(ncols)
    }

    pub(crate) fn ncols(&self) -> usize {
        self.ncols
    }

    pub(crate) fn swap_rows(&mut self, mut i: usize, mut j: usize) {
        if i == j {
            return;
        }
        if i > j {
            std::mem::swap(&mut i, &mut j);
        }
        let row_length = self.ncols();
        let elems = &mut self.elem[i * row_length..];
        let (first_row, rest) = elems.split_at_mut(row_length);
        let second_row_idx = (j - i - 1) * row_length;
        let second_row =
            &mut rest[second_row_idx..(second_row_idx + row_length)];
        first_row.swap_with_slice(second_row)
    }

    pub(crate) fn truncate_rows(&mut self, nrows: usize) {
        trace!("Truncating to {nrows} rows");
        self.elem.truncate(nrows * self.ncols());
        self.nrows = nrows;
    }
}

impl<T> Matrix<T>
where
    T: Copy
        + One
        + Zero
        + Inv<Output = T>
        + MulAssign
        + Mul<Output = T>
        + SubAssign,
{
    pub(crate) fn row_reduce(&mut self) {
        self.forward_elimination();
        self.backward_substitution();
    }

    fn forward_elimination(&mut self) {
        for nrow in 0..std::cmp::min(self.nrows(), self.ncols()) {
            let Some(pivot_col) = self.pivot(nrow) else {
                return;
            };
            let pivot =
                std::mem::replace(&mut self[(nrow, pivot_col)], T::one());
            let inv_pivot = pivot.inv();
            let row = self.row_mut(nrow);
            for e in &mut row[pivot_col + 1..] {
                *e *= inv_pivot;
            }
            for sub_row in (nrow + 1)..self.nrows() {
                let fact = std::mem::replace(
                    &mut self[(sub_row, pivot_col)],
                    T::zero(),
                );
                for col in (pivot_col + 1)..self.ncols() {
                    let sub = fact * self[(nrow, col)];
                    self[(sub_row, col)] -= sub;
                }
            }
        }
    }

    fn backward_substitution(&mut self) {
        for nrow in (1..self.nrows()).rev() {
            let Some(pivot_col) = self.first_nonzero_col(nrow, nrow) else {
                continue;
            };
            let pivot_col = nrow + pivot_col;
            debug_assert!(self[(nrow, pivot_col)].is_one());
            for sub_row in 0..nrow {
                let fact = std::mem::replace(
                    &mut self[(sub_row, pivot_col)],
                    T::zero(),
                );
                for col in (pivot_col + 1)..self.ncols() {
                    let sub = fact * self[(nrow, col)];
                    self[(sub_row, col)] -= sub;
                }
            }
        }
    }

    fn pivot(&mut self, nrow: usize) -> Option<usize> {
        if !self[(nrow, nrow)].is_zero() {
            return Some(nrow);
        }
        let (pivot_col, pivot_row) = (nrow..self.nrows())
            .filter_map(|n| self.first_nonzero_col(n, nrow).map(|c| (c, n)))
            .min()?;
        self.swap_rows(nrow, pivot_row);
        Some(pivot_col + nrow)
    }
}

impl<T: Zero> Matrix<T> {
    pub(crate) fn trim_end(&mut self) {
        let last_nonzero_row =
            self.rrows().position(|r| !r.iter().all(Zero::is_zero));
        if let Some(last_nonzero_row) = last_nonzero_row {
            let new_nrows = self.nrows - last_nonzero_row;
            self.truncate_rows(new_nrows);
        } else {
            trace!("All rows zero");
            self.elem.clear()
        }
    }

    fn first_nonzero_col(&self, nrow: usize, nskip: usize) -> Option<usize> {
        self.row(nrow).iter().skip(nskip).position(|e| !e.is_zero())
    }
}

impl<T> Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (row, col) = index;
        let row_length = self.ncols();
        &self.elem[row * row_length + col]
    }
}

impl<T> IndexMut<(usize, usize)> for Matrix<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (row, col) = index;
        let row_length = self.ncols();
        &mut self.elem[row * row_length + col]
    }
}

impl<T: Display> Display for Matrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.nrows() {
            0 => writeln!(f, "[]"),
            1 => {
                write!(f, "[")?;
                for e in &self.elem {
                    write!(f, " {e}")?;
                }
                writeln!(f, " ]")
            }
            _ => {
                for row in 0..self.nrows() {
                    write!(f, "|")?;
                    let start_idx = row * self.ncols;
                    for e in &self.elem[start_idx..(start_idx + self.ncols())] {
                        write!(f, " {e}")?;
                    }
                    writeln!(f, " |")?;
                }
                Ok(())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ffnt::Z64;

    #[test]
    fn solve_linear() {
        let _ = env_logger::builder().is_test(true).try_init();

        const P: u64 = 7;
        let one: Z64<P> = One::one();
        let mut sys = Matrix::from_vec(2, vec![one, -one, one, one, one, one]);
        sys.row_reduce();
        eprintln!("{sys}");
        // assert_eq!(x, [one, Zero::zero()]);
    }
}
