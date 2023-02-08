// TODO: there has to be something better
//       problems with swapping rows:
//       - ndarray is beyond clumsy to use
//       - nalgebra uses column-major storage
//       - grid can't do it

use std::ops::{Index, IndexMut};

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub(crate) struct Matrix<T> {
    nrows: usize,
    ncols: usize,
    elem: Vec<T>
}

impl<T> Matrix<T> {
    pub(crate) fn from_vec(nrows: usize, elem: Vec<T>) -> Self {
        assert_eq!(elem.len() % nrows, 0);
        let ncols = elem.len() / nrows;
        Self {
            nrows,
            ncols,
            elem
        }
    }

    pub(crate) fn nrows(&self) -> usize {
        self.nrows
    }

    pub(crate) fn row_mut(&mut self, r: usize) -> &mut [T] {
        let row_length = self.ncols();
        let start_idx = r * row_length;
        &mut self.elem[start_idx..(start_idx + row_length)]
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
        let mut second_row = &mut rest[second_row_idx..(second_row_idx + row_length)];
        first_row.swap_with_slice(&mut second_row)
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
