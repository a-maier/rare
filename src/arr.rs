use std::{
    fmt::{self, Display},
    ops::{
        Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
    },
};

// array class that supports basic element-wise arithmetics
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub(crate) struct Arr<T, const N: usize>(pub [T; N]);

impl<T, const N: usize> Index<usize> for Arr<T, N> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T, const N: usize> IndexMut<usize> for Arr<T, N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<T: AddAssign, const N: usize> AddAssign for Arr<T, N> {
    fn add_assign(&mut self, rhs: Self) {
        for (lhs, rhs) in self.0.iter_mut().zip(rhs.0.into_iter()) {
            *lhs += rhs;
        }
    }
}

impl<T: AddAssign, const N: usize> Add<Arr<T, N>> for Arr<T, N> {
    type Output = Self;

    fn add(mut self, rhs: Arr<T, N>) -> Self::Output {
        self += rhs;
        self
    }
}

impl<T: SubAssign, const N: usize> SubAssign for Arr<T, N> {
    fn sub_assign(&mut self, rhs: Self) {
        for (lhs, rhs) in self.0.iter_mut().zip(rhs.0.into_iter()) {
            *lhs -= rhs;
        }
    }
}

impl<T: SubAssign, const N: usize> Sub<Arr<T, N>> for Arr<T, N> {
    type Output = Self;

    fn sub(mut self, rhs: Arr<T, N>) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<T: Neg, const N: usize> Neg for Arr<T, N> {
    type Output = Arr<<T as Neg>::Output, N>;

    fn neg(self) -> Self::Output {
        Arr(self.0.map(|n| -n))
    }
}

impl<T: Copy + MulAssign, const N: usize> MulAssign<T> for Arr<T, N> {
    fn mul_assign(&mut self, rhs: T) {
        for n in &mut self.0 {
            *n *= rhs;
        }
    }
}

impl<'a, T: MulAssign<&'a T>, const N: usize> MulAssign<&'a T> for Arr<T, N> {
    fn mul_assign(&mut self, rhs: &'a T) {
        for n in &mut self.0 {
            *n *= rhs;
        }
    }
}

impl<T: Copy + Mul, const N: usize> Mul<T> for Arr<T, N> {
    type Output = Arr<<T as Mul>::Output, N>;

    fn mul(self, rhs: T) -> Self::Output {
        Arr(self.0.map(|n| n * rhs))
    }
}

impl<'a, T: Mul<&'a T>, const N: usize> Mul<&'a T> for Arr<T, N> {
    type Output = Arr<<T as Mul<&'a T>>::Output, N>;

    fn mul(self, rhs: &'a T) -> Self::Output {
        Arr(self.0.map(|n| n * rhs))
    }
}

impl<T: Display, const N: usize> Display for Arr<T, N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[")?;
        if let Some((first, rest)) = self.0.split_first() {
            write!(f, "{first}")?;
            for r in rest {
                write!(f, ", {r}")?;
            }
        }
        write!(f, "]")
    }
}
