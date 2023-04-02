macro_rules! count {
    () => {0usize};
    ($_head:tt $($tail:tt)*) => {1usize + count!($($tail)*)};
}

pub(crate) use count;
