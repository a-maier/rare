use seq_macro::seq;

pub(crate) fn slice_start<T, const N: usize, const M: usize>(
    a: &[T; N]
) -> &[T; M] {
    a[..M].try_into().unwrap()
}

pub(crate) const ALL_VARS: [&str; 100] = {
    let mut vars = [""; 100];
    seq!( N in 0..100 {
        vars[N] = concat!("x", stringify!(N));
    });
    vars
};

pub(crate) const ALL_VARS_Z: [&str; 100] = {
    let mut vars = [""; 100];
    seq!( N in 0..100 {
        vars[N] = concat!("z", stringify!(N));
    });
    vars
};
