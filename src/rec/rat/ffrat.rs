use ffnt::Z64;
use rug::{Integer, Rational};

use crate::{
    algebra::{
        poly::flat::{FlatMono, FlatPoly},
        rat::{NoneError, Rat},
    },
    rec::crt::{merge_crt, rat_reconstruct},
};

/// rational function over finite characteristic that does not necessarily fit in a `Z64`
#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub(crate) struct FFRat<const N: usize> {
    pub(crate) rat: Rat<FlatPoly<Integer, N>>,
    pub(crate) modulus: Integer,
}

impl<const N: usize> FFRat<N> {
    pub(crate) fn merge_crt<const P: u64>(
        &mut self,
        new_rat: Rat<FlatPoly<Z64<P>, N>>,
    ) {
        let (num, den) = std::mem::take(&mut self.rat).into_num_den();
        let (new_num, new_den) = new_rat.into_num_den();
        debug_assert_eq!(num.len(), new_num.len());
        debug_assert_eq!(den.len(), new_den.len());
        let mut num = num.into_terms();
        let mut den = den.into_terms();
        let new_num = new_num.into_terms();
        let new_den = new_den.into_terms();

        let terms = num.iter_mut().chain(den.iter_mut());
        let new_terms = new_num.into_iter().chain(new_den);
        for (term, new_term) in terms.zip(new_terms) {
            debug_assert_eq!(term.powers, new_term.powers);
            merge_crt(&mut term.coeff, new_term.coeff, &self.modulus);
        }
        let num = FlatPoly::from_raw_terms(num);
        let den = FlatPoly::from_raw_terms(den);
        self.rat = Rat::from_num_den_unchecked(num, den);
        self.modulus *= P;
    }
}

impl<const P: u64, const N: usize> From<Rat<FlatPoly<Z64<P>, N>>> for FFRat<N> {
    fn from(source: Rat<FlatPoly<Z64<P>, N>>) -> Self {
        let rat = source.into();
        Self {
            rat,
            modulus: P.into(),
        }
    }
}

impl<'a, const N: usize> TryFrom<&'a FFRat<N>> for Rat<FlatPoly<Rational, N>> {
    type Error = NoneError;

    fn try_from(source: &'a FFRat<N>) -> Result<Self, Self::Error> {
        let mut num = Vec::with_capacity(source.rat.num().len());
        for term in source.rat.num().terms() {
            let rat_coeff = rat_reconstruct(&term.coeff, &source.modulus)
                .ok_or(NoneError {})?;
            num.push(FlatMono::new(rat_coeff, term.powers));
        }
        let num = FlatPoly::from_raw_terms(num);

        let mut den = Vec::with_capacity(source.rat.den().len());
        for term in source.rat.den().terms().iter() {
            let rat_coeff = rat_reconstruct(&term.coeff, &source.modulus)
                .ok_or(NoneError {})?;
            den.push(FlatMono::new(rat_coeff, term.powers));
        }
        let den = FlatPoly::from_raw_terms(den);

        Ok(Rat::from_num_den_unchecked(num, den))
    }
}

impl<'a, const N: usize> TryFrom<&'a FFRat<N>> for Rat<FlatPoly<Integer, N>> {
    type Error = NoneError;

    fn try_from(source: &'a FFRat<N>) -> Result<Self, Self::Error> {
        Rat::<FlatPoly<Rational, N>>::try_from(source).map(|r| r.into())
    }
}
