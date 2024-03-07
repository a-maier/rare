use ffnt::Z64;
use log::{debug, trace};
use paste::paste;
use rand::Rng;
use rug::Integer;
use seq_macro::seq;

use crate::{
    algebra::{poly::flat::FlatPoly, rat::Rat},
    rec::{
        primes::LARGE_PRIMES,
        rat::{
            ffrat::FFRat,
            finite::{cuyt_lee::RatRecMod, linear::RecLinear},
        },
    },
    traits::{Rec, TryEval},
};

#[allow(deprecated)]
use crate::rec_rat::{normalise_coeff, RatRec};

seq!(N in 0..114 {
    paste! { const [<P N>]: u64 = LARGE_PRIMES[N]; }
});

#[allow(deprecated)]
impl<F, const N: usize> Rec<RatRec, [Integer; N]> for F
where
    F: Rec<
        RatRecMod,
        [[Z64<P0>; N]; 1],
        Output = Option<Rat<FlatPoly<Z64<P0>, N>>>,
    >,
    F: TryEval<[Z64<P1>; N], Output = Z64<P1>>,
    F: TryEval<[Z64<P2>; N], Output = Z64<P2>>,
    F: TryEval<[Z64<P3>; N], Output = Z64<P3>>,
    F: TryEval<[Z64<P4>; N], Output = Z64<P4>>,
    F: TryEval<[Z64<P5>; N], Output = Z64<P5>>,
    F: TryEval<[Z64<P6>; N], Output = Z64<P6>>,
    F: TryEval<[Z64<P7>; N], Output = Z64<P7>>,
    F: TryEval<[Z64<P8>; N], Output = Z64<P8>>,
    F: TryEval<[Z64<P9>; N], Output = Z64<P9>>,
    F: TryEval<[Z64<P10>; N], Output = Z64<P10>>,
    F: TryEval<[Z64<P11>; N], Output = Z64<P11>>,
    F: TryEval<[Z64<P12>; N], Output = Z64<P12>>,
    F: TryEval<[Z64<P13>; N], Output = Z64<P13>>,
    F: TryEval<[Z64<P14>; N], Output = Z64<P14>>,
    F: TryEval<[Z64<P15>; N], Output = Z64<P15>>,
    F: TryEval<[Z64<P16>; N], Output = Z64<P16>>,
    F: TryEval<[Z64<P17>; N], Output = Z64<P17>>,
    F: TryEval<[Z64<P18>; N], Output = Z64<P18>>,
    F: TryEval<[Z64<P19>; N], Output = Z64<P19>>,
    F: TryEval<[Z64<P20>; N], Output = Z64<P20>>,
    F: TryEval<[Z64<P21>; N], Output = Z64<P21>>,
    F: TryEval<[Z64<P22>; N], Output = Z64<P22>>,
    F: TryEval<[Z64<P23>; N], Output = Z64<P23>>,
    F: TryEval<[Z64<P24>; N], Output = Z64<P24>>,
    F: TryEval<[Z64<P25>; N], Output = Z64<P25>>,
    F: TryEval<[Z64<P26>; N], Output = Z64<P26>>,
    F: TryEval<[Z64<P27>; N], Output = Z64<P27>>,
    F: TryEval<[Z64<P28>; N], Output = Z64<P28>>,
    F: TryEval<[Z64<P29>; N], Output = Z64<P29>>,
    F: TryEval<[Z64<P30>; N], Output = Z64<P30>>,
    F: TryEval<[Z64<P31>; N], Output = Z64<P31>>,
    F: TryEval<[Z64<P32>; N], Output = Z64<P32>>,
    F: TryEval<[Z64<P33>; N], Output = Z64<P33>>,
    F: TryEval<[Z64<P34>; N], Output = Z64<P34>>,
    F: TryEval<[Z64<P35>; N], Output = Z64<P35>>,
    F: TryEval<[Z64<P36>; N], Output = Z64<P36>>,
    F: TryEval<[Z64<P37>; N], Output = Z64<P37>>,
    F: TryEval<[Z64<P38>; N], Output = Z64<P38>>,
    F: TryEval<[Z64<P39>; N], Output = Z64<P39>>,
    F: TryEval<[Z64<P40>; N], Output = Z64<P40>>,
    F: TryEval<[Z64<P41>; N], Output = Z64<P41>>,
    F: TryEval<[Z64<P42>; N], Output = Z64<P42>>,
    F: TryEval<[Z64<P43>; N], Output = Z64<P43>>,
    F: TryEval<[Z64<P44>; N], Output = Z64<P44>>,
    F: TryEval<[Z64<P45>; N], Output = Z64<P45>>,
    F: TryEval<[Z64<P46>; N], Output = Z64<P46>>,
    F: TryEval<[Z64<P47>; N], Output = Z64<P47>>,
    F: TryEval<[Z64<P48>; N], Output = Z64<P48>>,
    F: TryEval<[Z64<P49>; N], Output = Z64<P49>>,
    F: TryEval<[Z64<P50>; N], Output = Z64<P50>>,
    F: TryEval<[Z64<P51>; N], Output = Z64<P51>>,
    F: TryEval<[Z64<P52>; N], Output = Z64<P52>>,
    F: TryEval<[Z64<P53>; N], Output = Z64<P53>>,
    F: TryEval<[Z64<P54>; N], Output = Z64<P54>>,
    F: TryEval<[Z64<P55>; N], Output = Z64<P55>>,
    F: TryEval<[Z64<P56>; N], Output = Z64<P56>>,
    F: TryEval<[Z64<P57>; N], Output = Z64<P57>>,
    F: TryEval<[Z64<P58>; N], Output = Z64<P58>>,
    F: TryEval<[Z64<P59>; N], Output = Z64<P59>>,
    F: TryEval<[Z64<P60>; N], Output = Z64<P60>>,
    F: TryEval<[Z64<P61>; N], Output = Z64<P61>>,
    F: TryEval<[Z64<P62>; N], Output = Z64<P62>>,
    F: TryEval<[Z64<P63>; N], Output = Z64<P63>>,
    F: TryEval<[Z64<P64>; N], Output = Z64<P64>>,
    F: TryEval<[Z64<P65>; N], Output = Z64<P65>>,
    F: TryEval<[Z64<P66>; N], Output = Z64<P66>>,
    F: TryEval<[Z64<P67>; N], Output = Z64<P67>>,
    F: TryEval<[Z64<P68>; N], Output = Z64<P68>>,
    F: TryEval<[Z64<P69>; N], Output = Z64<P69>>,
    F: TryEval<[Z64<P70>; N], Output = Z64<P70>>,
    F: TryEval<[Z64<P71>; N], Output = Z64<P71>>,
    F: TryEval<[Z64<P72>; N], Output = Z64<P72>>,
    F: TryEval<[Z64<P73>; N], Output = Z64<P73>>,
    F: TryEval<[Z64<P74>; N], Output = Z64<P74>>,
    F: TryEval<[Z64<P75>; N], Output = Z64<P75>>,
    F: TryEval<[Z64<P76>; N], Output = Z64<P76>>,
    F: TryEval<[Z64<P77>; N], Output = Z64<P77>>,
    F: TryEval<[Z64<P78>; N], Output = Z64<P78>>,
    F: TryEval<[Z64<P79>; N], Output = Z64<P79>>,
    F: TryEval<[Z64<P80>; N], Output = Z64<P80>>,
    F: TryEval<[Z64<P81>; N], Output = Z64<P81>>,
    F: TryEval<[Z64<P82>; N], Output = Z64<P82>>,
    F: TryEval<[Z64<P83>; N], Output = Z64<P83>>,
    F: TryEval<[Z64<P84>; N], Output = Z64<P84>>,
    F: TryEval<[Z64<P85>; N], Output = Z64<P85>>,
    F: TryEval<[Z64<P86>; N], Output = Z64<P86>>,
    F: TryEval<[Z64<P87>; N], Output = Z64<P87>>,
    F: TryEval<[Z64<P88>; N], Output = Z64<P88>>,
    F: TryEval<[Z64<P89>; N], Output = Z64<P89>>,
    F: TryEval<[Z64<P90>; N], Output = Z64<P90>>,
    F: TryEval<[Z64<P91>; N], Output = Z64<P91>>,
    F: TryEval<[Z64<P92>; N], Output = Z64<P92>>,
    F: TryEval<[Z64<P93>; N], Output = Z64<P93>>,
    F: TryEval<[Z64<P94>; N], Output = Z64<P94>>,
    F: TryEval<[Z64<P95>; N], Output = Z64<P95>>,
    F: TryEval<[Z64<P96>; N], Output = Z64<P96>>,
    F: TryEval<[Z64<P97>; N], Output = Z64<P97>>,
    F: TryEval<[Z64<P98>; N], Output = Z64<P98>>,
    F: TryEval<[Z64<P99>; N], Output = Z64<P99>>,
    F: TryEval<[Z64<P100>; N], Output = Z64<P100>>,
    F: TryEval<[Z64<P101>; N], Output = Z64<P101>>,
    F: TryEval<[Z64<P102>; N], Output = Z64<P102>>,
    F: TryEval<[Z64<P103>; N], Output = Z64<P103>>,
    F: TryEval<[Z64<P104>; N], Output = Z64<P104>>,
    F: TryEval<[Z64<P105>; N], Output = Z64<P105>>,
    F: TryEval<[Z64<P106>; N], Output = Z64<P106>>,
    F: TryEval<[Z64<P107>; N], Output = Z64<P107>>,
    F: TryEval<[Z64<P108>; N], Output = Z64<P108>>,
    F: TryEval<[Z64<P109>; N], Output = Z64<P109>>,
    F: TryEval<[Z64<P110>; N], Output = Z64<P110>>,
    F: TryEval<[Z64<P111>; N], Output = Z64<P111>>,
    F: TryEval<[Z64<P112>; N], Output = Z64<P112>>,
    F: TryEval<[Z64<P113>; N], Output = Z64<P113>>,
{
    type Output = Option<Rat<FlatPoly<Integer, N>>>;

    fn rec_with_ran(&mut self, rec: RatRec, mut rng: impl Rng) -> Self::Output {
        let rec_mod = RatRecMod::new(rec.extra_pts());

        debug!("Trying rational reconstruction over characteristic {P0}");
        let mod_rec = Rec::<RatRecMod, [[Z64<P0>; N]; 1]>::rec_with_ran(
            self, rec_mod, &mut rng,
        )?;
        let ncoeff = mod_rec.num().len() + mod_rec.den().len() - 1;
        let mod_rec = normalise_coeff(mod_rec);
        debug!("Reconstructed {mod_rec}");
        let mut mod_rec = FFRat::from(mod_rec);
        let mut res: Result<Rat<FlatPoly<Integer, N>>, _> =
            (&mod_rec).try_into();

        seq!( M in 1..114 {{
            const P: u64 = paste!{ [<P M>] };
            debug!("Trying rational reconstruction over characteristic {P}");
            // next batch of sampling points
            let mut pts = Vec::with_capacity(ncoeff);
            while pts.len() < ncoeff {
                let pt: [Z64<P>; N] = [(); N].map(|_| rng.gen());
                if let Some(val) = self.try_eval(&pt) {
                    pts.push((pt, val))
                }
            }
            // if we can already reproduce the points we are done
            if let Ok(res_ref) = res.as_ref() {
                let sample_same = pts.iter()
                    .all(|(pt, val)| res_ref.try_eval(pt) == Some(*val));

                if sample_same {
                    return Some(res.unwrap());
                }
            }

            let next_mod_rec = mod_rec.rat.rec_linear(pts)?;
            let next_mod_rec = normalise_coeff(next_mod_rec);
            debug!("Reconstructed {next_mod_rec}");
            mod_rec.merge_crt(next_mod_rec);
            res  = (&mod_rec).try_into();
        }});

        debug!("Rational reconstruction failed");
        trace!("Final value: {res:?}");
        None
    }
}
