use ffnt::Z64;
use rand::Rng;

pub(crate) fn pt_iter<const P: u64>(
    mut rng: impl Rng,
) -> impl Iterator<Item = Z64<P>> {
    let start: Z64<P> = rng.gen();
    (0..P).map(move |i| start + unsafe { Z64::new_unchecked(i) })
}
