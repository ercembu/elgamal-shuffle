pub mod cipher {
    use curve25519_dalek::ristretto::{RistrettoPoint, CompressedRistretto};
    use rust_elgamal::{Scalar, Ciphertext};

    pub fn cipher2bytes(
        text: &Ciphertext,
    ) -> Vec<u8> {
        let tuple = text.inner();
        [tuple.0.compress().as_bytes().as_slice(), tuple.1.compress().as_bytes().as_slice()].concat()
    }
}

use rust_elgamal::{Scalar};
#[derive(Clone)]
pub struct Challenges {
    pub x: Scalar,
    pub y: Scalar,
    pub z: Scalar,
}
