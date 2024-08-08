#![allow(dead_code)]
use curve25519_dalek::ristretto::RistrettoPoint;
use rust_elgamal::{Scalar, Ciphertext};
pub enum EGInp {
    Rist(RistrettoPoint),
    Scal(Scalar),
}
