use curve25519_dalek::ristretto::RistrettoPoint;
use rust_elgamal::{Scalar};
pub enum EGInp {
    Rist(RistrettoPoint),
    Scal(Scalar),
}

