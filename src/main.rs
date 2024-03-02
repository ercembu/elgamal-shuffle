#![allow(non_snake_case)]
use rust_elgamal::{DecryptionKey, Scalar, GENERATOR_TABLE};
use rand::rngs::StdRng;
use rand::SeedableRng;

use curve25519_dalek::ristretto::RistrettoPoint;
use bulletproofs::PedersenGens;

mod arguers;
mod traits;
mod vec_traits;
mod vec_utils;
mod mat_utils;
mod enums;


fn main() {
    let m = 4;
    let n = 4;
    let mu = 2;

    let N = mu * n;

    let mut rng = StdRng::from_entropy();
    let mut cr = arguers::CommonRef::new(n, rng);


    //let q: RistrettoPoint = cr.
}
