#![allow(non_snake_case)]
use rust_elgamal::{DecryptionKey, Scalar, GENERATOR_TABLE, Ciphertext};
use rand::rngs::StdRng;
use rand::SeedableRng;

use curve25519_dalek::ristretto::RistrettoPoint;
use bulletproofs::PedersenGens;

use crate::enums::EGInp;

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

    let deck: Vec<Scalar> = (0..N).map(|card| Scalar::from(card as u64))
                                    .collect();
    let deck_r: Vec<Scalar> = (0..N).map(|_| cr.rand_scalar()).collect();
    let deck_ciphertext: Vec<Ciphertext> = deck.iter()
                                                .zip(deck_r)
                                                .map(|(card, r)| cr.encrypt(&EGInp::Scal(card.clone()), 
                                                                            &r
                                                                )
                                                )
                                            .collect();

    let permutation: Vec<usize> = cr.rand_perm(&(0..N).collect());

    println!("{:?}", permutation);


    //let q: RistrettoPoint = cr.
}
