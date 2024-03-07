#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};
use rand::rngs::StdRng;
use rand::SeedableRng;

use merlin::Transcript;

use crate::enums::EGInp;

mod arguers;
mod traits;
mod vec_traits;
mod vec_utils;
mod mat_utils;
mod enums;
mod prover;
mod transcript;
mod mat_traits;
mod mexp_prover;
mod prod_prover;


fn main() {
    let m: u64 = 4;
    let n: u64 = 4;
    let mu: u64 = 2;

    let N = m * n;

    assert!(m % mu == 0);

    //Setup rng and common reference key(public key for ElGamal,
    // commitment key for Pedersen)
    let mut rng = StdRng::from_entropy();
    let mut cr = arguers::CommonRef::new(N, rng);

    let deck: Vec<Scalar> = (0..N).map(|card| Scalar::from(card as u64))
                                    .collect();
    let deck_r: Vec<Scalar> = (0..N).map(|_| cr.rand_scalar()).collect();
    let C_deck: Vec<Ciphertext> = deck.iter()
                                                .zip(deck_r)
                                                .map(|(card, r)| cr.encrypt(&EGInp::Scal(card.clone()), 
                                                                            &r
                                                                )
                                                )
                                            .collect();

    let permutation: Vec<u64> = cr.rand_perm(&(0..N).collect());

    let rho: Vec<Scalar> = permutation.iter()
                                            .map(|_| cr.rand_scalar())
                                            .collect();


    let C_permd: Vec<Ciphertext> = permutation.iter()
                                                .map(|pi| C_deck[*pi as usize].clone())
                                                .collect();
    let C_: Vec<Ciphertext> = rho.iter()
                                    .map(|r| cr.encrypt(&EGInp::Scal(Scalar::from(1 as u128)), r))
                                    .collect();

    let mut prover_transcript = Transcript::new(b"ShuffleProof");
    let mut shuffle_prover = prover::ShuffleProver::new(
                            m,
                            n,
                            mu,
                            C_deck,
                            C_,
                            permutation,
                            rho,
                            cr
                        );
                                
    shuffle_prover.prove(&mut prover_transcript);

    //let q: RistrettoPoint = cr.
}
