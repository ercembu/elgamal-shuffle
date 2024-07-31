#![allow(non_snake_case)]
#![feature(iter_next_chunk)]
use rust_elgamal::{Scalar, Ciphertext};
use rand::rngs::StdRng;
use rand::SeedableRng;

use merlin::Transcript;

use crate::utils::enums::EGInp;

mod arguers;
mod traits;
mod utils;
mod provers;

use crate::traits::*;
use crate::utils::*;
use crate::provers::*;


fn main() {
    let m: usize = 4;
    let n: usize = 4;
    let mu: usize = 2;

    let N = m * n;

    assert!(m % mu == 0);

    //Setup rng and common reference key(public key for ElGamal,
    // commitment key for Pedersen)
    let mut rng = StdRng::from_entropy();
    let mut cr = arguers::CommonRef::new(N as u64, rng);

    let deck: Vec<Scalar> = (0..N).map(|card| Scalar::from(card as u64))
                                    .collect();
    let deck_r: Vec<Scalar> = (0..N).map(|_| cr.rand_scalar()).collect();
    let C_deck: Vec<Ciphertext> = deck.iter()
                                                .zip(deck_r.clone())
                                                .map(|(card, r)| cr.encrypt(&EGInp::Scal(card.clone()), 
                                                                            &r
                                                                )
                                                )
                                            .collect();

    let permutation: Vec<u64> = cr.rand_perm(&(1..=(N as u64)).collect());

    let rho: Vec<Scalar> = (0..N).map(|_| cr.rand_scalar()).collect();

    let C_permd: Vec<Ciphertext> = permutation.iter()
                                                .zip(rho.iter())
                                                .map(|(pi, r)| {
                                                    let card: Ciphertext = C_deck[*pi as usize - 1].clone();
                                                    card + cr.encrypt(&EGInp::Scal(Scalar::zero()), &r)

                                                    }
                                                )
                                                .collect();

    let mut prover_transcript = Transcript::new(b"ShuffleProof");
    let mut shuffle_prover = prover::ShuffleProver::new(
                            m,
                            n,
                            mu,
                            C_deck,
                            C_permd,
                            permutation,
                            rho,
                            cr
                        );
                                
    let mut shuffle_proof = shuffle_prover.prove(&mut prover_transcript);

    let mut verifier_transcript = Transcript::new(b"ShuffleProof");

    assert!(shuffle_prover
            .verify(&mut verifier_transcript, shuffle_proof)
            .is_ok());
    
}
