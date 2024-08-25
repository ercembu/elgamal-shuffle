//! Main file exampling setup and parameters
#![allow(warnings)]
#![allow(non_snake_case)]
#![allow(warnings)]
#![feature(iter_next_chunk)]
use rust_elgamal::{Scalar, Ciphertext};
use rand_chacha::ChaCha20Rng;
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
    let mut rng = ChaCha20Rng::from_entropy();
    let mut cr = arguers::CommonRef::new(N as u64, rng);

    //Create Open Deck from scalars
    let deck: Vec<Scalar> = (0..N).map(|card| Scalar::from(card as u64))
                                    .collect();
    //Get blinding values for the deck
    let deck_r: Vec<Scalar> = (0..N).map(|_| cr.rand_scalar()).collect();
    //Encrypt the decks in ElGamal
    let C_deck: Vec<Ciphertext> = deck.iter()
                                        .zip(deck_r.clone())
                                        .map(|(card, r)| 
                                             cr.encrypt(&EGInp::Scal(card.clone()), &r)
                                        )
                                .collect();

    //Get a random permutation
    let permutation: Vec<u64> = cr.rand_perm(&(1..=(N as u64)).collect());

    //Get more blinding values for hiding the permuted cards
    let rho: Vec<Scalar> = (0..N).map(|_| cr.rand_scalar()).collect();

    //Permute the cards while hiding them with another ElGamal Encryption
    let C_permd: Vec<Ciphertext> = permutation.iter()
                                                .zip(rho.iter())
                                                .map(|(pi, r)| {
                                                    let card: Ciphertext = C_deck[*pi as usize - 1].clone();
                                                    card + cr.encrypt(&EGInp::Scal(Scalar::zero()), &r)

                                                    }
                                                )
                                                .collect();

    //Create the main shuffle prover
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
                                
    //Recieve the provers used in the argument and the final proof
    let (mut zero_prover,
         mut sv_prover,
         mut hadam_prover,
         mut prod_prover,
         mut mexp_prover,
         mut shuffle_proof) = shuffle_prover.prove(&mut prover_transcript);

    //Create a clean transcript for the verification
    let mut verifier_transcript = Transcript::new(b"ShuffleProof");

    //Assert verification returns Ok
    assert!(shuffle_prover
            .verify(&mut verifier_transcript, shuffle_proof, 
                    zero_prover,
                    sv_prover,
                    hadam_prover,
                    prod_prover,
                    mexp_prover)
            .is_ok());
    
}
