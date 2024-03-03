#![allow(non_snake_case)]
use rust_elgamal::{DecryptionKey, Scalar, GENERATOR_TABLE, Ciphertext};
use rand::rngs::StdRng;
use rand::SeedableRng;

use curve25519_dalek::ristretto::RistrettoPoint;
use bulletproofs::PedersenGens;
use merlin::Transcript;

use crate::enums::EGInp;
use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;

use crate::traits::EGMult;

///Prover struct for Shuffle Argument
pub struct ShuffleProver {
    /// N = m * n, cards
    m: usize,
    n: usize,
    ///Unshuffled Deck
    c_deck: Vec<Ciphertext>,
    ///Shuffled(C_pi) Deck
    cp_deck: Vec<Ciphertext>,
    /// Permutation pi for the shuffle
    pi: Vec<usize>,
    /// Permutation blinding factors
    rho: Vec<Scalar>,
    /// Common reference key
    com_ref: CommonRef,
}

impl ShuffleProver {

    pub fn new(
        m: usize,
        n: usize,
        c_deck: Vec<Ciphertext>,
        cp_deck: Vec<Ciphertext>,
        pi: Vec<usize>,
        rho: Vec<Scalar>,
        com_ref: CommonRef,
    ) -> Self {
        Self {
            m: m,
            n: n,
            c_deck: c_deck,
            cp_deck: cp_deck,
            pi: pi,
            rho: rho,
            com_ref: com_ref,
        }
    }

    pub fn prove(&mut self, trans: &mut Transcript) 
    {
        let r: Vec<Scalar> = self.pi.iter()
            .map(|_| self.com_ref.rand_scalar())
            .collect();
        let a: Vec<Scalar> = self.pi.iter()
            .map(|value| Scalar::from(*value as u128))
            .collect();

        let c_a: Vec<RistrettoPoint> = self.com_ref.commit_vec(a.clone(), r.clone());

        let x = trans.challenge_scalar(b"x");
        let s: Vec<Scalar> = self.pi.iter()
            .map(|_| self.com_ref.rand_scalar())
            .collect();

        let b: Vec<Scalar> = self.pi.iter()
            .map(|value| x.pow(*value))
            .collect();
    }

}
