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
        ///Prover
        ///Commit permutation
        let r: Vec<Scalar> = self.pi.iter()
            .map(|_| self.com_ref.rand_scalar())
            .collect();
        let a: Vec<Scalar> = self.pi.iter()
            .map(|value| Scalar::from(*value as u128))
            .collect();

        let c_a: Vec<RistrettoPoint> = self.com_ref.commit_vec(a.clone(), r.clone());

        ///Challenge x
        let x = trans.challenge_scalar(b"x");

        ///Prover
        ///Commit exp permutation
        let s: Vec<Scalar> = self.pi.iter()
            .map(|_| self.com_ref.rand_scalar())
            .collect();

        let b: Vec<Scalar> = self.pi.iter()
            .map(|value| x.pow(*value))
            .collect();

        let c_b: Vec<RistrettoPoint> = self.com_ref.commit_vec(b.clone(), s.clone());

        ///Challenge y, z
        let y = trans.challenge_scalar(b"y");
        let z = trans.challenge_scalar(b"z");

        let _z: Vec<Scalar> = (0..self.pi.len()).map(|_| z.clone()).collect();
        let zeros: Vec<Scalar> = (0..self.pi.len()).map(|_| Scalar::from(0 as u128)).collect();
        let c_z: Vec<RistrettoPoint> = self.com_ref.commit_vec(_z, zeros);

        let c_d: Vec<RistrettoPoint> = c_a.iter()
            .zip(c_b.iter())
            .map(|(a, b)| a* y + b)
            .collect();

        ///Product Argument
        let d: Vec<Scalar> = a.iter()
            .zip(b.iter())
            .map(|(a_, b_)| a_* y + b_)
            .collect();
        let t: Vec<Scalar> = r.iter()
            .zip(s.iter())
            .map(|(r_, s_)| r_* y + s_)
            .collect();

        self.prod_prove(trans, d, t, z, c_d, c_z);

        ///Multi-Expo Argument
        let rho_: Vec<Scalar> = self.rho.iter()
            .zip(b.iter())
            .map(|(r, b_)| -r * b_)
            .collect();

        let x_: Vec<Scalar> = (1..=self.m*self.n).map(|e| x.pow(e as usize)).collect();
        self.mexp_prove(trans, b, s, rho_, x_, c_b);
    }

    fn prod_prove(
        &mut self,
        trans: &mut Transcript,
        d: Vec<Scalar>,
        t: Vec<Scalar>,
        z: Scalar,
        c_d: Vec<RistrettoPoint>,
        c_z: Vec<RistrettoPoint>,
    ) {
        ;
        
    }

    fn mexp_prove(
        &mut self,
        trans: &mut Transcript,
        b: Vec<Scalar>,
        s: Vec<Scalar>,
        rho: Vec<Scalar>,
        x: Vec<Scalar>,
        c_b: Vec<RistrettoPoint>
    ) {
        ;
    }


}
