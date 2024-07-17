#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::{RistrettoPoint, CompressedRistretto};
use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;

use crate::traits::{EGMult, InnerProduct};
use crate::mat_traits::MatTraits;
use crate::vec_utils::VecUtil::scalar_to_str;

use bulletproofs::ProofError;
use crate::mexp_prover::{MexpProver, MexpProof};
use crate::prod_prover::{ProdProver, ProdProof};

#[derive(Clone)]
pub struct ShuffleProof {
    pub(crate) c_A : Vec<RistrettoPoint>,
    pub(crate) c_B : Vec<RistrettoPoint>,
    pub(crate) mexp: MexpProof,
    pub(crate) mexp_prover: MexpProver,
    pub(crate) prod: ProdProof,
    pub(crate) prod_prover: ProdProver,
}

///Prover struct for Shuffle Argument
#[derive(Clone)]
pub struct ShuffleProver {
    /// N = m * n, cards
    m: u64,
    n: u64,
    /// mu for factorization of m
    mu: u64,
    ///Unshuffled Deck
    c_deck: Vec<Ciphertext>,
    ///Shuffled(C_pi) Deck
    cp_deck: Vec<Ciphertext>,
    /// Permutation pi for the shuffle
    pi: Vec<u64>,
    /// Permutation blinding factors
    rho: Vec<Scalar>,
    /// Common reference key
    com_ref: CommonRef,
}


impl ShuffleProver {

    pub fn new(
        m: u64,
        n: u64,
        mu: u64,
        c_deck: Vec<Ciphertext>,
        cp_deck: Vec<Ciphertext>,
        pi: Vec<u64>,
        rho: Vec<Scalar>,
        com_ref: CommonRef,
    ) -> Self {
        Self {
            m: m,
            n: n,
            mu: mu,
            c_deck: c_deck,
            cp_deck: cp_deck,
            pi: pi,
            rho: rho,
            com_ref: com_ref,
        }
    }

    pub fn prove(&mut self, trans: &mut Transcript) -> ShuffleProof
    {
        trans.shuffle_domain_sep(self.n, self.m);
        //Prover
        //Commit permutation
        let r: Vec<Scalar> = (0..self.m)
            .map(|_| self.com_ref.rand_scalar())
            .collect();
        let a: Vec<Scalar> = self.pi.iter()
            .map(|value| Scalar::from(*value as u128))
            .collect();

        let c_a: Vec<RistrettoPoint> = self.com_ref.commit_vec(a.clone(), r.clone());

        trans.append_point_vec(b"c_A", &c_a.iter().map(|point| point.compress()).collect::<Vec<CompressedRistretto>>());

        //Challenge x
        //let x = Scalar::from(2 as u128);;
        let x = trans.challenge_scalar(b"x");

        //Commit exp permutation
        let s: Vec<Scalar> = (0..self.m)
            .map(|_| self.com_ref.rand_scalar())
            .collect();

        let b: Vec<Scalar> = self.pi.iter()
            .map(|value| x.pow(*value))
            .collect();

        let c_b: Vec<RistrettoPoint> = self.com_ref.commit_vec(b.clone(), s.clone());
        trans.append_point_vec(b"c_B", &c_b.iter().map(|point| point.compress()).collect::<Vec<CompressedRistretto>>());

        //Multi-Expo Argument
        let rho_: Scalar = -self.rho.dot(&b);
        let x_: Vec<Scalar> = (1..=self.m*self.n).map(|e| x.pow(e as u64)).collect();
        //TODO: check correctness here
        let C_x: Ciphertext = self.c_deck.as_slice().pow(x_.as_slice());
        assert!(self.cp_deck.len() == (self.m * self.n) as usize);
        let mut cp_iter = self.cp_deck.clone().into_iter();
        let C_mat: Vec<Vec<Ciphertext>> =  (0..self.m).map(|_| (0..self.n).map(|_| cp_iter.next().unwrap())
                                                           .collect::<Vec<Ciphertext>>()
                                                            ).collect();

        let mut b_iter = b.clone().into_iter();
        let b_mat: Vec<Vec<Scalar>> =  (0..self.m).map(|_| (0..self.n).map(|_| b_iter.next().unwrap())
                                                           .collect::<Vec<Scalar>>()
                                                            ).collect();
        let mut mexp_prover = MexpProver::new(C_mat, C_x, c_b.clone(), b_mat, s.clone(), rho_, self.com_ref.clone());
        let mexp_proof = mexp_prover.prove(trans, x.clone());

        //Challenge y, z
        let y = trans.challenge_scalar(b"y");
        let z = trans.challenge_scalar(b"z");

        let _z: Vec<Scalar> = (0..self.pi.len()).map(|_| -z.clone()).collect();
        let zeros: Vec<Scalar> = (0..self.pi.len()).map(|_| Scalar::from(0 as u128)).collect();
        let c_z: Vec<RistrettoPoint> = self.com_ref.commit_vec(_z.clone(), zeros);

        let c_d: Vec<RistrettoPoint> = c_a.iter()
            .zip(c_b.iter())
            .map(|(a, b)| a* y + b)
            .collect();

        //Product Argument
        let d: Vec<Scalar> = a.iter()
            .zip(b.iter())
            .map(|(a_, b_)| a_* y + b_)
            .collect();
        let t: Vec<Scalar> = r.iter()
            .zip(s.iter())
            .map(|(r_, s_)| r_* y + s_)
            .collect();

        //Multiply cD * cZ
        let cd_cz: Vec<RistrettoPoint> = c_d.iter()
            .zip(c_z.iter())
            .map(|(d_, z_)| d_ + z_)
            .collect();//cd * c-z
        let d_z: Vec<Scalar> = d.iter()
            .zip(_z.iter())
            .map(|(d_, z_)| d_ + z_)
            .collect();
        let mut d_z_iter = d_z.clone().into_iter();
        let d_z: Vec<Vec<Scalar>> = (0..self.m).map(|_| 
                                                (0..self.n).map(|_| d_z_iter.next().unwrap()).collect::<Vec<Scalar>>()
                                                )
            .collect();
        //(yi + xi − z)

        let product: Scalar = (1..=self.m * self.n).fold(
            Scalar::from(1 as u128),
            |acc, i| 
            acc * (y * Scalar::from(i as u128) + (0..i).fold(
                    Scalar::from(1 as u128), 
                    |acc, _| 
                    acc * x) 
                - z)
            );
        let mut prod_prover = ProdProver::new(cd_cz, d_z, t, product, self.com_ref.clone());
        let prod_proof = prod_prover.prove(trans, x.clone());



        ShuffleProof {
            c_A: c_a,
            c_B: c_b,
            mexp: mexp_proof,
            mexp_prover: mexp_prover,
            prod: prod_proof,
            prod_prover: prod_prover,
        }
    }

    pub fn verify(
        &mut self,
        trans: &mut Transcript,
        mut proof: ShuffleProof,
    ) -> Result<(), ProofError> {
        trans.shuffle_domain_sep(self.n, self.m);

        trans.val_append_point_vec(b"c_A", &proof.c_A.iter().map(|point| point.compress()).collect::<Vec<CompressedRistretto>>())?;

        trans.val_append_point_vec(b"c_B", &proof.c_B.iter().map(|point| point.compress()).collect::<Vec<CompressedRistretto>>())?;

        let mexp_proof = proof.mexp;
        let mut mexp_prover = proof.mexp_prover;
        //mexp_prover.verify(mexp_proof, trans)?;
        //
        assert!(proof.prod.verify(trans, &mut self.com_ref.clone()).is_ok());

        Ok(())
    }


}
