#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::RistrettoPoint;
use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;
use crate::hadamard_prover::{HadamProof, HadamProver};

use crate::traits::{EGMult, InnerProduct};
use crate::mat_traits::MatTraits;
#[derive(Clone)]
pub struct ProdProof {
    pub(crate) c_b: RistrettoPoint,
    pub(crate) had_proof: HadamProof,
}
///Prover struct for Product Argument
#[derive(Clone)]
pub struct ProdProver {
    ///commitments to arguments
    c_A: Vec<RistrettoPoint>,

    /// Scalar matrix for arguments
    A: Vec<Vec<Scalar>>,

    ///Blinding factor for arguments
    r: Vec<Scalar>,

    /// Result Scalar
    b: Scalar,

    /// Common reference key
    com_ref: CommonRef,
}

impl ProdProver {
    pub fn new(
        c_A: Vec<RistrettoPoint>,
        A: Vec<Vec<Scalar>>,
        r: Vec<Scalar>,
        b: Scalar,
        com_ref: CommonRef,
    ) -> Self {
        Self { 
            c_A: c_A,
            A: A,
            r: r,
            b: b,
            com_ref: com_ref
        }
    }

    pub fn prove(
        &mut self,
        trans: &mut Transcript,
        x: Scalar,
    ) -> ProdProof {
        let s: Scalar = self.com_ref.rand_scalar();

        let m = self.A.len();
        let n = self.A[0].len();

        let a_vec: Vec<Scalar> = (0..m)
            .map(|i| {
                (0..n)
                    .fold(
                        Scalar::from(1 as u128),
                        |acc, j|
                        acc * self.A[j][i]
                        )
            })
        .collect();
        let c_b: RistrettoPoint = self.com_ref.commit(a_vec.clone(), s);

        //TODO: HADAMARD ARGUMENT
        let mut hadamard_prover: HadamProver = HadamProver::new(self.c_A.clone(), c_b.clone(), self.A.clone(), s.clone(), self.r.clone(), self.com_ref.clone());

        let hadamard_proof: HadamProof = hadamard_prover.prove(trans);
        //TODO: SINGLE_VALUE PRODUCT
        
        ProdProof {
            c_b: RistrettoPoint::random(&mut self.com_ref.rng),
            had_proof: hadamard_proof,
        }
    }

}
