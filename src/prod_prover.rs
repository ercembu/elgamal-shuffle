#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::RistrettoPoint;
use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;

use crate::traits::{EGMult, InnerProduct};
use crate::mat_traits::MatTraits;
#[derive(Clone)]
pub struct ProdProof {
    pub(crate) c_b: RistrettoPoint,
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
        
        ProdProof {
            c_b: RistrettoPoint::random(&mut self.com_ref.rng),    
        }
    }

}
