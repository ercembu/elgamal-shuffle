#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::RistrettoPoint;
use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;

use crate::traits::{EGMult, InnerProduct};
use crate::mat_traits::MatTraits;
///Prover struct for Product Argument
pub struct ProdProver<'a> {
    ///commitments to arguments
    c_A: &'a Vec<RistrettoPoint>,

    /// Scalar matrix for arguments
    A: Vec<Vec<&'a Scalar>>,

    ///Blinding factor for arguments
    r: &'a Vec<Scalar>,

    /// Result Scalar
    b: Scalar,

    /// Common reference key
    com_ref: &'a CommonRef,
}

impl<'a> ProdProver<'a> {
    pub fn new(
        c_A: &'a Vec<RistrettoPoint>,
        A: Vec<Vec<&'a Scalar>>,
        r: &'a Vec<Scalar>,
        b: Scalar,
        com_ref: &'a CommonRef,
    ) -> Self {
        Self { 
            c_A: c_A,
            A: A,
            r: r,
            b: b,
            com_ref: com_ref
        }
    }
}
