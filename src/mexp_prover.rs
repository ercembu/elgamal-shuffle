#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::RistrettoPoint;
use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;

use crate::traits::{EGMult, InnerProduct};
use crate::mat_traits::MatTraits;
///Prover struct for Multi-Expo Argument
pub struct MexpProver<'a> {
    
    ///m x n Cipher Matrix for C_m vectors
    C_mat: Vec<Vec<&'a Ciphertext>>,

    /// Result of the argued Multi-Exp Ciphertext
    C: Ciphertext,

    ///commitments to the open permutation
    c_A: &'a Vec<RistrettoPoint>,

    /// Scalar matrix for open permutations
    A: Vec<Vec<&'a Scalar>>,

    ///Blinding factor for open permutations
    r: &'a Vec<Scalar>,

    /// ElGamal blinding scalar rho
    rho: Scalar,

    /// Common reference key
    com_ref: &'a CommonRef,
}

impl<'a> MexpProver<'a> {
    pub fn new(
        C_mat: Vec<Vec<&'a Ciphertext>>,
        C: Ciphertext,
        c_A: &'a Vec<RistrettoPoint>,
        A: Vec<Vec<&'a Scalar>>,
        r: &'a Vec<Scalar>,
        rho: Scalar,
        com_ref: &'a CommonRef,
    ) -> Self {
        Self {    
            C_mat: C_mat,
            C: C,
            c_A: c_A,
            A: A,
            r: r,
            rho: rho,
            com_ref: com_ref,
        }
    }
    ///C_x = ElG(1;rho)C'^b
    ///c_b = com(b, s)
    pub(crate) fn prove(
        &mut self,
        trans: &mut Transcript,
    ) {

        //Format to matrix m x n: m * n = N
        println!("C MATRIX: {:?}", self.C_mat.size());
        println!("A MATRIX: {:?}", self.A.size());

    }
}
