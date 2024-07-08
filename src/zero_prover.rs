#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};
use std::iter;

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::traits::MultiscalarMul;
use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;

use crate::traits::{Hadamard, EGMult, InnerProduct, Multiplicat,
                    Addition};
use crate::mat_traits::MatTraits;

#[derive(Clone)]
pub struct ZeroProof {
    C: RistrettoPoint,
}

#[derive(Clone)]
pub struct ZeroProver {
    c_Ai: Vec<RistrettoPoint>,
    c_Bi: Vec<RistrettoPoint>,
    map: fn(Vec<Scalar>, Vec<Scalar>, Scalar)->Scalar,
    A: Vec<Vec<Scalar>>,
    r: Vec<Scalar>,
    B: Vec<Vec<Scalar>>,
    s: Vec<Scalar>,
    com_ref: CommonRef
}

impl ZeroProver {
    pub fn new(
        c_Ai: Vec<RistrettoPoint>,
        c_Bi: Vec<RistrettoPoint>,
        map: fn(Vec<Scalar>, Vec<Scalar>, Scalar)->Scalar,
        A: Vec<Vec<Scalar>>,
        r: Vec<Scalar>,
        B: Vec<Vec<Scalar>>,
        s: Vec<Scalar>,
        com_ref: CommonRef,
    ) -> Self {
        Self {
            c_Ai: c_Ai,
            c_Bi: c_Bi,
            map: map,
            A: A,
            r: r,
            B: B,
            s: s,
            com_ref: com_ref,
        }
    }

    pub fn prove(
        &mut self, 
        trans: &mut Transcript
    ) -> ZeroProof {
        ZeroProof{
            C: RistrettoPoint::random(&mut self.com_ref.rng)
        }
    }
}

