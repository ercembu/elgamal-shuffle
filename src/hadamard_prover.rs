#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};
use std::iter;

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::traits::MultiscalarMul;
use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;
use crate::errors::ProofError;

use crate::traits::{Hadamard, EGMult, InnerProduct, Multiplicat,
                    Addition};
use crate::mat_traits::MatTraits;

use crate::zero_prover::{ZeroProof, ZeroProver};

#[derive(Clone)]
pub struct HadamProof{
    c_Bi: Vec<RistrettoPoint>,
    c_Di: Vec<RistrettoPoint>,
    c_D: RistrettoPoint,
    c_1: RistrettoPoint,
    zero_proof: ZeroProof,
}

#[derive(Clone)]
pub struct HadamProver {
    c_A: Vec<RistrettoPoint>,
    c_b: RistrettoPoint,
    A: Vec<Vec<Scalar>>,
    s: Scalar,
    r: Vec<Scalar>,
    com_ref: CommonRef,
}

impl HadamProver {
    pub fn new(
        c_A: Vec<RistrettoPoint>,
        c_b: RistrettoPoint,
        A: Vec<Vec<Scalar>>,
        s: Scalar,
        r: Vec<Scalar>,
        com_ref: CommonRef
    ) -> Self {
        Self {
            c_A: c_A,
            c_b: c_b,
            A: A,
            s: s,
            r: r,
            com_ref: com_ref,
        }
    }

    pub fn exp_dot(
        a: Vec<Scalar>,
        d: Vec<Scalar>,
        y: Scalar
    ) -> Scalar {
        assert!(a.len() == d.len());
        (0..a.len())
            .fold(
                Scalar::zero(),
                |acc, j|
                acc + (a[j] * d[j] * y.pow(j as u64))
                )
    }

    pub fn prove(
        &mut self,
        trans: &mut Transcript,
    ) -> HadamProof {

        let m = self.A.len();
        let n = self.A[0].len();


        let b: Vec<Scalar> = self.A.iter()
            .fold(
                vec![Scalar::from(1 as u128); n],
                |acc, a_i|
                acc.hadamard(a_i.clone())
                );
        let B: Vec<Vec<Scalar>> = (0..m).map(
            |i| {
                match i {
                    i if i == m-1 => b.clone(),
                    _ => (0..i).fold(
                            vec![Scalar::from(1 as u128); n],
                            |acc, j|
                            acc.hadamard(self.A[j].clone())
                        ),
                    }
            }
            ).collect();

        let s_vec: Vec<Scalar> = (1..m-1).map(|_| self.com_ref.rand_scalar())
            .collect();

        let s_vec: Vec<Scalar> = iter::once(self.r[0]).chain(s_vec.into_iter())
            .chain(iter::once(self.s)).collect();

        let c_B: Vec<RistrettoPoint> = (0..B.len()).map(
            |i|
            match i {
                0 => self.c_A[0].clone(),
                i if i == m-1 => self.c_b.clone(),
                _ => self.com_ref.commit(B[i].clone(), s_vec[i].clone()),
            }
        ).collect();

        let x = trans.challenge_scalar(b"x");
        let y = trans.challenge_scalar(b"y");

        let x_pow: Vec<Scalar> = (0..c_B.len()).map(|i| x.pow(i as u64)).collect();
        let c_Di: Vec<RistrettoPoint> = (0..c_B.len()).map(
            |i|
            c_B[i] * x_pow[i]
        ).collect();

        let c_D: RistrettoPoint = RistrettoPoint::multiscalar_mul(&x_pow[0..m-1],
                                                                  &c_Di[1..m]);

        let c_1: RistrettoPoint = self.com_ref.commit(vec![-Scalar::from(1 as u128); m], 
                                      Scalar::zero());

        let D: Vec<Vec<Scalar>> = (0..m).map(
            |i|
            B[i].mult(&x.pow((i+1) as u64))
        ).collect();

        let t_: Vec<Scalar> = (0..m).map(
            |i|
            s_vec[i] * x.pow((i+1) as u64)
        ).collect();

        let d: Vec<Scalar> = (0..m-1).fold(
            vec![Scalar::zero(); B[0].len()],
            |acc, i|
            acc.add(&B[i+1].mult(&x.pow((i+1) as u64)))
            );

        let t: Scalar = (0..m-1).fold(
            Scalar::zero(),
            |acc, i|
            acc + x.pow((i+1) as u64) * s_vec[i+1]
            );

        let mut zero_prover: ZeroProver = ZeroProver::new(
            [&self.c_A[1..m], &[c_1].as_slice()].concat(),
            [&c_Di[0..m-1], &[c_D].as_slice()].concat(),
            HadamProver::exp_dot,
            self.A.clone(),
            self.r.clone(),
            D.clone(),
            t_.clone(),
            self.com_ref.clone());

        let zero_proof: ZeroProof = zero_prover.prove(trans);
        
        HadamProof {
            c_Bi: vec![],
            c_Di: vec![],
            c_D: c_D.clone(),
            c_1: RistrettoPoint::random(&mut self.com_ref.rng),
            zero_proof: zero_proof,
        }
    }
}

impl HadamProof {
    pub fn verify(
        &mut self,
        trans: &mut Transcript,
        com_ref: &mut CommonRef,
    ) -> Result<(), ProofError> {
        self.zero_proof.verify(trans, com_ref)?;
        Ok(())
    }
}
