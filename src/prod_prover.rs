#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::RistrettoPoint;
use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;
use crate::hadamard_prover::{HadamProof, HadamProver};
use crate::sv_prover::{SVProof, SVProver};
use crate::errors::ProofError;
use crate::utils::Challenges;

use crate::traits::{EGMult, InnerProduct};
use crate::mat_traits::MatTraits;
#[derive(Clone)]
pub struct ProdProof {
    pub(crate) c_b: RistrettoPoint,
    pub(crate) had_proof: HadamProof,
    pub(crate) had_prover: HadamProver,
    pub(crate) sv_proof: SVProof,
    pub(crate) sv_prover: SVProver,
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
    
    /// Challenges from oracle, purely random
    pub(crate) chall: Challenges
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
            com_ref: com_ref,
            chall: Challenges{x:Scalar::one(), y:Scalar::one(), z:Scalar::one()},
        }
    }

    pub fn prove(
        &mut self,
        trans: &mut Transcript,
    ) -> ProdProof {
        let s: Scalar = self.com_ref.rand_scalar();

        let m = self.A.len();
        let n = self.A[0].len();

        let a_vec: Vec<Scalar> = (0..n)
            .map(|i| {
                (0..m)
                    .fold(
                        Scalar::from(1 as u128),
                        |acc, j|
                        acc * self.A[j][i]
                        )
            })
        .collect();
        let c_b: RistrettoPoint = self.com_ref.commit(a_vec.clone(), s);

        //TODO: HADAMARD ARGUMENT
        let mut hadamard_prover: HadamProver = HadamProver::new(self.c_A.clone(), c_b.clone(), self.A.clone(), self.r.clone(), a_vec.clone(),  s.clone(), self.com_ref.clone());

        hadamard_prover.chall = self.chall.clone();
        let hadamard_proof: HadamProof = hadamard_prover.prove(trans);
        //TODO: SINGLE_VALUE PRODUCT
        //
        let mut sv_prover: SVProver = SVProver::new(c_b,
                                                    a_vec,
                                                    self.b,
                                                    s,
                                                    self.com_ref.clone()
                                                    );
        sv_prover.chall = self.chall.clone();

        let sv_proof: SVProof = sv_prover.prove(trans);
        
        ProdProof {
            c_b: RistrettoPoint::random(&mut self.com_ref.rng),
            had_proof: hadamard_proof,
            had_prover: hadamard_prover,
            sv_proof: sv_proof,
            sv_prover: sv_prover,
        }
    }

    pub fn verify(
        &mut self,
        trans: &mut Transcript,
        proof: ProdProof,
    ) -> Result<(), ProofError> {
        let mut had_prover = proof.had_prover;
        had_prover.verify(trans, proof.had_proof)?;

        let mut sv_prover = proof.sv_prover;
        sv_prover.verify(trans, proof.sv_proof)?;
        Ok(())
    }
}

#[test]
fn test_prod() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    
    let mut prover_transcript = Transcript::new(b"testProdProof");

    let mut rng = StdRng::seed_from_u64(2);//from_entropy();
    let m: usize = 4;
    let n: usize = 6;
    let mut com_ref = CommonRef::new((m*n) as u64, rng);

    let a: Vec<Vec<Scalar>> = vec![vec![com_ref.rand_scalar(); m]; n];
    let r: Vec<Scalar> = vec![com_ref.rand_scalar(); m];

    let b: Scalar = a.clone().into_iter().flatten()
        .fold(Scalar::one(),
        |acc, a_ij|
        acc * a_ij
        );

    let c_A: Vec<RistrettoPoint> = com_ref.commit_mat(a.clone(), r.clone());

    let mut prod_prover = ProdProver::new(
        c_A,
        a.to_col(),
        r,
        b,
        com_ref);
    let proof = prod_prover.prove(&mut prover_transcript);

    let mut verifier_transcript = Transcript::new(b"testProdProof");

    assert!(prod_prover.verify(&mut verifier_transcript, proof).is_ok());

}
