#![allow(non_snake_case)]
use std::mem;

use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::RistrettoPoint;
use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::provers::{hadamard_prover::{HadamProof, HadamProver},
                        sv_prover::{SVProof, SVProver},
                        zero_prover::{ZeroProver}};

use crate::traits::{traits::{Hadamard, 
                                Timeable,
                                EGMult, 
                                InnerProduct, 
                                Multiplicat,
                                Addition
                            }, 
                    mat_traits::MatTraits};

use crate::utils::{utils::Challenges,
                    transcript::TranscriptProtocol,
                    errors::ProofError};
#[derive(Clone)]
pub struct ProdProof {
    pub(crate) c_b: RistrettoPoint,
    pub(crate) had_proof: HadamProof,
    pub(crate) sv_proof: SVProof,
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

impl Timeable for ProdProver {
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
    ) -> (ZeroProver, SVProver, HadamProver, ProdProof) {
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

        //HADAMARD ARGUMENT
        let mut hadamard_prover: HadamProver = HadamProver::new(self.c_A.clone(), c_b.clone(), self.A.clone(), self.r.clone(), a_vec.clone(),  s.clone(), self.com_ref.clone());

        hadamard_prover.chall = self.chall.clone();
        let hadamard_time = hadamard_prover.start_time();
        let (zero_prover, hadamard_proof): (ZeroProver, HadamProof) = hadamard_prover.prove(trans);
        println!("\n");
        println!("Hadamard Proof Time:\t{}", hadamard_prover.elapsed(hadamard_time));
        println!("Hadamard Proof:\t{}", mem::size_of_val(&hadamard_proof));
        //SINGLE_VALUE PRODUCT
        //
        let mut sv_prover: SVProver = SVProver::new(c_b,
                                                    a_vec,
                                                    self.b,
                                                    s,
                                                    self.com_ref.clone()
                                                    );
        sv_prover.chall = self.chall.clone();

        let sv_time = sv_prover.start_time();
        let sv_proof: SVProof = sv_prover.prove(trans);
        println!("\n");
        println!("SV Product Proof Time:\t{}", sv_prover.elapsed(sv_time));
        println!("SV Product Proof Size:\t{}", mem::size_of_val(&sv_proof));
        
        (zero_prover,
         sv_prover,
         hadamard_prover,
        ProdProof {
            c_b: RistrettoPoint::random(&mut self.com_ref.rng),
            had_proof: hadamard_proof,
            sv_proof: sv_proof,
        }
        )
    }

    pub fn verify(
        &mut self,
        trans: &mut Transcript,
        proof: ProdProof,
        mut had_prover: HadamProver,
        mut zero_prover: ZeroProver,
        mut sv_prover: SVProver
    ) -> Result<(), ProofError> {
        let verify_time = had_prover.start_time();
        had_prover.verify(trans, proof.had_proof, zero_prover)?;
        println!("\n");
        println!("Hadamard Verify Time:\t{}", had_prover.elapsed(verify_time));

        let verify_time = sv_prover.start_time();
        sv_prover.verify(trans, proof.sv_proof)?;
        println!("\n");
        println!("SV Product Verify Time:\t{}", sv_prover.elapsed(verify_time));
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
    let (mut zero_prover, mut sv_prover, mut hadam_prover, proof) = prod_prover.prove(&mut prover_transcript);

    let mut verifier_transcript = Transcript::new(b"testProdProof");

    assert!(prod_prover.verify(&mut verifier_transcript, proof, hadam_prover,
                               zero_prover, sv_prover).is_ok());

}
