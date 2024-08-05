#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext, IsIdentity};
use std::iter;

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::traits::MultiscalarMul;
use merlin::Transcript;

use crate::arguers::CommonRef;
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
pub(crate) struct SVProof {
    c_d: RistrettoPoint,
    c_sig: RistrettoPoint,
    c_del: RistrettoPoint,
    a_vec: Vec<Scalar>,
    b_vec: Vec<Scalar>,
    r: Scalar,
    s: Scalar,
    x: Scalar,
}

#[derive(Clone)]
pub struct SVProver {
    c_A: RistrettoPoint,
    a_vec: Vec<Scalar>,
    b: Scalar,
    r: Scalar,
    com_ref: CommonRef,
    pub(crate) chall: Challenges
}

impl Timeable for SVProver {
}

impl SVProver {
    pub fn new(
        c_A: RistrettoPoint,
        a_vec: Vec<Scalar>,
        b: Scalar,
        r: Scalar,
        com_ref: CommonRef,
    ) -> Self {
        Self {
            c_A: c_A,
            a_vec: a_vec,
            b: b,
            r: r,
            com_ref: com_ref,
            chall: Challenges{x:Scalar::one(), y:Scalar::one(), z:Scalar::one()},
        }
    }

    pub fn prove(
        &mut self,
        trans: &mut Transcript
    ) -> SVProof {
        trans.append_message(b"dom-sep", b"SVProof");
        let n = self.a_vec.len();

        let b_vec: Vec<Scalar> = (0..n).map(|i|
                                            (0..=i).fold(
                                                Scalar::one(),
                                                |acc, j|
                                                acc * self.a_vec[j].clone())
                                            ).collect();

        let d_vec: Vec<Scalar> = vec![self.com_ref.rand_scalar(); n];
        let r_d: Scalar = self.com_ref.rand_scalar();

        let sig_vec: Vec<Scalar> = (0..n).map(|i|
                                              match i {
                                                0 => d_vec[0].clone(),
                                                i if i == (n-1) => Scalar::zero(),
                                                _ => self.com_ref.rand_scalar(),
                                              }
                                            ).collect();

        let (s_1, s_x) = (self.com_ref.rand_scalar(), self.com_ref.rand_scalar());

        let c_d = self.com_ref.commit(d_vec.clone(), r_d.clone());

        let open_sigma: Vec<Scalar> = (0..n-1).map(|i|
                                                   -sig_vec[i].clone() 
                                                   * d_vec[i+1].clone()
                                                ).collect();

        let c_sig = self.com_ref.commit(open_sigma, s_1.clone());

        let open_delta: Vec<Scalar> = (0..n-1).map(|i|
                                                   sig_vec[i+1].clone()
                                                   - self.a_vec[i+1].clone()
                                                   * sig_vec[i].clone()
                                                   - b_vec[i].clone()
                                                   * d_vec[i+1].clone()
                                                ).collect();
        let c_del = self.com_ref.commit(open_delta, s_x.clone());

        let x = self.chall.x.clone();

        let a_vec: Vec<Scalar> = (0..n).map(|i|
                                            self.a_vec[i] * x.clone() 
                                            + d_vec[i].clone()
                                        ).collect();

        let r: Scalar = x.clone() * self.r + r_d;

        let b_vec: Vec<Scalar> = (0..n).map(|i|
                                            b_vec[i] * x.clone() 
                                            + sig_vec[i].clone()
                                        ).collect();

        let s: Scalar = x.clone() * s_x + s_1;

        SVProof{
            c_d: c_d,
            c_sig: c_sig,
            c_del: c_del,
            a_vec: a_vec,
            b_vec: b_vec,
            r: r,
            s: s,
            x: x,
            
        }
    }

    pub fn verify(
        &mut self,
        trans: &mut Transcript,
        proof: SVProof,
    ) -> Result<(), ProofError> {
        trans.append_message(b"dom-sep", b"SVProof");
        let n = self.a_vec.len();

        let x = proof.x;
        assert!((self.c_A * x) + proof.c_d == self.com_ref.commit(proof.a_vec.clone(),
                                                               proof.r));

        let open_ds: Vec<Scalar> = (0..n-1).map(|i|
                                 proof.b_vec[i+1] * x
                                 - proof.b_vec[i] * proof.a_vec[i+1]
                                ).collect();
        assert!(proof.c_del * x + proof.c_sig == self.com_ref.commit(open_ds, proof.s));

        assert!(proof.b_vec[0] == proof.a_vec[0]);
        assert!(proof.b_vec[n-1] == x*self.b);

        Ok(())
        
    }
}
#[test]
fn test_sv_hiding_a() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;


    let mut rng = StdRng::seed_from_u64(2);//from_entropy();
    let m: usize = 4;
    let n: usize = 4;
    let mut com_ref = CommonRef::new((n*m) as u64, rng);

    let x = Scalar::one();
    let a: Vec<Scalar> = vec![com_ref.rand_scalar(); n];
    let r: Scalar = Scalar::one();

    let c_a: RistrettoPoint = com_ref.commit(a.clone(), r.clone());

    let d_vec: Vec<Scalar> = vec![Scalar::zero(); n];
    let r_d: Scalar = Scalar::zero();

    let c_d = com_ref.commit(d_vec.clone(), r_d.clone());

    let a_vec: Vec<Scalar> = (0..n).map(|i|
                                        a[i] * x.clone() 
                                        + d_vec[i].clone()
                                    ).collect();

    let r_: Scalar = x.clone() * r + r_d;

    assert!(c_a*x + c_d == com_ref.commit(a_vec.clone(), r_.clone()));

}

#[test]
fn test_sv_base() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    let mut prover_transcript = Transcript::new(b"testSVProof");

    let mut rng = StdRng::seed_from_u64(2);//from_entropy();
    let m: usize = 4;
    let n: usize = 4;
    let mut com_ref = CommonRef::new((n*m) as u64, rng);

    let a: Vec<Scalar> = vec![com_ref.rand_scalar(); n];
    let r: Scalar = com_ref.rand_scalar();

    let b: Scalar = a.iter().fold(Scalar::one(), |acc, i| acc * i.clone());

    let c_a: RistrettoPoint = com_ref.commit(a.clone(), r.clone());

    let mut sv_prover = SVProver::new(
        c_a,
        a,
        b,
        r,
        com_ref);

    let sv_proof = sv_prover.prove(&mut prover_transcript);

    let mut verifier_transcript = Transcript::new(b"testSVProof");
    assert!(sv_prover.verify(&mut verifier_transcript, sv_proof).is_ok());
}
