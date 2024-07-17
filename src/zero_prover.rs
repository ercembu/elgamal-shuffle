#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};
use std::iter;
use rand::rngs::StdRng;
use rand::SeedableRng;

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::traits::MultiscalarMul;
use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;
use crate::errors::ProofError;

use crate::traits::{Hadamard, EGMult, InnerProduct, Multiplicat,
                    Addition};
use crate::mat_traits::MatTraits;

#[derive(Clone)]
pub struct ZeroProof {
    c_A0: RistrettoPoint,
    c_Bm: RistrettoPoint,
    c_D: Vec<RistrettoPoint>,
    a_vec: Vec<Scalar>,
    b_vec: Vec<Scalar>,
    r: Scalar,
    s: Scalar,
    t: Scalar,
}

#[derive(Clone)]
pub struct ZeroProver {
    c_Ai: Vec<RistrettoPoint>,
    c_Bi: Vec<RistrettoPoint>,
    bi_map: fn(Vec<Scalar>, Vec<Scalar>, Scalar)->Scalar,
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
        bi_map: fn(Vec<Scalar>, Vec<Scalar>, Scalar)->Scalar,
        A: Vec<Vec<Scalar>>,
        r: Vec<Scalar>,
        B: Vec<Vec<Scalar>>,
        s: Vec<Scalar>,
        com_ref: CommonRef,
    ) -> Self {
        Self {
            c_Ai: c_Ai,
            c_Bi: c_Bi,
            bi_map: bi_map,
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
        let n: usize = self.A[0].len();
        let m: usize = self.A.len();
        let m_0: usize = self.r.len();
        let a_0: Vec<Scalar> = (0..n).map(|_| self.com_ref.rand_scalar()).collect();
        let b_m: Vec<Scalar> = (0..n).map(|_| self.com_ref.rand_scalar()).collect();

        let r_0: Scalar = self.com_ref.rand_scalar();
        let s_m: Scalar = self.com_ref.rand_scalar();

        let c_A0: RistrettoPoint = self.com_ref.commit(a_0.clone(), r_0.clone());
        let c_Bm: RistrettoPoint = self.com_ref.commit(b_m.clone(), s_m.clone());

        let y = trans.challenge_scalar(b"y");

        let blind_B: Vec<Vec<Scalar>> = [&self.B.clone()[..], 
                                            &[b_m].as_slice()]
                                        .concat();
        let blind_A: Vec<Vec<Scalar>> = [&[a_0].as_slice(),
                                            &self.A.clone()[..]]
                                        .concat();

        let blind_r: Vec<Scalar> = [&[r_0].as_slice(),
                                    &self.r.clone()[..]]
                                    .concat();
        let blind_s: Vec<Scalar> = [&self.s.clone()[..],
                                    &[s_m].as_slice()]
                                    .concat();
        let mut d_k: Vec<Scalar> = vec![Scalar::zero(); 2*m + 1];

        for i in 0..=m {
            for j in 0..=m {
                let k = m + i -j;

                d_k[k] = d_k[k] + (self.bi_map)(blind_A[i].clone(), blind_B[j].clone(), y.clone());
            }
        }

        let t : Vec<Scalar> = (0..=2*m).map(|i| match i {
            i if i == m+1 => Scalar::zero(),
            _ => self.com_ref.rand_scalar(),
        }).collect();

        println!("{:#?}", d_k[m-2]);
        println!("{:#?}", Scalar::zero());
        
        let c_D: Vec<RistrettoPoint> = self.com_ref.commit_vec(d_k.clone(), t.clone());

        let x = trans.challenge_scalar(b"x");

        let a : Vec<Scalar> = (0..=m).fold(
            vec![Scalar::zero(); n],
            |acc, i|
            acc.add(&blind_A[i].mult(&x.pow(i as u64)))
            );

        let r : Scalar = (0..=m).fold(
            Scalar::zero(),
            |acc, i|
            acc + (blind_r[i] * x.pow(i as u64))
            );

        let b: Vec<Scalar> = (0..=m).fold(
            vec![Scalar::zero(); n],
            |acc, j|
            acc.add(&blind_B[j].mult(&x.pow((m-j) as u64)))
            );
        let s: Scalar = (0..=m).fold(
            Scalar::zero(),
            |acc, j|
            acc + (blind_s[j] * x.pow((m-j) as u64))
            );

        let t_val: Scalar = (0..=2*m).fold(
            Scalar::zero(),
            |acc, k|
            acc + (t[k] * x.pow(k as u64))
            );


        
        ZeroProof{
            c_A0: c_A0,
            c_Bm: c_Bm,
            c_D: c_D,
            a_vec: a,
            b_vec: b,
            r: r,
            s: s,
            t: t_val,
        }
    }
}
impl ZeroProof{

    pub fn verify(
        &mut self,
        trans: &mut Transcript,
        com_ref: &mut CommonRef,
    ) -> Result<(), ProofError> {
        let m = (self.c_D.len() - 1)/ 2;
        let n = self.a_vec.len();

        assert!(self.c_D[m+1] == com_ref.commit(vec![Scalar::zero()], Scalar::zero()));

        Ok(())
    }
}

#[test]
fn test_base() {
    use crate::hadamard_prover::HadamProver;

    let mut prover_transcript = Transcript::new(b"testZeroProof");

    let mut rng = StdRng::from_entropy();
    let n: u64 = 4;
    let m: u64 = 2;
    let mut com_ref = CommonRef::new(n, rng);
    let a: Vec<Vec<Scalar>> = vec![vec![Scalar::zero(); 2]; 4];
    let r: Vec<Scalar> = vec![Scalar::zero(); 4];
    let b: Vec<Vec<Scalar>> = vec![vec![Scalar::one(); 2]; 4];
    let s: Vec<Scalar> = vec![Scalar::zero(); 4];


    let c_A: Vec<RistrettoPoint> = com_ref.commit_mat(a.clone(), r.clone());
    let c_B: Vec<RistrettoPoint> = com_ref.commit_mat(b.clone(), s.clone());

    let mut zero_prover: ZeroProver = ZeroProver::new(
        c_A,
        c_B,
        HadamProver::exp_dot,
        a,
        r,
        b,
        s,
        com_ref
        );

    let mut zero_proof: ZeroProof = zero_prover.prove(&mut prover_transcript);
        

}

