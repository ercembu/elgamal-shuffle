#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT;

use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;

use crate::enums::EGInp;
use crate::errors::ProofError;

use crate::traits::{EGMult, InnerProduct, Addition, Multiplicat};
use crate::mat_traits::MatTraits;
pub struct MexpProof {
    pub(crate) c_A0: RistrettoPoint,
    pub(crate) c_Bk: Vec<RistrettoPoint>,
    pub(crate) Ek  : Vec<Ciphertext>,
    pub(crate) a_  : Vec<Scalar>,
    pub(crate) r   : Scalar,
    pub(crate) b   : Scalar,
    pub(crate) s   : Scalar,
    pub(crate) tau : Scalar,
}
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
    com_ref: &'a mut CommonRef,
}

impl<'a> MexpProver<'a> {
    pub fn new(
        C_mat: Vec<Vec<&'a Ciphertext>>,
        C: Ciphertext,
        c_A: &'a Vec<RistrettoPoint>,
        A: Vec<Vec<&'a Scalar>>,
        r: &'a Vec<Scalar>,
        rho: Scalar,
        com_ref: &'a mut CommonRef,
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
    ) -> MexpProof {

        //Format to matrix m x n: m * n = N
        let (m, n) = self.C_mat.size();

        trans.mexp_domain_sep(m.clone() as u64, (m/2).try_into().unwrap());

        println!("C MATRIX: {:?}", self.C_mat.size());
        println!("A MATRIX: {:?}", self.A.size());

        let a_0: Vec<Scalar> = (0..n).map(|_| self.com_ref.rand_scalar()).collect();
        let r_0: Scalar = self.com_ref.rand_scalar();

        let mut b_: Vec<Scalar> = vec![];
        let mut s_: Vec<Scalar> = vec![];
        let mut tau_: Vec<Scalar> = vec![];

        for k in 0..=(2*m - 1) {
            match k {
                m => {
                    b_.push(self.com_ref.rand_scalar());
                    s_.push(self.com_ref.rand_scalar());
                    tau_.push(self.com_ref.rand_scalar());
                },
                _ => {
                    b_.push(Scalar::zero());
                    s_.push(Scalar::zero());
                    tau_.push(self.rho.clone());
                }
            }
        }

        let k = 2*m;

        let c_A0: RistrettoPoint = self.com_ref.commit(a_0.clone(), r_0);

        let c_Bk: Vec<RistrettoPoint> = self.com_ref.commit_vec(b_.clone(), s_.clone());

        let mut Gbk: Vec<Ciphertext> = b_.iter()
                                        .zip(tau_.iter())
                                        .map(|(_b, _tau)| {
                                            self.com_ref.encrypt(&EGInp::Rist(RISTRETTO_BASEPOINT_POINT * _b), 
                                                                 _tau)
                                        }).collect();

        for i in 1..m {
            for j in 0..m {
                let k = m + j - i;

                Gbk[k] = Gbk[k] + self.C_mat[i].as_slice().pow(self.A[j].as_slice());
            }
        }

        let Ek: Vec<Ciphertext> = Gbk;

        println!("{:?}", Ek);
        //Send: cA0, {cBk}2m−1k=0 , {Ek}2m−1k=0
        //Challenge: x ← Z∗q.

        let x: Scalar = trans.challenge_scalar(b"x");

        let x_: Vec<Scalar> = (1..=m).map(|exp| x.pow(exp.try_into().unwrap())).collect();
        println!("{:?}", x_);

        let Ax: Vec<Scalar> = self.A.mult(&x_);
        let a_: Vec<Scalar> = a_0.add(&Ax);

        let r: Scalar = r_0 + self.r.clone().dot(&x_);

        let b: Scalar = b_[0] + (1..m*2 -1).map(|k| 
                                               b_[k] * x.pow(k.try_into().unwrap()))
                                            .reduce(|acc, bx| acc + bx).unwrap();

        let s: Scalar = s_[0] + (1..m*2 -1).map(|k| 
                                               s_[k] * x.pow(k.try_into().unwrap()))
                                            .reduce(|acc, sx| acc + sx).unwrap();
        let tau: Scalar = tau_[0] + (1..m*2 -1).map(|k| 
                                               tau_[k] * x.pow(k.try_into().unwrap()))
                                            .reduce(|acc, taux| acc + taux).unwrap();
        //Send: a_, r, b, s, tau
        MexpProof{
            c_A0: c_A0,
            c_Bk: c_Bk,
            Ek  : Ek,
            a_  : a_,
            r   : r,
            b   : b,
            s   : s,
            tau : tau,
        }


    }

    pub fn verify(
        &mut self,
        proof: MexpProof,
        trans: &mut Transcript,
    ) -> Result<(), ProofError> {
        Ok(())
    }

}
