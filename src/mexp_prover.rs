#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT;

use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;

use crate::enums::EGInp;

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
    ) {

        //Format to matrix m x n: m * n = N
        let (m, n) = self.C_mat.size();
        println!("C MATRIX: {:?}", self.C_mat.size());
        println!("A MATRIX: {:?}", self.A.size());

        let a_0: Vec<Scalar> = (0..n).map(|_| self.com_ref.rand_scalar()).collect();
        let r_0: Scalar = self.com_ref.rand_scalar();

        let mut b: Vec<Scalar> = vec![];
        let mut s: Vec<Scalar> = vec![];
        let mut tau: Vec<Scalar> = vec![];

        for k in 0..=(2*m - 1) {
            match k {
                m => {
                    b.push(self.com_ref.rand_scalar());
                    s.push(self.com_ref.rand_scalar());
                    tau.push(self.com_ref.rand_scalar());
                },
                _ => {
                    b.push(Scalar::zero());
                    s.push(Scalar::zero());
                    tau.push(self.rho.clone());
                }
            }
        }

        let k = 2*m;

        let c_A0: RistrettoPoint = self.com_ref.commit(a_0, r_0);

        let c_Bk: Vec<RistrettoPoint> = self.com_ref.commit_vec(b.clone(), s);

        let Gbk: Vec<Ciphertext> = b.iter()
                                        .zip(tau.iter())
                                        .map(|(b_, tau_)| {
                                            self.com_ref.encrypt(&EGInp::Rist(RISTRETTO_BASEPOINT_POINT * b_), 
                                                                 tau_)
                                        }).collect();

        let mut E_k: Vec<Ciphertext> = Vec::with_capacity(k);
        unsafe {E_k.set_len(k);}

        for i in 1..m {
            for j in 0..m {
                let k = m + j - i;

                E_k[k] = Gbk[k] + self.C_mat[i].as_slice().pow(self.A[j].as_slice())
            }
        }
        println!("{:?}", E_k);


    }

}
