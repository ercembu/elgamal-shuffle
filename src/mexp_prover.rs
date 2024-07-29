#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT;

use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;

use crate::enums::EGInp;
use crate::utils::Challenges;
use bulletproofs::ProofError;

use crate::traits::{EGMult, InnerProduct, Addition, Multiplicat};
use crate::mat_traits::MatTraits;
#[derive(Clone, Default)]
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
#[derive(Clone)]
pub struct MexpProver {
    
    ///m x n Cipher Matrix for C_m vectors
    C_mat: Vec<Vec<Ciphertext>>,

    /// Result of the argued Multi-Exp Ciphertext
    C: Ciphertext,

    ///commitments to the open permutation
    c_A: Vec<RistrettoPoint>,

    /// Scalar matrix for open permutations
    A: Vec<Vec<Scalar>>,

    ///Blinding factor for open permutations
    r: Vec<Scalar>,

    /// ElGamal blinding scalar rho
    rho: Scalar,

    /// Common reference key
    com_ref: CommonRef,

    chall: Challenges,
}

impl MexpProver {
    pub fn new(
        C_mat: Vec<Vec<Ciphertext>>,
        C: Ciphertext,
        c_A: Vec<RistrettoPoint>,
        A: Vec<Vec<Scalar>>,
        r: Vec<Scalar>,
        rho: Scalar,
        com_ref: CommonRef,
    ) -> Self {
        Self {    
            C_mat: C_mat,
            C: C,
            c_A: c_A,
            A: A,
            r: r,
            rho: rho,
            com_ref: com_ref,
            chall: Challenges::default(),
        }
    }
    ///C_x = ElG(1;rho)C'^b
    ///c_b = com(b, s)
    pub(crate) fn prove(
        &mut self,
        trans: &mut Transcript,
        x: Scalar,
    ) -> MexpProof {

        //Format to matrix m x n: m * n = N
        let (m, n) = self.C_mat.size();

        trans.mexp_domain_sep(m.clone() as u64, (m/2).try_into().unwrap());

        let a_0: Vec<Scalar> = (0..n).map(|_| self.com_ref.rand_scalar()).collect();
        let r_0: Scalar = self.com_ref.rand_scalar();

        let mut b_: Vec<Scalar> = vec![];
        let mut s_: Vec<Scalar> = vec![];
        let mut tau_: Vec<Scalar> = vec![];

        for k in 0..=(2*m - 1) {
            match k==m {
                true => {
                    b_.push(Scalar::zero());
                    s_.push(Scalar::zero());
                    tau_.push(self.rho.clone());
                }
                false => {
                    b_.push(self.com_ref.rand_scalar());
                    s_.push(self.com_ref.rand_scalar());
                    tau_.push(self.com_ref.rand_scalar());
                }
            }
        }


        let c_A0: RistrettoPoint = self.com_ref.commit(a_0.clone(), r_0);

        let c_Bk: Vec<RistrettoPoint> = b_.clone().into_iter()
                                            .zip(s_.clone().into_iter())
                                            .map(|(b_k, s_k)| self.com_ref.commit(vec![b_k], s_k))
                                            .collect();

        let mut Gbk: Vec<Ciphertext> = b_.iter()
                                        .zip(tau_.iter())
                                        .map(|(_b, _tau)| {
                                            self.com_ref.encrypt(&EGInp::Scal(*_b) 
                                                                  ,_tau)
                                        }).collect();
        println!("{:#?}", Gbk[m]);
        println!("{:#?}", self.com_ref.encrypt(&EGInp::Scal(Scalar::zero()),
                            &self.rho.clone()));

        for i in 0..m {
            for j in 0..=m {
                let k = m + j - i - 1;

                Gbk[k] = Gbk[k] + self.C_mat[i].as_slice().pow(
                    if j == 0 {a_0.as_slice()} else {self.A[j-1].as_slice()}
                    );
            }
        }

        let Ek: Vec<Ciphertext> = Gbk;

        //Send: cA0, {cBk}2m−1k=0 , {Ek}2m−1k=0
        trans.append_point(b"c_A0", &c_A0.compress());
        trans.append_point_vec(b"c_Bk", &c_Bk.iter()
                                        .map(|p| p.compress()).collect());
        trans.append_cipher_vec(b"Ek", &Ek);
        //Challenge: x ← Z∗q.


        let x_: Vec<Scalar> = (1..=m).map(|exp| x.pow(exp.try_into().unwrap())).collect();

        let Ax: Vec<Scalar> = self.A.to_col().mult(&x_);
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
        trans.append_scalar_vec(b"a_", &a_);
        trans.append_scalar(b"r", &r);
        trans.append_scalar(b"b", &b);
        trans.append_scalar(b"s", &s);
        trans.append_scalar(b"tau", &tau);

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
        let (m, n) = self.C_mat.size();
        trans.mexp_domain_sep(m.clone() as u64, (m/2).try_into().unwrap());
        //TODO:maybe not needed
        //trans.val_append_point_vec(b"c_Bk", &proof.c_Bk.iter()
        //                                .map(|p| p.compress()).collect())?;

        assert!(proof.c_Bk[m] == self.com_ref.commit(vec![Scalar::zero()], Scalar::zero()));
        //TODO: problematic, check for encryption correctness
        println!("{}", proof.Ek.len());
        assert!(proof.Ek[m] == self.C);

        trans.append_scalar_vec(b"a_", &proof.a_);
        trans.append_scalar(b"r", &proof.r);
        trans.append_scalar(b"b", &proof.b);
        trans.append_scalar(b"s", &proof.s);
        trans.append_scalar(b"tau", &proof.tau);
        Ok(())
    }
}

#[test]
fn test_mexp_base() {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use crate::EGInp;
    
    let mut rng = StdRng::seed_from_u64(2);//from_entropy();
    let m: usize = 6;
    let n: usize = 4;

    let mut cr = CommonRef::new((m*n) as u64, rng);

    let deck: Vec<Scalar> = (0..(m*n)).map(|card| Scalar::from(card as u64))
                                    .collect();
    let deck_r: Vec<Scalar> = (0..(m*n)).map(|_| cr.rand_scalar()).collect();
    let C_deck: Vec<Ciphertext> = deck.iter()
                                        .zip(deck_r.clone())
                                        .map(|(card, r)| cr.encrypt(&EGInp::Scal(card.clone()), 
                                                                    &r
                                                        )
                                        )
                                    .collect();
    let mut c_iter = C_deck.clone().into_iter();
    let C_mat: Vec<Vec<Ciphertext>> = (0..m).map(|_| (0..n).map(|_| c_iter.next().unwrap()).collect::<Vec<Ciphertext>>()).collect();

    let A: Vec<Vec<Scalar>> = vec![vec![cr.rand_scalar(); m];n];
    let r: Vec<Scalar> = vec![cr.rand_scalar(); m];

    let c_A = cr.commit_mat(A.clone(), r.clone());

    let A = A.to_col();

    let rho: Scalar = cr.rand_scalar();
    let base: Ciphertext = cr.encrypt(&EGInp::Scal(Scalar::zero()), &rho);

    let C: Ciphertext = (0..m).fold(base.clone(), |acc, i| acc + 
                                                    (&C_mat[i][..]).pow(&A[i][..])
                                                    );

    let mut mexp_prover = MexpProver::new(
            C_mat,
            C,
            c_A,
            A,
            r,
            rho,
            cr.clone()
        );

    let x: Scalar = Scalar::one();
    let mut prover_transcript = Transcript::new(b"testMexpProof");
    let mexp_proof = mexp_prover.prove(&mut prover_transcript, x.clone());
    let mut verifier_transcript = Transcript::new(b"testMexpProof");

    assert!(mexp_prover.verify(mexp_proof, &mut verifier_transcript).is_ok());

}
