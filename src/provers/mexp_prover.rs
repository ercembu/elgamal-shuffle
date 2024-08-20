//! Struct Proof and Prover for Multi-Exponentiation Argument
#![allow(non_snake_case)]
use std::mem;
use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT;

use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::traits::{traits::{Hadamard, 
                                HeapSize,
                                EasySize,
                                EGMult, 
                                InnerProduct, 
                                Multiplicat,
                                Addition,
                                Timeable
                            }, 
                    mat_traits::MatTraits};

use crate::utils::{utils::Challenges,
                    transcript::TranscriptProtocol,
                    errors::ProofError,
                    enums::EGInp};

///Data struct for the non optimized mexp arguments
#[derive(Clone, Default)]
pub struct MexpProof {
    ///Commitment to the beggining blinding vector
    pub(crate) c_A0: RistrettoPoint,
    ///Commitment to the blinding values
    pub(crate) c_Bk: Vec<RistrettoPoint>,
    ///Message products of the blinded diagonals
    pub(crate) Ek  : Vec<Ciphertext>,
    ///Blinded column vectors
    pub(crate) a_  : Vec<Scalar>,
    ///Blinding value for the column vectors
    pub(crate) r   : Scalar,
    ///X Challenged blinding values
    pub(crate) b   : Scalar,
    ///X Challenged blinding values for b values
    pub(crate) s   : Scalar,
    ///X Challenged ElGamal Blinding Values
    pub(crate) tau : Scalar,
}

impl HeapSize for MexpProof {
    fn heap_size(&self) -> usize {
        self.c_A0.ez_size()
            + self.c_Bk.ez_size()
            + self.Ek.ez_size()
            + self.a_.ez_size()
            + self.r.ez_size()
            +self.b.ez_size()
            +self.s.ez_size()
            +self.tau.ez_size()
    }

}

///Data Struct for Optimized final Mexp Proof
///Optimization yields square root N smaller proof size while taking square root N more time
#[derive(Clone)]
pub struct MexpOptimProof {
    ///Committed Blinding values for the ElGamal messages
    pub(crate) c_b: Vec<RistrettoPoint>,
    ///Result of the Blinded Diagonal Product
    pub(crate) E_k: Vec<Ciphertext>,
    ///Open Blinding values for elGamal messages
    pub(crate) b: Scalar,
    ///Hiding value for the commitment of the blinding values
    pub(crate) s: Scalar,
    ///X Blinded Messages
    pub(crate) open_C: Ciphertext,
    ///Diagonal Message challenged by x and hidden via b values
    pub(crate) C_: Ciphertext,
}
impl HeapSize for MexpOptimProof {
    fn heap_size(&self) -> usize {
        self.c_b.ez_size()
            + self.E_k.ez_size()
            +self.b.ez_size()
            +self.s.ez_size()
            +self.open_C.ez_size()
            +self.C_.ez_size()
    }

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

    ///Global Challenges
    pub(crate) chall: Challenges,
}

impl Timeable for MexpProver {
}

impl MexpProver {
    ///Base Constructor
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

    ///Prover method for a non-optimized proof
    pub(crate) fn prove(
        &mut self,
        trans: &mut Transcript,
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

        for i in 0..m {
            for j in 0..=m {
                let k = m + j - i - 1;

                Gbk[k] = Gbk[k] + self.C_mat[i].as_slice().pow(
                    if j == 0 {a_0.as_slice()} else {self.A[j-1].as_slice()}
                    );
            }
        }

        let Ek: Vec<Ciphertext> = Gbk;



        let x = self.chall.x.clone();
        let x_: Vec<Scalar> = (1..=m).map(|exp| x.pow(exp.try_into().unwrap())).collect();

        let Ax: Vec<Scalar> = self.A.mult(&x_);
        let a_: Vec<Scalar> = a_0.add(&Ax);

        let r: Scalar = r_0 + self.r.clone().dot(&x_);

        let b: Scalar = (1..=2*m-1).fold(b_[0], |acc, k| acc + b_[k] * x.pow(k as u64));
        let s: Scalar = (1..=2*m-1).fold(s_[0], |acc, k| acc + s_[k] * x.pow(k as u64));
        let tau: Scalar = (1..=2*m-1).fold(tau_[0], |acc, k| acc + tau_[k] * x.pow(k as u64));
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

    ///Verifier method for a non optimized proof
    pub fn verify(
        &mut self,
        proof: MexpProof,
        trans: &mut Transcript,
    ) -> Result<(), ProofError> {
        let (m, n) = self.C_mat.size();
        trans.mexp_domain_sep(m.clone() as u64, (m/2).try_into().unwrap());

        assert!(proof.c_Bk[m] == self.com_ref.commit(vec![Scalar::zero()], Scalar::zero()));
        assert!(proof.Ek[m] == self.C);

        let x_: Vec<Scalar> = (1..=m).map(|exp| self.chall.x.pow(exp.try_into().unwrap())).collect();

        let open_A: RistrettoPoint = self.c_A.iter()
            .zip(x_.iter())
            .fold(proof.c_A0, |acc, (a, x)| acc + (a * x));

        let commit_A: RistrettoPoint = self.com_ref.commit(proof.a_.clone(),
                                                            proof.r.clone());

        let x = self.chall.x.clone();

        assert!(open_A == commit_A);

        let open_B: RistrettoPoint = (1..=2*m-1).fold(proof.c_Bk[0], |acc, i|
                                                    acc + (proof.c_Bk[i] 
                                                    * x.pow((i) as u64)));
        let commit_B: RistrettoPoint = self.com_ref.commit(vec![proof.b.clone()],
                                                           proof.s.clone());

        assert!(open_B == commit_B);

        let open_E: Ciphertext = (1..=2*m-1).fold(proof.Ek[0], |acc, k|
                                                  acc + proof.Ek[k] * x.pow(k as u64));

        let commit_E: Ciphertext = (0..m).fold(
            self.com_ref.encrypt(&EGInp::Scal(proof.b), &proof.tau),
            |acc, i|
            acc + (&self.C_mat[i][..]).pow(&proof.a_.mult(&x.pow((m-i-1) as u64))[..])
        );

        assert!(open_E == commit_E);

        Ok(())
    }

    ///Prover method for an optimized Mexp Argument
    pub fn prove_optim(
        &mut self,
        trans: &mut Transcript,
        mu: usize,
    ) -> MexpOptimProof {
        let (m, n) = self.C_mat.size();
        let m_: usize = m / mu;

        let x = self.chall.x.clone();

        trans.mexp_domain_sep(m.clone() as u64, (m/2).try_into().unwrap());

        //Create blinding vectors
        let mut b_: Vec<Scalar> = vec![self.com_ref.rand_scalar(); 2*mu -1];
        let mut s_: Vec<Scalar> = vec![self.com_ref.rand_scalar(); 2*mu -1];
        let mut tau_: Vec<Scalar> = vec![self.com_ref.rand_scalar(); 2*mu -1];

        b_[mu-1] = Scalar::zero();
        s_[mu-1] = Scalar::zero();
        tau_[mu-1] = self.rho.clone();

        //Commitments to b vector
        let c_bk: Vec<RistrettoPoint> = (0..=2*mu -2).map(|k|
                                                          self.com_ref.commit(vec![b_[k].clone()], s_[k].clone())
                                                          ).collect();
        //Base Ciphertexts for diagonal products
        let mut Gbk: Vec<Ciphertext> = b_.iter()
                                        .zip(tau_.iter())
                                        .map(|(_b, _tau)| {
                                            self.com_ref.encrypt(&EGInp::Scal(*_b) 
                                                                  ,_tau)
                                        }).collect();

        //Diagonal Operations
        for i in 1..=mu {
            for j in 1..=mu {
                let k = j + mu - i - 1 ;
                Gbk[k] = (0..=m_ - 1).fold(Gbk[k].clone(),
                                            |acc, l|
                                                acc + self.C_mat[mu * l + i - 1].as_slice()
                                                                .pow(self.A[mu * l + j - 1].as_slice())
                                           );
            }
        }

        let E_k = Gbk;
        


        //Exponentiatied challenges
        let x_: Vec<Scalar> = (0..=2*mu -2).map(|exp| x.pow(exp.try_into().unwrap())).collect();

        //Challenged blinding values
        let b: Scalar = b_.dot(&x_.clone());
        let s: Scalar = s_.dot(&x_.clone());
        let rho_: Scalar = tau_.dot(&x_.clone());

        //Challenging the initial arguments
        let mut a_: Vec<Vec<Scalar>> = vec![];
        let mut r_: Vec<Scalar> = vec![];

        for l in 1..=m_ {
            a_.push((1..=mu).fold(vec![Scalar::zero(); self.A[0].len()],
                                    |acc, j|
                                    acc.add(&self.A[mu * (l - 1) + j - 1]
                                                            .mult(&x.pow((j-1) as u64))
                                            )
                                    ));
            r_.push((1..=mu).fold(Scalar::zero(),
                                    |acc, j|
                                    acc + (&self.r[mu * (l - 1) + j - 1]
                                                            * x.pow((j-1) as u64)
                                            )
                                    ));
        }

        let mut C_l: Vec<Vec<Ciphertext>> = vec![];
        let mut c_A_prime: Vec<RistrettoPoint> = vec![];

        //Base Message for diagonal products
        let G_b: Ciphertext = self.com_ref.encrypt(&EGInp::Scal(-b.clone()) 
                                                              ,&Scalar::zero()
                                                  );
        //Challenged Diagonal Messages
        let E_x: Ciphertext = E_k.clone().as_slice().pow(x_.clone().as_slice());

        let C_prime: Ciphertext = G_b + E_x;


        //Various challengings
        for l in 1..=m_ {
            C_l.push(
                (1..mu).fold(self.C_mat[mu*l -1].clone(),
                    |acc, i| {
                        let rhs: Vec<Ciphertext> = self.C_mat[mu*(l - 1) + i - 1]
                            .iter()
                            .zip(acc.iter())
                            .map(|(c, a)| a + c * x.pow((mu - i) as u64))
                            .collect();
                        rhs
                    })
                );
            c_A_prime.push(
                (2..=mu).fold(self.c_A[mu*l -1].clone(),
                    |acc, j| {
                        let rhs: RistrettoPoint = self.c_A[mu*(l - 1) + j - 1]
                            * x.pow((j-1) as u64);
                        acc + rhs
                    })
                );
        }
        let base: Ciphertext = self.com_ref.encrypt(&EGInp::Scal(Scalar::zero()) 
                                                              ,&rho_.clone()
                                                  );

        //Mexp of the initial arguments
        let open_C: Ciphertext = C_l.iter()
            .zip(a_.iter())
            .fold(base, 
                  |acc, (c, a)| acc + c.as_slice().pow(a.as_slice())
                  );

        MexpOptimProof {
            c_b: c_bk,
            E_k: E_k,
            b: b,
            s: s,
            open_C: open_C,
            C_: C_prime,
        }
    }
    
    ///Verifier for the optimized Mexp Proof
    pub fn verify_optim(
        &mut self,
        proof: MexpOptimProof,
        trans: &mut Transcript,
        mu: usize,
    ) -> Result<(), ProofError> {
        assert!(proof.c_b[mu-1] == self.com_ref.commit(vec![Scalar::zero()],
                                                        Scalar::zero())
                );

        assert!(proof.E_k[mu-1] == self.C);

        let x = self.chall.x.clone();
        let x_: Vec<Scalar> = (0..=2*mu -2).map(|exp| x.pow(exp.try_into().unwrap())).collect();
        let commit_b: RistrettoPoint = proof.c_b.iter()
            .zip(x_.iter())
            .map(|(b_, x_)| b_ * x_).sum();
        assert!(commit_b == self.com_ref.commit(vec![proof.b], proof.s));

        assert!(proof.C_ == proof.open_C);
        Ok(())
    }
}

#[test]
fn test_mexp_base_obs() {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use crate::utils::enums::EGInp;
    use std::time::SystemTime;

    let now = SystemTime::now();
    
    let mut rng = StdRng::seed_from_u64(2);//from_entropy();
    let m: usize = 8;
    let n: usize = 8;

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

    let x: Scalar = cr.rand_scalar();
    mexp_prover.chall.x = x;

    let now = mexp_prover.start_time();
    let mut prover_transcript = Transcript::new(b"testMexpProof");
    let mexp_proof = mexp_prover.prove(&mut prover_transcript);
    let mut verifier_transcript = Transcript::new(b"testMexpProof");

    println!("Base Mexp Proof Size:\t{}", &mexp_proof.heap_size());


    match now.elapsed() {
       Ok(elapsed) => {
           // it prints '2'
           println!("Base Mexp Proof Time:\t{}", elapsed.as_millis());
       }
       Err(e) => {
           // an error occurred!
           println!("Error: {e:?}");
       }
   }
    let now = mexp_prover.start_time();
    assert!(mexp_prover.verify(mexp_proof, &mut verifier_transcript).is_ok());
    match now.elapsed() {
       Ok(elapsed) => {
           // it prints '2'
           println!("Base Mexp Verify Time:\t{}", elapsed.as_millis());
       }
       Err(e) => {
           // an error occurred!
           println!("Error: {e:?}");
       }
   }

}
#[test]
fn test_mexp_optim() {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use crate::utils::enums::EGInp;
    use std::time::SystemTime;

    let now = SystemTime::now();
    
    let mut rng = StdRng::seed_from_u64(2);//from_entropy();
    let m: usize = 16;
    let n: usize = 8;
    let mu: usize = m / 1;

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

    let x: Scalar = cr.rand_scalar();
    mexp_prover.chall.x = x.clone();
    let mut prover_transcript = Transcript::new(b"testMexpProof");
    let mexp_proof = mexp_prover.prove_optim(&mut prover_transcript,
                                                mu);
    let mut verifier_transcript = Transcript::new(b"testMexpProof");

    assert!(mexp_prover.verify_optim(mexp_proof, &mut verifier_transcript, mu).is_ok());

    match now.elapsed() {
       Ok(elapsed) => {
           println!("Optimized, Mexp(mu={}): {}", mu, elapsed.as_millis());
       }
       Err(e) => {
           println!("Error: {e:?}");
       }
   }
}
