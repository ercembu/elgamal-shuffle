#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::{RistrettoPoint, CompressedRistretto};
use merlin::Transcript;

use crate::arguers::CommonRef;
use crate::transcript::TranscriptProtocol;

use crate::traits::{EGMult, InnerProduct};
use crate::mat_traits::MatTraits;
use crate::vec_utils::VecUtil::scalar_to_str;

use bulletproofs::ProofError;
use crate::mexp_prover::{MexpProver, MexpProof};
use crate::prod_prover::{ProdProver, ProdProof};
use crate::utils::Challenges;

#[derive(Clone)]
pub struct ShuffleProof {
    pub(crate) c_A : Vec<RistrettoPoint>,
    pub(crate) c_B : Vec<RistrettoPoint>,
    pub(crate) mexp: MexpProof,
    pub(crate) mexp_prover: MexpProver,
    pub(crate) prod: ProdProof,
    pub(crate) prod_prover: ProdProver,
    pub(crate) y: Scalar,
}

///Prover struct for Shuffle Argument
#[derive(Clone)]
pub struct ShuffleProver {
    /// N = m * n, cards
    m: usize,
    n: usize,
    /// mu for factorization of m
    mu: usize,
    ///Unshuffled Deck
    c_deck: Vec<Ciphertext>,
    ///Shuffled(C_pi) Deck
    cp_deck: Vec<Ciphertext>,
    /// Permutation pi for the shuffle
    pi: Vec<u64>,
    /// Permutation blinding factors
    rho: Vec<Scalar>,
    /// Common reference key
    com_ref: CommonRef,
    /// Challenges from oracle, purely random
    chall: Challenges
}


impl ShuffleProver {

    pub fn new(
        m: usize,
        n: usize,
        mu: usize,
        c_deck: Vec<Ciphertext>,
        cp_deck: Vec<Ciphertext>,
        pi: Vec<u64>,
        rho: Vec<Scalar>,
        com_ref: CommonRef,
    ) -> Self {
        Self {
            m: m,
            n: n,
            mu: mu,
            c_deck: c_deck,
            cp_deck: cp_deck,
            pi: pi,
            rho: rho,
            com_ref: com_ref,
            chall: Challenges{x:Scalar::one(), y:Scalar::one(), z:Scalar::one()},
        }
    }

    pub fn prove(&mut self, trans: &mut Transcript) -> ShuffleProof
    {
        trans.shuffle_domain_sep(self.n as u64, self.m as u64);
        //Prover
        //Commit permutation
        let r: Vec<Scalar> = (0..self.m)
            .map(|_| self.com_ref.rand_scalar())
            .collect();
        let a: Vec<Scalar> = self.pi.iter()
            .map(|value| Scalar::from(value.clone() as u128))
            .collect();

        let c_a: Vec<RistrettoPoint> = self.com_ref.commit_vec(a.clone(), r.clone());


        //Challenge x
        //let x = Scalar::from(2 as u128);;
        let x = trans.challenge_scalar(b"x");

        //Commit exp permutation
        let s: Vec<Scalar> = (0..self.m)
            .map(|_| self.com_ref.rand_scalar())
            .collect();

        let b: Vec<Scalar> = self.pi.iter()
            .map(|value| x.pow(value.clone()))
            .collect();

        let c_b: Vec<RistrettoPoint> = self.com_ref.commit_vec(b.clone(), s.clone());

        //Multi-Expo Argument
        let rho_: Scalar = -self.rho.dot(&b);
        let x_: Vec<Scalar> = (1..=self.m*self.n).map(|e| x.pow(e as u64)).collect();
        //TODO: check correctness here
        let C_x: Ciphertext = self.c_deck.as_slice().pow(x_.as_slice());
        assert!(self.cp_deck.len() == (self.m * self.n) as usize);
        let mut cp_iter = self.cp_deck.clone().into_iter();
        let C_mat: Vec<Vec<Ciphertext>> =  (0..self.m).map(|_| (0..self.n).map(|_| cp_iter.next().unwrap())
                                                           .collect::<Vec<Ciphertext>>()
                                                            ).collect();

        let mut b_iter = b.clone().into_iter();
        let b_mat: Vec<Vec<Scalar>> =  (0..self.m).map(|_| (0..self.n).map(|_| b_iter.next().unwrap())
                                                           .collect::<Vec<Scalar>>()
                                                            ).collect();
        let mut mexp_prover = MexpProver::new(C_mat, C_x, c_b.clone(), b_mat, s.clone(), rho_, self.com_ref.clone());
        let mexp_proof = MexpProof::default(); //mexp_prover.prove(trans, x.clone());

        //Challenge y, z
        let y = trans.challenge_scalar(b"y");
        let z = trans.challenge_scalar(b"z");

        self.chall.x = x.clone();
        self.chall.y = y.clone();
        self.chall.z = z.clone();

        let _z: Vec<Scalar> = (0..self.pi.len()).map(|_| -z.clone()).collect();
        let zeros: Vec<Scalar> = (0..self.m).map(|_| Scalar::zero()).collect();
        let c_z: Vec<RistrettoPoint> = self.com_ref.commit_vec(_z.clone(), zeros);

        let c_d: Vec<RistrettoPoint> = c_a.iter()
            .zip(c_b.iter())
            .map(|(a, b)| a* y + b)
            .collect();

        //Product Argument
        let d: Vec<Scalar> = a.iter()
            .zip(b.iter())
            .map(|(a_, b_)| a_* y + b_)
            .collect();
        let t: Vec<Scalar> = r.iter()
            .zip(s.iter())
            .map(|(r_, s_)| r_* y + s_)
            .collect();

        //Multiply cD * cZ
        let cd_cz: Vec<RistrettoPoint> = c_d.iter()
            .zip(c_z.iter())
            .map(|(d_, z_)| d_ + z_)
            .collect();//cd * c-z
        let d_z: Vec<Scalar> = d.iter()
            .zip(_z.iter())
            .map(|(d_, z_)| d_ + z_)
            .collect();
        let mut d_z_iter = d_z.clone().into_iter();
        let d_z: Vec<Vec<Scalar>> = (0..self.m).map(|_| 
                                                (0..self.n).map(|_| d_z_iter.next().unwrap()).collect::<Vec<Scalar>>()
                                                )
            .collect();
        //(yi + x^i − z)

        let product: Scalar = (1..=self.m * self.n).fold(
            Scalar::one(),
            |acc, i| 
            acc * (y.clone() * Scalar::from(i as u128) + x.pow(i as u64) - z.clone())
            );
        let mut prod_prover = ProdProver::new(cd_cz, 
                                              d_z, 
                                              t, 
                                              product, 
                                              self.com_ref.clone());
        prod_prover.chall = self.chall.clone();
        let prod_proof = prod_prover.prove(trans);



        ShuffleProof {
            c_A: c_a,
            c_B: c_b,
            mexp: mexp_proof,
            mexp_prover: mexp_prover,
            prod: prod_proof,
            prod_prover: prod_prover,
            y: y,
        }
    }

    pub fn verify(
        &mut self,
        trans: &mut Transcript,
        mut proof: ShuffleProof,
    ) -> Result<(), ProofError> {
        trans.shuffle_domain_sep(self.n as u64, self.m as u64);

        trans.val_append_point_vec(b"c_A", &proof.c_A.iter().map(|point| point.compress()).collect::<Vec<CompressedRistretto>>())?;

        trans.val_append_point_vec(b"c_B", &proof.c_B.iter().map(|point| point.compress()).collect::<Vec<CompressedRistretto>>())?;

        let mexp_proof = proof.mexp;
        let mut mexp_prover = proof.mexp_prover;
        //mexp_prover.verify(mexp_proof, trans)?;
        let mut prod_prover = proof.prod_prover;
        assert!(prod_prover.verify(trans, proof.prod).is_ok());

        Ok(())
    }
}
#[test]
fn test_prover_product() {
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
    let pi: Vec<u64> = cr.rand_perm(&(1..=((m*n) as u64)).collect());
    let rho: Vec<Scalar> = vec![cr.rand_scalar(); m*n];

    let base: Vec<Ciphertext> = cr.encrypt_vec(vec![Scalar::one(); m*n].as_slice(), &rho.clone());

    let C_pid: Vec<Ciphertext> = pi.iter()
        .map(|i|
             C_deck[i.clone() as usize - 1].clone()
             ).collect();

    let C_prime: Vec<Ciphertext> = base.into_iter()
        .zip(C_pid.into_iter())
        .map(|(b_, c_p)|
             b_ + c_p
             ).collect();

    //let pi: Vec<Scalar> = pi.into_iter().map(|p| Scalar::from(p)).collect();;
    let r: Vec<Scalar> = (0..m)
        .map(|_| cr.rand_scalar())
        .collect();
    let a: Vec<Scalar> = pi.iter()
        .map(|value| Scalar::from(value.clone() as u128))
        .collect();

    let c_a: Vec<RistrettoPoint> = cr.commit_vec(a.clone(), r.clone());
    let x = cr.rand_scalar();

    //Commit exp permutation
    let s: Vec<Scalar> = (0..m)
        .map(|_| cr.rand_scalar())
        .collect();

    let b: Vec<Scalar> = pi.iter()
        .map(|value| x.pow(value.clone()))
        .collect();

    let c_b: Vec<RistrettoPoint> = cr.commit_vec(b.clone(), s.clone());
    let y = cr.rand_scalar();
    let z = cr.rand_scalar();

    let _z: Vec<Scalar> = (0..pi.len()).map(|_| -z.clone()).collect();
    let zeros: Vec<Scalar> = (0..m).map(|_| Scalar::zero()).collect();
    let c_z: Vec<RistrettoPoint> = cr.commit_vec(_z.clone(), zeros);
    println!("{}", c_z.len());

    let c_d: Vec<RistrettoPoint> = c_a.iter()
        .zip(c_b.iter())
        .map(|(a, b)| a* y + b)
        .collect();

    //Product Argument
    let d: Vec<Scalar> = a.iter()
        .zip(b.iter())
        .map(|(a_, b_)| a_* y + b_)
        .collect();
    let t: Vec<Scalar> = r.iter()
        .zip(s.iter())
        .map(|(r_, s_)| r_* y + s_)
        .collect();

    //Multiply cD * cZ
    let cd_cz: Vec<RistrettoPoint> = c_d.iter()
        .zip(c_z.iter())
        .map(|(d_, z_)| d_ + z_)
        .collect();//cd * c-z
    let d_z: Vec<Scalar> = d.iter()
        .zip(_z.iter())
        .map(|(d_, z_)| d_ + z_)
        .collect();
    let mut d_z_iter = d_z.clone().into_iter();
    let d_z_mat: Vec<Vec<Scalar>> = (0..m).map(|_| 
                                            (0..n).map(|_| d_z_iter.next().unwrap()).collect::<Vec<Scalar>>()
                                            )
        .collect();
    //let d_z = d_z.to_col();
    //(yi + x^i − z)

    let product: Scalar = (1..=m * n).fold(
        Scalar::one(),
        |acc, i| 
        acc * (y.clone() * Scalar::from(i as u128) + x.pow(i as u64) - z.clone())
        );


    let mut prod_prover = ProdProver::new(cd_cz, 
                                          d_z_mat, 
                                          t, 
                                          product, 
                                          cr.clone());
    prod_prover.chall = Challenges{x:x, y:y, z:z};
    let mut trans = Transcript::new(b"testShuffleProof");
    let prod_proof = prod_prover.prove(&mut trans);
    let mut trans = Transcript::new(b"testShuffleProof");

    assert!(prod_prover.verify(&mut trans, prod_proof).is_ok());
}

#[test]
fn test_prover() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use crate::EGInp;
    
    let mut rng = StdRng::seed_from_u64(2);//from_entropy();
    let m: usize = 8;
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
    let pi: Vec<u64> = cr.rand_perm(&(1..=((m*n) as u64)).collect());
    let rho: Vec<Scalar> = vec![cr.rand_scalar(); m*n];

    let base: Vec<Ciphertext> = cr.encrypt_vec(vec![Scalar::one(); m*n].as_slice(), &rho.clone());

    let C_pid: Vec<Ciphertext> = pi.iter()
        .map(|i|
             C_deck[i.clone() as usize - 1].clone()
             ).collect();

    let C_prime: Vec<Ciphertext> = base.into_iter()
        .zip(C_pid.into_iter())
        .map(|(b_, c_p)|
             b_ + c_p
             ).collect();

    //let pi: Vec<Scalar> = pi.into_iter().map(|p| Scalar::from(p)).collect();;

    let mut prover_transcript = Transcript::new(b"testShuffleProof");
    let mut shuffle_prover = ShuffleProver::new(
                            m,
                            n,
                            m * n,
                            C_deck,
                            C_prime,
                            pi,
                            rho,
                            cr
                        );
                                
    let mut shuffle_proof = shuffle_prover.prove(&mut prover_transcript);

    let mut verifier_transcript = Transcript::new(b"testShuffleProof");

    assert!(shuffle_prover
            .verify(&mut verifier_transcript, shuffle_proof)
            .is_ok());
}

