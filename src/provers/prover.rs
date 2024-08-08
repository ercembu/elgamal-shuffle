#![allow(non_snake_case)]
use std::mem;

use rust_elgamal::{Scalar, Ciphertext};

use curve25519_dalek::ristretto::{RistrettoPoint, CompressedRistretto};
use merlin::Transcript;

use crate::arguers::CommonRef;

use crate::vec_utils::VecUtil::scalar_to_str;

use crate::provers::{mexp_prover::{MexpProof, MexpOptimProof, MexpProver},
                        sv_prover::{SVProver},
                        zero_prover::{ZeroProver},
                        hadamard_prover::{HadamProver},
                        prod_prover::{ProdProof, ProdProver}};

use crate::traits::{traits::{Timeable,
                                HeapSize,
                                EasySize,
                                Hadamard, 
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
pub struct ShuffleProof {
    pub(crate) c_A : Vec<RistrettoPoint>,
    pub(crate) c_B : Vec<RistrettoPoint>,
    pub(crate) mexp: MexpOptimProof,
    pub(crate) prod: ProdProof,
    pub(crate) y: Scalar,
}

impl HeapSize for ShuffleProof {
    fn heap_size(&self) -> usize {
        self.c_A.ez_size()
            + self.c_B.ez_size()
            + self.mexp.heap_size()
            + self.prod.heap_size()
            + self.y.ez_size()
    }
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

impl Timeable for ShuffleProver {

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
        println!("m: {}, n: {}, mu: {}", m, n, mu);
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

    pub fn prove(&mut self, trans: &mut Transcript) 
        -> (ZeroProver, SVProver, HadamProver, ProdProver, MexpProver, ShuffleProof)
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
        let mut mexp_prover = MexpProver::new(C_mat, 
                                              C_x, 
                                              c_b.clone(), 
                                              b_mat, 
                                              s.clone(), 
                                              rho_, 
                                              self.com_ref.clone());

        let mexp_time = mexp_prover.start_time();
        let mexp_proof = mexp_prover.prove_optim(trans, x.clone(), self.mu);

        println!("\n");
        println!("Opt Mexp Proof Time:\t{}", mexp_prover.elapsed(mexp_time));
        println!("Opt Mexp Proof Size:\t{}", mexp_proof.heap_size());
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
        let product_time = prod_prover.start_time();
        let (zero_prover, 
             sv_prover, 
             hadamard_prover, 
             prod_proof) = prod_prover.prove(trans);
        println!("\n");
        println!("Product Proof Time:\t{}", prod_prover.elapsed(product_time));
        println!("Product Proof Size:\t{}", mem::size_of_val(&prod_proof));



        (zero_prover,
            sv_prover,
            hadamard_prover,
            prod_prover,
            mexp_prover,
            ShuffleProof {
                c_A: c_a,
                c_B: c_b,
                mexp: mexp_proof,
                prod: prod_proof,
                y: y,
            }
        )
    }

    pub fn verify(
        &mut self,
        trans: &mut Transcript,
        mut proof: ShuffleProof,
        mut zero_prover: ZeroProver,
        mut sv_prover: SVProver,
        mut hadamard_prover: HadamProver,
        mut prod_prover: ProdProver,
        mut mexp_prover: MexpProver,
    ) -> Result<(), ProofError> {
        trans.shuffle_domain_sep(self.n as u64, self.m as u64);

        let mexp_proof = proof.mexp;
        let mexp_time = mexp_prover.start_time();
        mexp_prover.verify_optim(mexp_proof, trans, self.mu)?;
        println!("\n");
        println!("Opt Mexp Verify Time:\t{}", mexp_prover.elapsed(mexp_time));
        let prod_time = prod_prover.start_time();
        assert!(prod_prover.verify(trans, proof.prod, 
                                   hadamard_prover,
                                   zero_prover,
                                   sv_prover).is_ok());
        println!("\n");
        println!("Product Verify Time:\t{}", prod_prover.elapsed(prod_time));

        Ok(())
    }
}
#[test]
fn test_product_prover() {
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
    let (zero_prover, 
         sv_prover, 
         hadam_prover,
         prod_proof) = prod_prover.prove(&mut trans);
    let mut trans = Transcript::new(b"testShuffleProof");

    assert!(prod_prover.verify(&mut trans, prod_proof,
                               hadam_prover,
                               zero_prover,
                               sv_prover).is_ok());
}

#[test]
fn test_prover_obs() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use crate::EGInp;
    use std::time::SystemTime;
    
    let mut rng = StdRng::seed_from_u64(2);//from_entropy();
    let m: usize = 8;
    let n: usize = 8;
    
    let mu: usize = 1; 

    let setup_time = SystemTime::now();
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
    println!("Setup Time:\t{}", setup_time.elapsed().unwrap().as_millis());

    let perm_time = SystemTime::now();
    let pi: Vec<u64> = cr.rand_perm(&(1..=((m*n) as u64)).collect());
    let rho: Vec<Scalar> = vec![cr.rand_scalar(); m*n];

    let base: Vec<Ciphertext> = cr.encrypt_vec(vec![Scalar::zero(); m*n].as_slice(), &rho.clone());

    let C_pid: Vec<Ciphertext> = pi.iter()
        .map(|i|
             C_deck[i.clone() as usize - 1].clone()
             ).collect();

    let C_prime: Vec<Ciphertext> = base.into_iter()
        .zip(C_pid.into_iter())
        .map(|(b_, c_p)|
             b_ + c_p
             ).collect();

    println!("Perm Time:\t{}", perm_time.elapsed().unwrap().as_millis());

    let mut prover_transcript = Transcript::new(b"testShuffleProof");
    let mut shuffle_prover = ShuffleProver::new(
                            m,
                            n,
                            mu,
                            C_deck,
                            C_prime,
                            pi,
                            rho,
                            cr
                        );
    let proof_time = shuffle_prover.start_time();
                                
    let (mut zero_prover, 
             mut sv_prover, 
             mut hadam_prover, 
             mut prod_prover, 
             mut mexp_prover, 
             mut shuffle_proof) = shuffle_prover.prove(&mut prover_transcript);
    

    println!("\n");
    println!("Shuffle Proof Time:\t{}", shuffle_prover.elapsed(proof_time));
    println!("Shuffle Proof Size:\t{}", shuffle_proof.heap_size());

    let mut verifier_transcript = Transcript::new(b"testShuffleProof");

    let verify_time = shuffle_prover.start_time();
    assert!(shuffle_prover
            .verify(&mut verifier_transcript, shuffle_proof,
                    zero_prover,
                    sv_prover,
                    hadam_prover,
                    prod_prover,
                    mexp_prover)
            .is_ok());
    println!("\n");
    println!("Shuffle Verify Time:\t{}", shuffle_prover.elapsed(verify_time));
}

