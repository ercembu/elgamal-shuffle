#![allow(non_snake_case)]
#![allow(dead_code)]
use rust_elgamal::{EncryptionKey, 
                    DecryptionKey, 
                    Scalar, 
                    Ciphertext};

use rand::rngs::StdRng;
use rand::prelude::SliceRandom;

use std::vec::Vec;
use std::iter;

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::traits::MultiscalarMul;

use crate::enums::EGInp;

#[derive(Clone)]
pub struct CommonRef {
    pub pk : EncryptionKey,
    pub dk : DecryptionKey,
    pub ck : Vec<RistrettoPoint>,
    pub rng: StdRng,
}

impl CommonRef {
    pub fn new(n: u64, mut rng: StdRng) -> Self {
        let dk = DecryptionKey::new(&mut rng);
        let pk = dk.encryption_key();

        let pd_gen: Vec<RistrettoPoint> = vec![RistrettoPoint::random(&mut rng); (n as usize) + 1];
        Self{pk: *pk, ck: pd_gen, dk: dk, rng: rng}
    }

    //Correct form: |a| = N; N = m * n; |r| = m;
    pub fn commit_vec(
        &mut self,
        a: Vec<Scalar>,
        r: Vec<Scalar>,
    ) -> Vec<RistrettoPoint> {
        assert!(a.len() % r.len() == 0);
        let n = a.len() / r.len();
        let mut A_iter = a.into_iter();
        r.into_iter()
            .map(|rand| {
                let a: Vec<Scalar> = (0..n).map(|_| A_iter.next().unwrap())
                    .collect();
                self.commit(a, rand)
            })
            .collect()
    }
    pub fn commit(
        &mut self,
        a: Vec<Scalar>,
        r: Scalar,
    ) -> RistrettoPoint {
        let (Gs, H) = self.ck.split_at(self.ck.len() - 1);
        let H = H[0];

        let a_len = a.len();
        let padded: Vec<Scalar> = a.into_iter()
                                    .chain(vec![Scalar::zero(); Gs.len() - a_len]
                                           .into_iter()
                                        ).collect();

        assert!(padded.len() == Gs.len());

        RistrettoPoint::multiscalar_mul(iter::once(r)
                                            .chain(padded.into_iter())
                                            .collect::<Vec<Scalar>>(), 
                                        iter::once(H)
                                        .chain(Gs.into_iter()
                                               .map(|G| *G))
                                        .collect::<Vec<RistrettoPoint>>())

        
        
        //self.ck.commit(*msg, Scalar::random(&mut self.rng))
    }
    
    pub fn commit_mat(
        &mut self,
        A: Vec<Vec<Scalar>>,
        r: Vec<Scalar>,
    ) -> Vec<RistrettoPoint> {
        assert!(A[0].len() == r.len());
        let n = A.len();
        let m = A[0].len();
        r.into_iter()
            .zip(0..m)
            .map( |(rand, a)| {
                let col: Vec<Scalar> = (0..n).map(|i| A[i][a]).collect();
                self.commit(col, rand)
            })
            .collect()

    }

    pub fn encrypt(
        &mut self,
        m: &EGInp,
        r: &Scalar,
    ) -> Ciphertext {
        match m {
            EGInp::Rist(x) => self.pk.encrypt_with(x.clone(), r.clone()),
            EGInp::Scal(x) => self.pk.exp_encrypt_with(x.clone(), r.clone()),
        }
    }

    pub fn encrypt_vec(
        &mut self,
        M: &[Scalar],
        R: &[Scalar],
    ) -> Ciphertext {
        let m_prod = M.iter().copied().reduce(|s1, s2| (s1 * s2)).unwrap();
        let r_prod = R.iter().copied().reduce(|s1, s2| (s1 * s2)).unwrap();
        self.pk.exp_encrypt_with(m_prod, r_prod) 
    }

    pub fn decrypt(
        &self,
        ctext: Ciphertext
    ) -> RistrettoPoint {
        self.dk.decrypt(ctext)
    }

    pub fn rand_scalar(
        &mut self
    ) -> Scalar {
        Scalar::random(&mut self.rng)
    }

    pub fn rand_perm(
        &mut self,
        range: &Vec<u64>,
    ) -> Vec<u64> {
        let mut res: Vec<u64> = (*range.clone()).to_vec();
        res.shuffle(&mut self.rng);
        res
    }
}

#[test]
fn test_common() {

}
