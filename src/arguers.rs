#![allow(non_snake_case)]
use rust_elgamal::{EncryptionKey, 
                    DecryptionKey, 
                    Scalar, 
                    GENERATOR_TABLE, 
                    Ciphertext};

use rand::rngs::StdRng;
use rand::SeedableRng;

use std::vec::Vec;
use std::iter;

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::traits::MultiscalarMul;
use bulletproofs::PedersenGens;

#[derive(Clone)]
pub struct CommonRef {
    pub pk : EncryptionKey,
    pub dk : DecryptionKey,
    pub ck : Vec<RistrettoPoint>,
    pub rng: StdRng,
}

impl CommonRef {
    pub fn new(n: usize, mut rng: StdRng) -> Self {
        let dk = DecryptionKey::new(&mut rng);
        let pk = dk.encryption_key();

        let pd_gen: Vec<RistrettoPoint> = vec![RistrettoPoint::random(&mut rng); n+1];
        Self{pk: *pk, ck: pd_gen, dk: dk, rng: rng}
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

    pub fn encrypt(
        &mut self,
        m: &RistrettoPoint,
        r: &Scalar,
    ) -> Ciphertext {
        self.pk.encrypt_with(*m, r.clone())
    }

    pub fn decrypt(
        &self,
        ctext: Ciphertext
    ) -> RistrettoPoint {
        self.dk.decrypt(ctext)
    }
}

#[test]
fn test_common() {
    let m = 4;
    let n = 4;
    let mu = 2;

    let N = mu * n;
    let mut rng = StdRng::from_entropy();
    let mut cr = CommonRef::new(n, rng);

    let mut rng = StdRng::from_entropy();
    let a = cr.commit(vec![Scalar::zero(); n], Scalar::random(&mut rng));


    let message: RistrettoPoint =  &Scalar::from(5u32) * &GENERATOR_TABLE;
    let encrypted = cr.encrypt(&message, &Scalar::random(&mut rng));
    let decrypted = cr.decrypt(encrypted);
    assert_eq!(message, decrypted);
    println!("done");

}
