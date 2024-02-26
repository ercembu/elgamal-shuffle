use rust_elgamal::{EncryptionKey, 
                    DecryptionKey, 
                    Scalar, 
                    GENERATOR_TABLE, 
                    Ciphertext};

use rand::rngs::StdRng;
use rand::SeedableRng;

use curve25519_dalek::ristretto::RistrettoPoint;
use bulletproofs::PedersenGens;

#[derive(Clone)]
pub struct CommonRef {
    pub pk : EncryptionKey,
    pub dk : DecryptionKey,
    pub ck : PedersenGens,
    pub rng: StdRng,
}

impl CommonRef {
    pub fn new(mut rng: StdRng) -> Self {
        let dk = DecryptionKey::new(&mut rng);
        let pk = dk.encryption_key();

        let pd_gen: PedersenGens = PedersenGens::default();
        Self{pk: *pk, ck: pd_gen, dk: dk, rng: rng}
    }

    pub fn commit(
        &mut self,
        msg: &Scalar,
    ) -> RistrettoPoint {
        self.ck.commit(*msg, Scalar::random(&mut self.rng))
    }

    pub fn encrypt(
        &mut self,
        msg: &RistrettoPoint,
    ) -> Ciphertext {
        self.pk.encrypt(*msg, &mut self.rng)
    }

    pub fn decrypt(
        &self,
        ctext: Ciphertext
    ) -> RistrettoPoint {
        self.dk.decrypt(ctext)
    }
}
