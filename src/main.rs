use rust_elgamal::{DecryptionKey, Scalar, GENERATOR_TABLE};
use rand::rngs::StdRng;
use rand::SeedableRng;

use curve25519_dalek::ristretto::RistrettoPoint;
use bulletproofs::PedersenGens;

mod arguers;


fn main() {
    let mut rng = StdRng::from_entropy();

    let mut cr = arguers::CommonRef::new(rng);


    let message: RistrettoPoint =  &Scalar::from(5u32) * &GENERATOR_TABLE;
    let encrypted = cr.encrypt(&message);
    let decrypted = cr.decrypt(encrypted);
    assert_eq!(message, decrypted);
    println!("done");
}
