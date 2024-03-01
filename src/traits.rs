use std::vec::Vec;
use std::iter;
use std::borrow::Borrow;

use rust_elgamal::{Scalar};
use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::traits::MultiscalarMul;

use crate::vec_utils::VecUtil;

pub trait Hadamard {
    type Msg;

    fn hadamard<I>(&self, rhs: I) -> Vec<Self::Msg>
    where 
        I: IntoIterator,
        I::Item: Borrow<Self::Msg>;

}


#[test]
fn test_hadamard() {
    let n = 4;
    let a: Vec<Scalar> = vec![Scalar::from(2 as u64); n];
    let b: Vec<Scalar> = vec![Scalar::from(3 as u64); n];

    let res: Vec<Scalar> = a.hadamard(&b);
    assert!(res.len() == a.len());

    println!("{}", VecUtil::scalar_to_str(&res));
}
