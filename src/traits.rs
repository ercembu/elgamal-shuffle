use std::vec::Vec;
use std::iter;
use std::borrow::Borrow;

use rust_elgamal::{Scalar, Ciphertext};
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

pub trait EGMult<I> {
    type Output;
    fn pow(&self, exp: I) -> Self::Output;
}

impl EGMult<usize> for Scalar {
    type Output = Scalar;
    fn pow(&self, exp: usize) -> Scalar {
        let mut result = Scalar::one();
        for _ in 0..exp {
            result *= self;
        }
        result
    }
}

impl EGMult<Scalar> for Ciphertext {
    type Output = Ciphertext;

    fn pow(&self, exp: Scalar) -> Ciphertext {
        self * exp
    }
}

impl EGMult<&[Scalar]> for &[Ciphertext] {
    type Output = Ciphertext;

    fn pow(&self, exp: &[Scalar]) -> Ciphertext {
        assert!(self.len() == exp.len(), "arguments missized");
        self.iter()
            .zip(exp)
            .map(|(c, a)| c * a)
            .collect::<Vec<Ciphertext>>()
            .into_iter()
            .reduce(|s1, s2| s1 + s2)
            .unwrap()
    }
}

impl EGMult<&[Vec<Scalar>]> for &[Ciphertext] {
    type Output = Vec<Ciphertext>;

    fn pow(&self, exp: &[Vec<Scalar>]) -> Vec<Ciphertext> {
        assert!(self.len() == exp.len());
        exp.iter().map(|a| self.pow(&a[..])).collect()
    }
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
