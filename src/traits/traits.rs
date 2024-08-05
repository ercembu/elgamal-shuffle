use std::vec::Vec;
use std::borrow::Borrow;

use crate::utils::vec_utils::VecUtil;

use rust_elgamal::{Scalar, Ciphertext};

use std::time::SystemTime;

pub trait Timeable {
    
    fn start_time(&self) -> SystemTime {
        SystemTime::now()
    }

    fn elapsed(&self, time: SystemTime) -> u128 {
        time.elapsed().unwrap().as_millis()
    }
}

pub trait Hadamard {
    type Msg;

    fn hadamard<I>(&self, rhs: I) -> Vec<Self::Msg>
    where 
        I: IntoIterator,
        I::Item: Borrow<Self::Msg>;

}

pub trait InnerProduct<I> {
    type Output;

    fn dot(&self, rhs: &I) -> Self::Output
    where 
        I: IntoIterator,
        I::Item: Borrow<Self::Output>;
}

pub trait Multiplicat {
    type RHS;
    type Out;

    fn mult(&self, rhs: &Self::RHS) -> Self::Out;

}
pub trait Addition {
    type Output;
    fn add(&self, rhs: &Self::Output) -> Self::Output;
}


pub trait EGMult<I> {
    type Output;
    fn pow(&self, exp: I) -> Self::Output;
}

impl EGMult<u64> for Scalar {
    type Output = Scalar;
    fn pow(&self, exp: u64) -> Scalar {
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
impl EGMult<&[Scalar]> for &[&Ciphertext] {
    type Output = Ciphertext;

    fn pow(&self, exp: &[Scalar]) -> Ciphertext {
        assert!(self.len() == exp.len(), "arguments missized");
        self.iter()
            .zip(exp)
            .map(|(c, a)| *c * a)
            .collect::<Vec<Ciphertext>>()
            .into_iter()
            .reduce(|s1, s2| s1 + s2)
            .unwrap()
    }
}
impl EGMult<&[&Scalar]> for &[&Ciphertext] {
    type Output = Ciphertext;

    fn pow(&self, exp: &[&Scalar]) -> Ciphertext {
        assert!(self.len() == exp.len(), "arguments missized");
        self.iter()
            .zip(exp)
            .map(|(c, a)| *c * *a)
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

}
