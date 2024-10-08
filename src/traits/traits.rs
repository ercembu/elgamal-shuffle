use std::mem;
use std::vec::Vec;
use std::borrow::Borrow;
use std::iter::IntoIterator;

use crate::utils::vec_utils::VecUtil;

use rust_elgamal::{Scalar, Ciphertext, RistrettoPoint};

use std::time::SystemTime;

pub trait Timeable {
    
    fn start_time(&self) -> SystemTime {
        SystemTime::now()
    }

    fn elapsed(&self, time: SystemTime) -> u128 {
        time.elapsed().unwrap().as_millis()
    }
}

pub trait EasySize {
    fn ez_size(&self) -> usize{
        mem::size_of_val(self)
    }
}

impl<T> EasySize for Vec<T> {
    fn ez_size(&self) -> usize {
        self.capacity() * mem::size_of::<T>()
    }
}

impl EasySize for Scalar {
}

impl EasySize for RistrettoPoint {
    fn ez_size(&self) -> usize{
        mem::size_of_val(&self.compress())
    }
}

impl EasySize for Ciphertext{
    fn ez_size(&self) -> usize{
        let inner = self.inner();
        mem::size_of_val(&inner.0.compress()) + mem::size_of_val(&inner.1.compress())
    }
}

pub trait HeapSize {
    fn heap_size(&self) -> usize;
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
