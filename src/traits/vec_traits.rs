use std::vec::Vec;
use std::borrow::Borrow;

use rust_elgamal::{Scalar};

use crate::traits::traits::{Hadamard, InnerProduct, Addition, Multiplicat};

impl Hadamard for Vec<Scalar> {
    type Msg = Scalar;

    fn hadamard<I>(&self, rhs: I) -> Vec<Scalar>
    where
        I: IntoIterator, 
        I::Item: Borrow<Scalar>,
    {
        self.iter().zip(rhs).map(|(a, b)| a * b.borrow()).collect()
    }
}

impl InnerProduct<Vec<Scalar>> for Vec<Scalar> {
    type Output = Scalar;

    fn dot(&self, rhs: &Vec<Scalar>) -> Scalar 
    {
        assert!(self.len() == rhs.len());

        
        self.iter()
            .zip(rhs.iter())
            .fold(Scalar::zero(), |acc, (l, r)| acc + (l * r))
    }
}

impl Addition for Vec<Scalar> {
    type Output = Vec<Scalar>;

    fn add(&self, rhs: &Vec<Scalar>) -> Vec<Scalar> {
        self.iter().zip(rhs.iter())
            .map(|(x, y)| x + y)
            .collect::<Vec<Scalar>>()
    }
}
impl Multiplicat for Vec<Scalar> {
    type RHS = Scalar;
    type Out = Vec<Scalar>;

    fn mult(&self, rhs: &Scalar) -> Vec<Scalar> {
        self.iter()
            .map(|x| x * rhs)
            .collect::<Vec<Scalar>>()
    }
}
impl Multiplicat for &Vec<Scalar> {
    type RHS = Scalar;
    type Out = Vec<Scalar>;

    fn mult(&self, rhs: &Scalar) -> Vec<Scalar> {
        self.iter()
            .map(|x| x * rhs)
            .collect::<Vec<Scalar>>()
    }
}
impl Multiplicat for Vec<&Scalar> {
    type RHS = Scalar;
    type Out = Vec<Scalar>;

    fn mult(&self, rhs: &Scalar) -> Vec<Scalar> {
        self.iter()
            .map(|x| *x * rhs)
            .collect::<Vec<Scalar>>()
    }
}
impl Multiplicat for &Vec<&Scalar> {
    type RHS = Scalar;
    type Out = Vec<Scalar>;

    fn mult(&self, rhs: &Scalar) -> Vec<Scalar> {
        self.iter()
            .map(|x| *x * rhs)
            .collect::<Vec<Scalar>>()
    }
}
