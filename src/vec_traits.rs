use std::vec::Vec;
use std::borrow::Borrow;

use rust_elgamal::{Scalar};

use crate::traits::{Hadamard, InnerProduct};

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
