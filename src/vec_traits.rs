use std::vec::Vec;
use std::iter;
use std::borrow::Borrow;
use std::fmt;

use rust_elgamal::{Scalar};

use crate::traits::Hadamard;

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

