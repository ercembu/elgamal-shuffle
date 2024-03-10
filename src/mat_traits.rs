use std::ops::Index;
use len_trait::len::*;

use rust_elgamal::{Scalar};
use crate::traits::{Addition, Multiplicat};

pub trait MatTraits<Idx=usize>
where
    Self: IntoIterator,
    Self: Index<Idx>,
    Self::Output: Len,
    Self: Len,
    Self::Item: Len,
{

    fn size(&self) -> (usize, usize);
}

impl<I> MatTraits<usize> for Vec<Vec<I>> {
    fn size(&self) -> (usize, usize) {
        (self.len(), self.index(0).len())
    }


}
impl Multiplicat for Vec<Vec<Scalar>>
{
    type RHS = Vec<Scalar>;
    type Out = Vec<Scalar>;

    fn mult(&self, rhs: &Vec<Scalar>) -> Vec<Scalar> {
        assert!(self.index(0).len() == rhs.len());
        self.iter()
            .zip(rhs.iter())
            .map(|(row, x)| row.mult(x))
            .reduce(|acc, x| acc.add(&x))
            .unwrap()
    }
}
impl Multiplicat for Vec<Vec<&Scalar>>
{
    type RHS = Vec<Scalar>;
    type Out = Vec<Scalar>;

    fn mult(&self, rhs: &Vec<Scalar>) -> Vec<Scalar> {
        assert!(self.index(0).len() == rhs.len());
        self.iter()
            .zip(rhs.iter())
            .map(|(row, x)| row.mult(x))
            .reduce(|acc, x| acc.add(&x))
            .unwrap()
    }
}
