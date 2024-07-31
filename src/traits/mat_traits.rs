use std::ops::Index;
use len_trait::len::*;

use rust_elgamal::{Scalar};
use crate::traits::traits::{Addition, Multiplicat};

pub trait MatTraits<I, Idx=usize>
where
    Self: IntoIterator,
    Self: Index<Idx>,
    Self: Len,
    Self::Output: Len,
    Self: Len,
    Self::Item: Len,
{

    fn size(&self) -> (usize, usize);
    fn to_col(&self) -> Vec<Vec<I>>;
}

impl<I> MatTraits<I, usize> for Vec<Vec<I>> 
where I: Clone,
      I: Copy
{
    fn size(&self) -> (usize, usize) {
        (self.len(), self.index(0).len())
    }


    fn to_col(&self) -> Vec<Vec<I>> 
    {
        let n = self.len();
        let m = self.index(0).len();

        let result: Vec<Vec<I>> = (0..m).map(|i| (0..n)
            .map(|j| self.index(j).index(i).clone())
                .collect()
            ).collect();

        result

    }


}
impl Multiplicat for Vec<Vec<Scalar>>
{
    type RHS = Vec<Scalar>;
    type Out = Vec<Scalar>;

    fn mult(&self, rhs: &Vec<Scalar>) -> Vec<Scalar> {
        assert!(self.len() == rhs.len());
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
        assert!(self.len() == rhs.len());
        self.iter()
            .zip(rhs.iter())
            .map(|(row, x)| row.mult(x))
            .reduce(|acc, x| acc.add(&x))
            .unwrap()
    }
}
