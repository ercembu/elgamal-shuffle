use std::ops::Index;
use len_trait::len::*;

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
