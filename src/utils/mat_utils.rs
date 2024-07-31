#![allow(dead_code)]
use rust_elgamal::{Scalar};
use crate::utils::vec_utils::VecUtil;

pub mod MatUtil {
    use super::*;
    pub fn print_scalar_mat(m: &Vec<Vec<Scalar>>) -> String {
        let mut result: String = String::from("[");
        for v in m {
            result.push_str(VecUtil::scalar_to_str(&v).as_str());
            result.push_str(",\n");
        }
        result.push_str("]");

        result

    }
}
