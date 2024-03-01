use ethnum::I256;
use rust_elgamal::{Scalar};

pub struct VecUtil;
impl VecUtil {
    pub fn scalar_to_str(v: &Vec<Scalar>) -> String {
        let mut result: String = String::from("[");
        for scalar in v {
            let str_result: String;
            let sc_str = I256::from_le_bytes(*scalar.as_bytes());
            if sc_str.to_string().len() > 10 {
                let m_one = I256::from_le_bytes((-Scalar::one().reduce()).to_bytes());
                str_result = (sc_str - (m_one + 1)).to_string();
            } else {str_result = sc_str.to_string();}
            result += &str_result;
            result.push_str(", ");
        }
        result.push_str("]");

        result

    }
}
