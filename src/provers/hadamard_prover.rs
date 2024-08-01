#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext};
use std::iter;

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::traits::MultiscalarMul;
use merlin::Transcript;

use crate::arguers::CommonRef;

use crate::traits::{traits::{Hadamard, 
                                EGMult, 
                                InnerProduct, 
                                Multiplicat,
                                Addition
                            }, 
                    mat_traits::MatTraits};

use crate::utils::{utils::Challenges,
                    transcript::TranscriptProtocol,
                    errors::ProofError};

use crate::provers::zero_prover::{ZeroProof, ZeroProver};

///Data struct for final proof arguments 
///to be sent to the verification method
#[derive(Clone)]
pub struct HadamProof{
    ///Commitments to the result : B<sub>i</sub> 
    ///= Hadamard (A<sub>0</sub>..A<sub>i</sub>)
    c_Bi: Vec<RistrettoPoint>,
    ///Blinded commitments of c<sub>Bi</sub> 
    c_D: RistrettoPoint,
    ///Commited `-1` vector value
    c_1: RistrettoPoint,
    ///Zero Argument Proof
    zero_proof: ZeroProof,
    ///Zero Argument Prover to be used in verification
    zero_prover: ZeroProver,
}

///Struct for initial Hadamard Proof Arguments
///
///Main objective is to show HadamardSum(A) == b
#[derive(Clone)]
pub struct HadamProver {
    ///Commitments to com<sub>ck</sub>(A, r)
    c_A: Vec<RistrettoPoint>,
    ///Commitments to com<sub>ck</sub>(b, s)
    c_b: RistrettoPoint,
    ///Open value matrices, and their blinding factors
    ///A: Z<sup>nxm</sup>
    A: Vec<Vec<Scalar>>,
    ///r: Z<sup>m</sup>
    r: Vec<Scalar>,
    ///B: Z<sup>nxm</sup>
    B: Vec<Vec<Scalar>>,
    ///s: Z<sup>m</sup>
    s: Vec<Scalar>,
    com_ref: CommonRef,
    /// Challenges from oracle, purely random
    pub(crate) chall: Challenges
}

impl HadamProver {
    ///Base Contstructor
    pub fn new(
        c_A: Vec<RistrettoPoint>,
        c_b: RistrettoPoint,
        A: Vec<Vec<Scalar>>,
        r: Vec<Scalar>,
        b: Vec<Scalar>,
        s: Scalar,
        com_ref: CommonRef
    ) -> Self {
        Self {
            c_A: c_A,
            c_b: c_b,
            A: A,
            s: s,
            b: b,
            r: r,
            com_ref: com_ref,
            chall: Challenges{x:Scalar::one(), y:Scalar::one(), z:Scalar::one()},
        }
    }

    ///Bilinear map to be used in the Zero Argument
    ///
    ///Essentially a blinded hadamard sum
    pub fn exp_dot(
        a: Vec<Scalar>,
        d: Vec<Scalar>,
        y: Scalar
    ) -> Scalar {
        assert!(a.len() == d.len());
        (0..a.len())
            .fold(
                Scalar::zero(),
                |acc, j|
                acc + (a[j] * d[j] * y.pow((j+1) as u64))
                )
    }

    ///prove method that creates a HadamardProof
    ///
    ///Main method used works over the mathmathical expression
    ///
    ///0 = Sum<sub>i=1</sub><sup>m-1</sup>(A<sub>i+1</sub> * 
    ///d<sub>i</sub> - `1` * d)
    ///
    ///where d and d<sub>i</sub> suggest 
    ///the last and rest of the openings to c<sub>D</sub> 
    ///
    ///reduction to a zero argument via 
    ///blinded products of hadamard sums
    pub fn prove(
        &mut self,
        trans: &mut Transcript,
    ) -> HadamProof {
        trans.append_message(b"dom_sep", b"HadamardProof");

        let m = self.A.len();
        let n = self.A[0].len();


        let B: Vec<Vec<Scalar>> = (0..m).map(
            |i| {
                match i {
                    i if i == m-1 => self.b.clone(),
                    _ => (0..i+1).fold(
                            vec![Scalar::one(); n],
                            |acc, j|
                            acc.hadamard(self.A[j].clone())
                        ),
                    }
            }
            ).collect();

        let s_vec: Vec<Scalar> = (0..m-2).map(|_| self.com_ref.rand_scalar())
            .collect();

        let s_vec: Vec<Scalar> = iter::once(self.r[0]).chain(s_vec.into_iter())
            .chain(iter::once(self.s)).collect();

        let c_B: Vec<RistrettoPoint> = (0..B.len()).map(
            |i|
            match i {
                0 => self.c_A[0].clone(),
                i if i == m-1 => self.c_b.clone(),
                _ => self.com_ref.commit(B[i].clone(), s_vec[i].clone()),
            }
        ).collect();

        let x = self.chall.x.clone();
        let y = self.chall.y.clone();

        let x_pow: Vec<Scalar> = (0..c_B.len()).map(|i| x.pow((i+1) as u64)).collect();
        let c_Di: Vec<RistrettoPoint> = (0..c_B.len()-1).map(
            |i|
            c_B[i] * x_pow[i].clone()
        ).collect();

        let c_D: RistrettoPoint = RistrettoPoint::multiscalar_mul(&x_pow[0..m-1],
                                                                  &c_B[1..m]);

        let c_1: RistrettoPoint = self.com_ref.commit(vec![-Scalar::one(); n], 
                                      Scalar::zero());

        let D: Vec<Vec<Scalar>> = (0..m-1).map(
            |i|
            B[i].mult(&x.pow((i+1) as u64))
        ).collect();

        let t_: Vec<Scalar> = (0..m-1).map(
            |i|
            s_vec[i] * x.pow((i+1) as u64)
        ).collect();

        let d: Vec<Scalar> = (0..m-1).fold(
            vec![Scalar::zero(); B[0].len()],
            |acc, i|
            acc.add(&B[i+1].mult(&x.pow((i+1) as u64)))
            );

        let t: Scalar = (0..m-1).fold(
            Scalar::zero(),
            |acc, i|
            acc + x.pow((i+1) as u64) * s_vec[i+1]
            );

        let mut zero_prover: ZeroProver = ZeroProver::new(
            [&self.c_A[1..m], &[c_1].as_slice()].concat(),
            [&c_Di[0..m-1], &[c_D].as_slice()].concat(),
            HadamProver::exp_dot,
            y.clone(),
            [&self.A.clone()[1..m], &[vec![-Scalar::one(); n]]].concat(),
            [&self.r.clone()[1..m], &[Scalar::zero()]].concat(),
            [D.clone().as_slice(), &[d]].concat(),
            [t_.clone().as_slice(), &[t]].concat(),
            self.com_ref.clone());
        zero_prover.chall = self.chall.clone();

        let zero_proof: ZeroProof = zero_prover.prove(trans);
        
        HadamProof {
            c_Bi: c_B,
            c_D: c_D,
            c_1: c_1,
            zero_proof: zero_proof,
            zero_prover: zero_prover,
        }
    }
    
    ///verify method that verifies a HadamardProof
    pub fn verify(
        &mut self,
        trans: &mut Transcript,
        proof: HadamProof,
    ) -> Result<(), ProofError> {
        trans.append_message(b"dom_sep", b"HadamardProof");
        let m = proof.c_Bi.len();
        let n = self.A[0].len();
        assert!(proof.c_Bi[0] == self.c_A[0]);
        assert!(proof.c_Bi[m-1] == self.c_b);
        assert!(proof.c_1 == self.com_ref.commit(vec![-Scalar::one();n],
                                                    Scalar::zero()));
        let mut zero_prover = proof.zero_prover;
        zero_prover.verify(trans, proof.zero_proof)?;
        Ok(())
    }
}

#[test]
fn test_base() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    let mut prover_transcript = Transcript::new(b"testHadamProof");

    let mut rng = StdRng::seed_from_u64(2);//from_entropy();
    let m: usize = 6;
    let n: usize = 4;
    let mut com_ref = CommonRef::new((m*n) as u64, rng);
    let a: Vec<Vec<Scalar>> = vec![vec![com_ref.rand_scalar(); m]; n];
    let r: Vec<Scalar> = vec![com_ref.rand_scalar(); m];
    let s: Scalar = com_ref.rand_scalar();

    let b: Vec<Scalar> = a.to_col().iter()
        .fold(
            vec![Scalar::one(); n],
            |acc, a_i|
            acc.hadamard(a_i.clone())
            );

    let c_A: Vec<RistrettoPoint> = com_ref.commit_mat(a.clone(), r.clone());
    let c_B: RistrettoPoint = com_ref.commit(b.clone(), s.clone());

    let mut hadam_prover =  HadamProver::new(
        c_A,
        c_B,
        a.to_col(),
        r,
        b,
        s,
        com_ref.clone()
    );

    let mut hadam_proof = hadam_prover.prove(&mut prover_transcript);

    let mut verifier_transcript = Transcript::new(b"testHadamProof");
    assert!(hadam_prover.verify(&mut verifier_transcript, hadam_proof).is_ok());


}
