#![allow(non_snake_case)]
use rust_elgamal::{Scalar, Ciphertext, IsIdentity};
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

///Data struct for final proof arguments 
///to be sent to the verification method
#[derive(Clone)]
pub struct ZeroProof {
    ///Commitment to blinding a<sub>0</sub>: com<sub>ck</sub>(a<sub>0</sub>, r)
    c_A0: RistrettoPoint,
    ///Commitment to blinding b<sub>m</sub>: com<sub>ck</sub>(b<sub>m</sub>, s)
    c_Bm: RistrettoPoint,
    ///Commitment to mapped blinded column vectors
    ///
    ///com<sub>ck</sub>(A<sub>blind</sub> * B<sub>blind</sub>, t) 
    c_D: Vec<RistrettoPoint>,

    ///Openings to the commitment values
    a_vec: Vec<Scalar>,
    b_vec: Vec<Scalar>,
    r: Scalar,
    s: Scalar,
    t: Scalar,
    x: Scalar,
}

///Struct for initial Zero Proof Arguments
///
///Main goal of the arguments is to show bi_map(A, B) == 0
#[derive(Clone)]
pub struct ZeroProver {
    ///Commitments to A column vectors
    c_Ai: Vec<RistrettoPoint>,
    
    ///Commitments to B column vectors
    c_Bi: Vec<RistrettoPoint>,

    ///Bilinear mapping fuction Z^n X Z^n -> Z
    ///
    ///In our case it is a weighted dot sum for column vectors
    ///injected from HadamardProver
    bi_map: fn(Vec<Scalar>, Vec<Scalar>, Scalar)->Scalar,

    ///Challenge scalar y to be used in the bi_map
    y: Scalar,

    ///Open value matrices, and their hiding factors
    A: Vec<Vec<Scalar>>,
    r: Vec<Scalar>,
    B: Vec<Vec<Scalar>>,
    s: Vec<Scalar>,

    ///Common reference data
    com_ref: CommonRef,

    /// Challenges from oracle, purely random
    pub(crate) chall: Challenges
}

impl ZeroProver {
    ///Implements base constructor
    pub fn new(
        c_Ai: Vec<RistrettoPoint>,
        c_Bi: Vec<RistrettoPoint>,
        bi_map: fn(Vec<Scalar>, Vec<Scalar>, Scalar)->Scalar,
        y: Scalar,
        A: Vec<Vec<Scalar>>,
        r: Vec<Scalar>,
        B: Vec<Vec<Scalar>>,
        s: Vec<Scalar>,
        com_ref: CommonRef,
    ) -> Self {
        Self {
            c_Ai: c_Ai,
            c_Bi: c_Bi,
            bi_map: bi_map,
            y: y,
            A: A,
            r: r,
            B: B,
            s: s,
            com_ref: com_ref,
            chall: Challenges{x:Scalar::one(), y:Scalar::one(), z:Scalar::one()},
        }
    }


    ///prove method that creates a ZeroProof
    pub fn prove(
        &mut self, 
        trans: &mut Transcript
    ) -> ZeroProof {
        trans.append_message(b"dom-sep", b"ZeroProof");
        let n: usize = self.A[0].len();
        let m: usize = self.A.len();

        //Get random blinding scalar vectors for the arguments
        let a_0: Vec<Scalar> = (0..n).map(|_| self.com_ref.rand_scalar()).collect();
        let b_m: Vec<Scalar> = (0..n).map(|_| self.com_ref.rand_scalar()).collect();

        let r_0: Scalar = self.com_ref.rand_scalar();;
        let s_m: Scalar = self.com_ref.rand_scalar();

        //Commit the blinding vectors
        let c_A0: RistrettoPoint = self.com_ref.commit(a_0.clone(), r_0.clone());
        let c_Bm: RistrettoPoint = self.com_ref.commit(b_m.clone(), s_m.clone());

        //Blind the argument matrices by appending the vectors, 
        //and their hiding scalars as well
        let blind_A: Vec<Vec<Scalar>> = iter::once(a_0)
            .chain(self.A.clone().into_iter())
            .collect();

        let blind_B: Vec<Vec<Scalar>> = self.B.clone().into_iter()
            .chain(iter::once(b_m))
            .collect();

        let blind_r: Vec<Scalar> = iter::once(r_0)
            .chain(self.r.clone().into_iter())
            .collect();

        let blind_s: Vec<Scalar> = self.s.clone().into_iter()
            .chain(iter::once(s_m))
            .collect();

        //Prepare bilinear map solution Z^2m+1
        let mut d_k: Vec<Scalar> = vec![Scalar::zero(); 2*m + 1];

        //This loop scheme allows diagonal sum of the mapped values,
        //efficiently?
        for i in 0..=m {
            for j in 0..=m {
                let k = m + i -j;

                d_k[k] += (self.bi_map)(
                    blind_A[i].clone(), blind_B[j].clone(), self.y.clone()
                    );
            }
        }

        //Get hiding scalars for the mapped scalars, set m+1 to ZERO
        //Because the m+1th diagonal is the main one(hence the total sum)
        let t : Vec<Scalar> = (0..=2*m).map(
            |i| 
            match i {
                i if i == m+1 => Scalar::zero(),
                _ => self.com_ref.rand_scalar(),
            }).collect();

        //Commit mapped scalars
        let c_D: Vec<RistrettoPoint> = self.com_ref.commit_vec(d_k.clone(), t.clone());

        let x = self.chall.x.clone();

        //Computing blinded openings to end commitments
        let a : Vec<Scalar> = (0..=m).fold(
            vec![Scalar::zero(); n],
            |acc, i|
            acc.add(&blind_A[i].mult(&x.pow(i as u64)))
            );

        let r : Scalar = (0..=m).fold(
            Scalar::zero(),
            |acc, i|
            acc + (blind_r[i] * x.pow(i as u64))
            );

        let b: Vec<Scalar> = (0..=m).fold(
            vec![Scalar::zero(); n],
            |acc, j|
            acc.add(&blind_B[j].mult(&x.pow((m-j) as u64)))
            );

        let s: Scalar = (0..=m).fold(
            Scalar::zero(),
            |acc, j|
            acc + (blind_s[j] * x.pow((m-j) as u64))
            );

        let t_val: Scalar = (0..=2*m).fold(
            Scalar::zero(),
            |acc, k|
            acc + (t[k] * x.pow(k as u64))
            );

        ZeroProof{
            c_A0: c_A0,
            c_Bm: c_Bm,
            c_D: c_D,
            a_vec: a,
            b_vec: b,
            r: r,
            s: s,
            t: t_val,
            x: x,
        }
    }

    ///verify method that verifies the ZeroProof
    pub fn verify(
        &mut self,
        trans: &mut Transcript,
        proof: ZeroProof,
    ) -> Result<(), ProofError> {
        trans.append_message(b"dom-sep", b"ZeroProof");
        let m = (proof.c_D.len() - 1)/ 2;
        let n = proof.a_vec.len();

        assert!((proof.c_D[m+1] == 
                 self.com_ref.commit(vec![Scalar::zero()], Scalar::zero())));

        let x = self.chall.x.clone();
        let x_pow: Vec<Scalar> = (0..=m).map(|i| x.pow((i) as u64)).collect();

        // proof of commitments to A
        let mut commit_A: RistrettoPoint = proof.c_A0;

        for i in 1..=m {
            commit_A += self.c_Ai[i-1] * x_pow[i].clone();
        }
        let open_A = self.com_ref.commit(proof.a_vec.clone(), proof.r);
        assert!((commit_A - open_A).is_identity());

        // proof of commitments to B
        let mut commit_B: RistrettoPoint 
            = (0..=m-1).fold(proof.c_Bm,
                               |acc, j|
                               acc + self.c_Bi[j] * x_pow[m-j].clone()
                            );
        let open_B = self.com_ref.commit(proof.b_vec.clone(), proof.s);
        assert!((commit_B - open_B).is_identity());

        let mut commit_D: RistrettoPoint
            = (1..=2*m).fold(proof.c_D[0],
                            |acc, k|
                            acc + proof.c_D[k] * x.pow((k) as u64)
                            );


        let open_values: Scalar = (self.bi_map)(proof.a_vec, proof.b_vec, self.y.clone());

        let open_D = self.com_ref.commit(vec![open_values], proof.t);

        assert!(commit_D == open_D);
        

        Ok(())
    }
}

#[test]
fn test_c_A() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    
    let mut rng = StdRng::seed_from_u64(2);//from_entropy();
    let m: usize = 13;
    let n: usize = 4;
    let mut com_ref = CommonRef::new(n as u64, rng);
    let a: Vec<Vec<Scalar>> = vec![vec![com_ref.rand_scalar(); m]; n];
    let r: Vec<Scalar> = vec![com_ref.rand_scalar(); m];

    let a_col = a.to_col();

    let c_A: Vec<RistrettoPoint> = com_ref.commit_mat(a.clone(), r.clone());

    let a_0: Vec<Scalar> = (0..n).map(|_| com_ref.rand_scalar()).collect();
    let r_0: Scalar = com_ref.rand_scalar();
    let c_A0: RistrettoPoint = com_ref.commit(a_0.clone(), r_0.clone());

    let blind_A: Vec<Vec<Scalar>> = [&[a_0].as_slice(),
                                        &a_col.clone()[..]]
                                    .concat();

    let blind_r: Vec<Scalar> = [&[r_0].as_slice(),
                                &r.clone()[..]]
                                .concat();

    let x = com_ref.rand_scalar();
    let a_p : Vec<Scalar> = (0..=m).fold(
        vec![Scalar::zero(); n],
        |acc, i|
        acc.add(&blind_A[i].mult(&x.pow(i as u64)))
        );

    let r_p : Scalar = (0..=m).fold(
        Scalar::zero(),
        |acc, i|
        acc + (blind_r[i] * x.pow(i as u64))
        );

    let x_pow: Vec<Scalar> = (0..=m).map(|i| x.pow((i) as u64)).collect();

    let mut commit_A: RistrettoPoint = c_A0;

    for i in 1..=m {
        commit_A += c_A[i-1] * x_pow[i].clone();
    }
    let open_A = com_ref.commit(a_p, r_p);
    assert!((commit_A - open_A).is_identity());
}
#[test]
fn test_base() {
    use crate::provers::hadamard_prover::HadamProver;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    let mut prover_transcript = Transcript::new(b"testZeroProof");

    let mut rng = StdRng::seed_from_u64(2);//from_entropy();
    let m: usize = 4;
    let n: usize = 13;
    let mut com_ref = CommonRef::new((n*m) as u64, rng);
    let a: Vec<Vec<Scalar>> = vec![vec![Scalar::zero(); m]; n];
    let r: Vec<Scalar> = vec![com_ref.rand_scalar(); m];
    let b: Vec<Vec<Scalar>> = vec![vec![com_ref.rand_scalar(); m]; n];
    let s: Vec<Scalar> = vec![com_ref.rand_scalar(); m];

    let y: Scalar = Scalar::random(&mut com_ref.rng);

    //Convert to column matrices for randomness
    let mut a_col = a.to_col();
    let mut b_col = b.to_col();

    a_col[0] = vec![com_ref.rand_scalar(); n];
    b_col[0] = vec![Scalar::zero(); n];

    a_col[m-1] = vec![com_ref.rand_scalar(); n];
    b_col[m-1] = vec![Scalar::zero(); n];

    //Convert back to row matrices
    let a = a_col.to_col();
    let b = b_col.to_col();


    let c_A: Vec<RistrettoPoint> = com_ref.commit_mat(a.clone(), r.clone());
    let c_B: Vec<RistrettoPoint> = com_ref.commit_mat(b.clone(), s.clone());

    let mut zero_prover: ZeroProver = ZeroProver::new(
        c_A,
        c_B,
        HadamProver::exp_dot,
        y,
        a.to_col(),
        r,
        b.to_col(),
        s,
        com_ref.clone()
        );

    let mut zero_proof: ZeroProof = zero_prover.prove(&mut prover_transcript);

    let mut verifier_transcript = Transcript::new(b"testZeroProof");
    assert!(zero_prover.verify(&mut verifier_transcript, zero_proof).is_ok());
        

}

