use curve25519_dalek::ristretto::CompressedRistretto;
use curve25519_dalek::scalar::Scalar;
use merlin::Transcript;

use bulletproofs::ProofError;

pub trait TranscriptProtocol {
    /// Append a domain seperator for an m*n = N card, shuffle proof
    fn shuffle_domain_sep(&mut self, n: u64, m: u64);
    /// Append a domain seperator for an m = mu * m' reduction, Multi-expo proof
    fn mexp_domain_sep(&mut self, m: u64, mu: u64);
    /// Append a domain seperator for an nxm sized, product proof
    fn prod_domain_sep(&mut self, n: u64, m: u64);
    /// Append a `scalar` with the given `label`.

    fn append_scalar(&mut self, label: &'static [u8], scalar: &Scalar);
    fn append_scalar_vec(&mut self, label: &'static [u8], scalar: &Vec<Scalar>);

    /// Append a `point` with the given `label`.
    fn append_point(&mut self, label: &'static [u8], point: &CompressedRistretto);
    fn append_point_vec(&mut self, label: &'static [u8], points: &Vec<CompressedRistretto>);

    fn val_append_point_vec(
        &mut self,
        label: &'static [u8],
        points: &Vec<CompressedRistretto>,
    ) -> Result<(), ProofError>;
    /// Check that a point is not the identity, then append it to the
    /// transcript.  Otherwise, return an error.
    fn validate_and_append_point(
        &mut self,
        label: &'static [u8],
        point: &CompressedRistretto,
    ) -> Result<(), ProofError>;

    /// Compute a `label`ed challenge variable.
    fn challenge_scalar(&mut self, label: &'static [u8]) -> Scalar;

}

impl TranscriptProtocol for Transcript {
    fn shuffle_domain_sep(&mut self, n: u64, m: u64) {
        self.append_message(b"dom-sep", b"Shuffleproof");
        self.append_u64(b"n", n);
        self.append_u64(b"m", m);
    }

    fn mexp_domain_sep(&mut self, m: u64, mu: u64) {
        self.append_message(b"dom-sep", b"MExpproof");
        self.append_u64(b"m_mexp", m);
        self.append_u64(b"mu_mexp", mu);
    }

    fn prod_domain_sep(&mut self, n: u64, m: u64) {
        self.append_message(b"dom-sep", b"Productproof");
        self.append_u64(b"n_prod", n);
        self.append_u64(b"m_prod", m);
    }
     fn append_scalar(&mut self, label: &'static [u8], scalar: &Scalar) {
        self.append_message(label, scalar.as_bytes());
    }
     fn append_scalar_vec(&mut self, label: &'static [u8], scalar: &Vec<Scalar>) {
        let mut res_vec = vec![];
        let res_bytes = &scalar.iter()
                                    .fold(res_vec,
                                          |acc:Vec<u8>, point:&Scalar| 
                                          [acc.as_slice(), point.as_bytes()].concat()
                                    );

        println!("{:?}", res_bytes);
        self.append_message(label, res_bytes);
    }

    fn append_point(&mut self, label: &'static [u8], point: &CompressedRistretto) {
        self.append_message(label, point.as_bytes());
    }

    fn append_point_vec(&mut self, label: &'static [u8], points: &Vec<CompressedRistretto>) {
        let mut res_vec = vec![];
        let res_bytes = &points.iter()
                                    .fold(res_vec,
                                          |acc:Vec<u8>, point:&CompressedRistretto| 
                                          [acc.as_slice(), point.as_bytes()].concat()
                                    );

        self.append_message(label, res_bytes);
    }
    fn val_append_point_vec(
        &mut self,
        label: &'static [u8],
        points: &Vec<CompressedRistretto>,
    ) -> Result<(), ProofError> {
        use curve25519_dalek::traits::IsIdentity;

        for point in points.iter() {
            if point.is_identity() {
                return Err(ProofError::VerificationError);
            }
        }

        self.append_point_vec(label, points);
        Ok(())
    }

    fn validate_and_append_point(
        &mut self,
        label: &'static [u8],
        point: &CompressedRistretto,
    ) -> Result<(), ProofError> {
        use curve25519_dalek::traits::IsIdentity;

        if point.is_identity() {
            Err(ProofError::VerificationError)
        } else {
            Ok(self.append_message(label, point.as_bytes()))
        }
    }

    fn challenge_scalar(&mut self, label: &'static [u8]) -> Scalar {
        let mut buf = [0u8; 64];
        self.challenge_bytes(label, &mut buf);

        Scalar::from_bytes_mod_order_wide(&buf)
    }

}
