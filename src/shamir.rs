// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.
use rand;
use numtheory::*;


#[derive(Debug)]
pub struct ShamirSecretSharing {
    pub threshold: usize,  // TODO this is really the reconstruction limit, not threshold
    pub parts: usize,
    pub prime: i64,
}

pub static SHAMIR_5_20: ShamirSecretSharing = ShamirSecretSharing {
    threshold: 5,
    parts: 20,
    prime: 41,
};

impl ShamirSecretSharing {

    pub fn share(&self, secret: i64) -> Vec<i64> {
        let poly = self.sample_polynomial(secret);
        self.evaluate_polynomial(&poly)
    }

    pub fn reconstruct(&self, indices: &[usize], shares: &[i64]) -> i64 {
        assert!(shares.len() == indices.len());
        assert!(shares.len() >= self.threshold);
        // add one to indices to get points
        let points: Vec<i64> = indices.iter().map(|&i| (i as i64) + 1i64).collect();
        lagrange_interpolation_at_zero(&*points, &shares, self.prime)
    }

    fn sample_polynomial(&self, zero_value: i64) -> Vec<i64> {
        // fix the first coefficient (corresponding to the evaluation at zero)
        let mut coefficients = vec![zero_value];
        // sample the remaining coefficients randomly
        //  - use secure randomness as per https://doc.rust-lang.org/rand/rand/index.html#cryptographic-security
        use rand::distributions::Sample;
        let mut range = rand::distributions::range::Range::new(0, self.prime - 1);
        let mut rng = rand::OsRng::new().unwrap();
        let random_coefficients: Vec<i64> = (1..self.threshold).map(|_| range.sample(&mut rng)).collect();
        coefficients.extend(random_coefficients);
        // return
        coefficients
    }

    fn evaluate_polynomial(&self, coefficients: &[i64]) -> Vec<i64> {
        // evaluate at all points
        (1..self.parts + 1)
            .map(|point| mod_evaluate_polynomial(coefficients, point as i64, self.prime))
            .collect()
    }

}


#[test]
fn test_evaluate_polynomial() {
    let ref tss = SHAMIR_5_20;
    let poly = vec![1,2,0];
    let values = tss.evaluate_polynomial(&poly);
    assert_eq!(*values, [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 0]);
}

#[test]
fn wikipedia_example() {
    let tss = ShamirSecretSharing {
        threshold: 3,
        parts: 6,
        prime: 1613
    };

    let shares = tss.evaluate_polynomial(&[1234, 166, 94]);
    assert_eq!(&*shares, &[1494, 329, 965, 176, 1188, 775]);

    assert_eq!(tss.reconstruct(&[0, 1, 2], &shares[0..3]), 1234);
    assert_eq!(tss.reconstruct(&[1, 2, 3], &shares[1..4]), 1234);
    assert_eq!(tss.reconstruct(&[2, 3, 4], &shares[2..5]), 1234);
}

#[test]
fn test_shamir() {
    let tss = ShamirSecretSharing {
        threshold: 3,
        parts: 6,
        prime: 41
    };
    let secret = 1;
    let shares = tss.share(secret);
    assert_eq!(tss.reconstruct(&[0, 1, 2],    &shares[0..3]), secret);
    assert_eq!(tss.reconstruct(&[1, 2, 3],    &shares[1..4]), secret);
    assert_eq!(tss.reconstruct(&[2, 3, 4, 5], &shares[2..6]), secret);
}
