use rand;

use super::{Secret, Share};
use numtheory::*;


#[derive(Debug)]
pub struct SingleThresholdSecretSharing {
    threshold: usize,
    parts: usize,
    prime: u64,
    exps: Vec<Vec<u64>>,
}


impl SingleThresholdSecretSharing {

    pub fn new(parts: usize, threshold: usize, prime: u64) -> SingleThresholdSecretSharing {
        let exps = (1..parts + 1)
            .map(|x| {
                let mut column = Vec::with_capacity(threshold);
                column.push(1);
                while column.len() < threshold {
                    let previous = *column.last().unwrap();
                    column.push((previous * x as u64) % prime);
                }
                column
            })
            .collect();
        SingleThresholdSecretSharing {
            threshold: threshold,
            parts: parts,
            prime: prime,
            exps: exps,
        }
    }

    pub fn share(&self, secret: Secret) -> Vec<Share> {
        let coefficients = self.sample_polynomial(secret);
        self.evaluate_polynomial(coefficients)
    }

    pub fn reconstruct(&self, shares: &[Share]) -> Secret {
        let (shares, _) = shares.split_at(self.threshold);
        let prime = self.prime as i64;
        let mut acc = 0i64;
        for j in 0..self.threshold {
            let xj = shares[j].0 as i64;
            let mut num = 1i64;
            let mut denum = 1i64;
            for m in 0..self.threshold {
                if j != m {
                    let xm = shares[m].0 as i64;
                    num = (num * xm) % prime;
                    denum = (denum * (xm - xj)) % prime;
                }
            }
            acc = (prime + acc +
                   (shares[j].1 as i64 * num) % prime * mod_inverse(denum, prime)) %
                  prime;
        }
        acc as u64
    }

    fn sample_polynomial(&self, zero: u64) -> Vec<u64> {
        // fix the first coefficient
        let mut coefficients = vec![zero];
        // sample the remaining coefficients randomly
        use rand::distributions::Sample;
        let mut range = rand::distributions::range::Range::new(0, self.prime - 1);
        let mut rng = rand::thread_rng();
        let random_coefficients: Vec<u64> = (1..self.threshold).map(|_| range.sample(&mut rng)).collect();
        coefficients.extend(random_coefficients);
        // return
        return coefficients;
    }

    fn evaluate_polynomial(&self, coefficients: Vec<u64>) -> Vec<(u64, u64)> {
        // TODO optimise with Horner's rule
        (1..self.parts + 1)
            .map(|x| x as u64)
            .map(|x| {
                (x,
                 coefficients.iter()
                    .enumerate()
                    .map(|(deg, coef)| coef * self.optimized_pow(x, deg as u64))
                    .fold(0, |a, b| (a + b) % self.prime))
            })
            .collect()
    }

    #[allow(dead_code)]
    fn pow(x: u64, e: u64, prime: u64) -> u64 {
        (0..e).fold(1, |a, _| (a * x) % prime)
    }

    fn optimized_pow(&self, x: u64, e: u64) -> u64 {
        if e == 0 {
            1
        } else {
            self.exps[x as usize - 1][e as usize]
        }
    }

}




#[test]
fn single() {
    let tss = SingleThresholdSecretSharing::new(20, 5, 41);

    let secret = 1;
    let shares = tss.share(secret);
    let recon_secret = tss.reconstruct(&*shares);
    println!("{:?}", recon_secret);
}

// #[test]
// fn multi() {
//
//     use super::MultiThresholdSecretSharing;
//
//     let prime = 41;
//     let threshold = 5;
//     let shares = 20;
//     let tss = ThresholdSecretSharing::new(shares, threshold, prime);
//
//
//     let secrets = vec![1, 2, 3];
//     let shares = tss.share(secrets);
//
//     let recon_secrets = tss.reconstruct(shares);
//
//     assert_eq!(secrets, recon_secrets);
// }










#[test]
fn wikipedia_example() {
    let secret = 1234;
    let prime = 1613;
    let threshold = 3;
    let parts = 6;

    let algo = SingleThresholdSecretSharing::new(parts, threshold, prime);

    let shares = algo.evaluate_polynomial(vec![1234, 166, 94]);
    assert_eq!(&*shares,
               &[(1, 1494), (2, 329), (3, 965), (4, 176), (5, 1188), (6, 775)]);

    assert_eq!(algo.reconstruct(&shares[0..3]), secret);
    assert_eq!(algo.reconstruct(&shares[1..4]), secret);
    assert_eq!(algo.reconstruct(&shares[2..5]), secret);
}

#[test]
fn test_pow() {
    assert_eq!(SingleThresholdSecretSharing::pow(2, 4, 11), 5);
}

#[test]
fn test_optimized_pow() {
    let tss = SingleThresholdSecretSharing::new(10, 5, 8191);
    assert_eq!(SingleThresholdSecretSharing::pow(2, 4, 8191), tss.optimized_pow(2, 4));
}
