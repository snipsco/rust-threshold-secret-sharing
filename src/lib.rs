extern crate rand;

#[derive(Debug)]
pub struct ThresholdSecretSharing {
    threshold: usize,
    parts: usize,
    prime: u64,
    exps: Vec<Vec<u64>>,
}

impl ThresholdSecretSharing {
    pub fn new(parts: usize, threshold: usize, prime: u64) -> ThresholdSecretSharing {
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
        ThresholdSecretSharing {
            threshold: threshold,
            parts: parts,
            prime: prime,
        exps: exps,
    }
    }

    #[allow(dead_code)]
    fn pow(x: u64, e: u64, prime: u64) -> u64 {
        (0..e).fold(1, |a, _| (a * x) % prime)
    }

    fn optimized_pow(&self, x: u64, e: u64) -> u64 {
        self.exps[x as usize - 1][e as usize]
    }

    fn encrypt_with_noise(&self, secret: u64, noise: Vec<u64>) -> Vec<(u64, u64)> {
        (1..self.parts + 1)
            .map(|x| x as u64)
            .map(|x| {
                (x,
                 noise.iter()
                    .enumerate()
                    .map(|(deg, coef)| coef * self.optimized_pow(x, deg as u64 + 1))
                    .fold(secret, |a, b| (a + b) % self.prime))
            })
            .collect()
    }

    pub fn encrypt(&self, secret: u64) -> Vec<(u64, u64)> {
        use rand::distributions::Sample;
        let mut range = rand::distributions::range::Range::new(0, self.prime - 1);
        let mut rng = rand::thread_rng();
        let noise = (1..self.threshold).map(|_| range.sample(&mut rng)).collect();
        self.encrypt_with_noise(secret, noise)
    }

    fn gcd(a: i64, b: i64) -> (i64, i64, i64) {
        if b == 0 {
            (a, 1, 0)
        } else {
            let n = a / b;
            let c = a % b;
            let r = ThresholdSecretSharing::gcd(b, c);
            (r.0, r.2, r.1 - r.2 * n)
        }
    }

    fn mod_inverse(k: i64, prime: i64) -> i64 {
        let k2 = k % prime;
        let r = if k2 < 0 {
            -Self::gcd(prime, -k2).2
        } else {
            Self::gcd(prime, k2).2
        };
        (prime + r) % prime
    }

    pub fn decrypt(&self, shares: &[(u64, u64)]) -> u64 {
        assert_eq!(shares.len(), self.threshold);
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
                   (shares[j].1 as i64 * num) % prime * Self::mod_inverse(denum, prime)) %
                  prime;
        }
        acc as u64
    }
}

#[test]
fn wikipedia_example() {
    let secret = 1234;
    let prime = 1613;
    let threshold = 3;
    let parts = 6;

    let algo = ThresholdSecretSharing::new(parts, threshold, prime);

    let shares = algo.encrypt_with_noise(1234, vec![166, 94]);
    assert_eq!(&*shares,
               &[(1, 1494), (2, 329), (3, 965), (4, 176), (5, 1188), (6, 775)]);

    assert_eq!(algo.decrypt(&shares[0..3]), secret);
    assert_eq!(algo.decrypt(&shares[1..4]), secret);
    assert_eq!(algo.decrypt(&shares[2..5]), secret);
}

#[test]
fn test_pow() {
    assert_eq!(ThresholdSecretSharing::pow(2, 4, 11), 5);
}

#[test]
fn test_optimized_pow() {
    let sss = ThresholdSecretSharing::new(10, 5, 8191);
    assert_eq!(ThresholdSecretSharing::pow(2, 4, 8191),
               sss.optimized_pow(2, 4));
}

#[test]
fn test_gcd() {
    assert_eq!(ThresholdSecretSharing::gcd(12, 16), (4, -1, 1));
}

#[test]
fn test_mod_inverse() {
    assert_eq!(ThresholdSecretSharing::mod_inverse(3, 7), 5);
}
