
// use super::{Secret, Share};
use numtheory::{mod_pow, fft2_inverse, fft3};
use rand;


#[derive(Debug)]
pub struct PackedSecretSharing {
    threshold: usize,
    share_count: usize,
    secret_count: usize,
    reconstruct_limit: usize, // = n
    n: usize,  // n = secret_count + threshold + 1
    m: usize,  // m = share_count + 1
    prime: i64,
    omega_n: i64,
    omega_m: i64
}

pub static PSS_4_8_3: PackedSecretSharing = PackedSecretSharing {
    threshold: 4,
    share_count: 8,
    secret_count: 3,
    reconstruct_limit: 8,
    n: 8, // 3 + 4 + 1
    m: 9, // 8 + 1
    prime: 433,
    omega_n: 354,
    omega_m: 150,
};

pub static PSS_4_26_3: PackedSecretSharing = PackedSecretSharing {
    threshold: 4,
    share_count: 26,
    secret_count: 3,
    reconstruct_limit: 8,
    n: 8,  // 3 + 4 + 1
    m: 27, // 26 + 1
    prime: 433,
    omega_n: 354,
    omega_m: 17,
};

pub static PSS_155_728_100: PackedSecretSharing = PackedSecretSharing {
    threshold: 155,
    share_count: 728,
    secret_count: 100,
    reconstruct_limit: 256,
    n: 256, // 100 + 155 + 1
    m: 729, // 728 + 1
    prime: 746497,
    omega_n: 95660,
    omega_m: 610121,
};

pub static PSS_155_19682_100: PackedSecretSharing = PackedSecretSharing {
    threshold: 155,
    share_count: 19682,
    secret_count: 100,
    reconstruct_limit: 256,
    n: 256,   // 100 + 155 + 1
    m: 19683, // 19682 + 1
    prime: 5038849,
    omega_n: 4318906,
    omega_m: 1814687,
};


impl PackedSecretSharing {

    pub fn share(&self, secrets: &Vec<i64>) -> Vec<(i64, i64)> {
        assert_eq!(secrets.len(), self.secret_count);
        // sample polynomial
        let mut poly = self.sample_polynomial(secrets);
        // .. extend it
        poly.extend( vec![0; self.m - self.n] );
        // .. evaluate it
        let values = self.evaluate_polynomial(poly);
        // .. zip the resulting values with the evaluation points
        let points = (0..self.m).map(|e| mod_pow(self.omega_m, (e as u32), self.prime));
        let shares: Vec<_> = points.zip(values).collect();
        // .. and return the pairs
        assert_eq!(shares.len(), self.m);
        shares // TODO ignore first share (it's always zero)
    }

    fn sample_polynomial(&self, secrets: &Vec<i64>) -> Vec<i64> {
        // sample randomness
        use rand::distributions::Sample;
        let mut range = rand::distributions::range::Range::new(0, self.prime - 1);
        let mut rng = rand::thread_rng();  // TODO use OsRng as dictated here https://doc.rust-lang.org/rand/rand/index.html#cryptographic-security
        let randomness: Vec<i64> = (0..self.threshold).map(|_| range.sample(&mut rng) as i64).collect();
        // recover polynomial
        let coefficients = self.recover_polynomial(secrets, randomness);
        coefficients
    }

    fn recover_polynomial(&self, secrets: &Vec<i64>, randomness: Vec<i64>) -> Vec<i64> {
        // fix the value corresponding to point 1
        let mut values: Vec<i64> = vec![0];
        // let the subsequent values correspond to the secrets
        values.extend(secrets);
        // fill in with random values
        values.extend(randomness);
        // run backward FFT to recover polynomial in coefficient representation
        assert_eq!(values.len(), self.n);
        let coefficients = fft2_inverse(values, self.omega_n, self.prime);
        coefficients
    }

    fn evaluate_polynomial(&self, coefficients: Vec<i64>) -> Vec<i64> {
        assert_eq!(coefficients.len(), self.m);
        let points = fft3(coefficients, self.omega_m, self.prime);
        points
    }

    pub fn reconstruct(&self, shares: &[(i64, i64)]) -> Vec<i64> {
        assert!(shares.len() >= self.reconstruct_limit);
        // TODO better way to do unzip?
        let shares_points: Vec<i64> = shares.iter().map(|x| x.0).collect();
        let shares_values: Vec<i64> = shares.iter().map(|x| x.1).collect();
        // interpolate using Newton's method
        use numtheory::{newton_interpolation_general, newton_evaluate};
        // TODO optimise by using Newton-equally-space variant
        let poly = newton_interpolation_general(&shares_points, &shares_values, self.prime);
        // evaluate at omega_n points to recover secrets
        // TODO optimise to avoid re-computation of power
        let secrets = (1..self.n)
            .map(|e| mod_pow(self.omega_n, e as u32, self.prime))
            .map(|point| newton_evaluate(&poly, point, self.prime))
            .take(self.secret_count)
            .collect();
        secrets
    }

}


#[test]
fn test_recover_polynomial() {
    let ref pss = PSS_4_8_3;
    let secrets = vec![1,2,3];
    let randomness = vec![8,8,8,8];  // use fixed randomness
    let poly = pss.recover_polynomial(&secrets, randomness);
    assert_eq!(poly, vec![113, -382, -172, 267, -325, 432, 388, -321]);
}

#[test]
fn test_evaluate_polynomial() {
    let ref pss = PSS_4_26_3;
    let poly = vec![113, 51, 261, 267, 108, 432, 388, 112, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let points = pss.evaluate_polynomial(poly);
    assert_eq!(points, vec![0, 77, 230, 91, 286, 179, 337, 83, 212, 88, 406, 58, 425, 345, 350, 336, 430, 404, 51, 60, 305, 395, 84, 156, 160, 112, 422]);
}

#[test]
fn test_share() {
    let ref pss = PSS_4_26_3;

    // do sharing
    let secrets = vec![5,6,7];
    let shares = pss.share(&secrets);

    // manually recover secrets
    let values: Vec<i64> = shares.iter().map({|&(_, v)| v}).collect();
    use numtheory::{fft3_inverse, mod_evaluate_polynomial};
    let poly = fft3_inverse(values, PSS_4_26_3.omega_m, PSS_4_26_3.prime);
    let recovered_secrets: Vec<i64> = (1..secrets.len()+1)
        .map(|i| mod_evaluate_polynomial(&poly, mod_pow(PSS_4_26_3.omega_n, i as u32, PSS_4_26_3.prime), PSS_4_26_3.prime))
        .collect();

    use numtheory::positivise;
    assert_eq!(positivise(&recovered_secrets, pss.prime), secrets);
}

#[test]
fn test_large_share() {
    let ref pss = PSS_155_19682_100;
    let secrets = vec![5 ; pss.secret_count];
    let shares = pss.share(&secrets);
    assert_eq!(shares.len(), pss.share_count+1); // TODO remove 1
}

#[test]
fn test_share_reconstruct() {
    let ref pss = PSS_4_26_3;
    let secrets = vec![5,6,7];
    let shares = pss.share(&secrets);

    use numtheory::positivise;

    // reconstruction must work for all shares
    let all_shares = shares.clone();
    let recovered_secrets = pss.reconstruct(&all_shares);
    assert_eq!(positivise(&recovered_secrets, pss.prime), secrets);

    // .. and for only sufficient shares
    let sufficient_shares: Vec<(i64, i64)> = shares.clone().into_iter().take(pss.reconstruct_limit).collect();
    let recovered_secrets = pss.reconstruct(&sufficient_shares);
    assert_eq!(positivise(&recovered_secrets, pss.prime), secrets);
}
