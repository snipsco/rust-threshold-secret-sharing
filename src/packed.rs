
// use super::{Secret, Share};
use numtheory::{mod_pow, fft2_inverse, fft3};
use rand;


#[derive(Debug,Copy,Clone)]
pub struct PackedSecretSharing {
    pub prime: i64,
    pub threshold: usize,
    pub share_count: usize,
    pub secret_count: usize,
    pub reconstruct_limit: usize, // = n
    n: usize,  // n = secret_count + threshold + 1
    m: usize,  // m = share_count + 1
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

    pub fn share(&self, secrets: &[i64]) -> Vec<i64> {
        assert_eq!(secrets.len(), self.secret_count);
        // sample polynomial
        let mut poly = self.sample_polynomial(secrets);
        // .. and extend it
        poly.extend( vec![0; self.m - self.n] );
        // evaluate polynomial to generate shares
        let mut shares = self.evaluate_polynomial(poly);
        // .. but remove first element since it should not be used as a share (it's always 1)
        shares.remove(0);
        // return
        assert_eq!(shares.len(), self.share_count);
        shares
    }

    fn sample_polynomial(&self, secrets: &[i64]) -> Vec<i64> {
        // sample randomness
        //  - for cryptographic use we should use OsRng as dictated here
        //    https://doc.rust-lang.org/rand/rand/index.html#cryptographic-security
        use rand::distributions::Sample;
        let mut range = rand::distributions::range::Range::new(0, self.prime - 1);
        let mut rng = rand::OsRng::new().unwrap();
        let randomness: Vec<i64> = (0..self.threshold).map(|_| range.sample(&mut rng) as i64).collect();
        // recover polynomial
        let coefficients = self.recover_polynomial(secrets, randomness);
        coefficients
    }

    fn recover_polynomial(&self, secrets: &[i64], randomness: Vec<i64>) -> Vec<i64> {
        // fix the value corresponding to point 1
        let mut values: Vec<i64> = vec![0];
        // let the subsequent values correspond to the secrets
        values.extend(secrets);
        // fill in with random values
        values.extend(randomness);
        // run backward FFT to recover polynomial in coefficient representation
        assert_eq!(values.len(), self.n);
        let coefficients = fft2_inverse(&values, self.omega_n, self.prime);
        coefficients
    }

    fn evaluate_polynomial(&self, coefficients: Vec<i64>) -> Vec<i64> {
        assert_eq!(coefficients.len(), self.m);
        let points = fft3(&coefficients, self.omega_m, self.prime);
        points
    }

    pub fn reconstruct(&self, indices: &[usize], shares: &[i64]) -> Vec<i64> {
        assert_eq!(shares.len(), indices.len());
        assert!(shares.len() >= self.reconstruct_limit);
        let shares_points: Vec<i64> = indices.iter().map(|&x| mod_pow(self.omega_m, x as u32 + 1, self.prime)).collect();
        // interpolate using Newton's method
        use numtheory::{newton_interpolation_general, newton_evaluate};
        // TODO optimise by using Newton-equally-space variant
        let poly = newton_interpolation_general(&shares_points, &shares, self.prime);
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
    let mut shares = pss.share(&secrets);

    // manually recover secrets
    use numtheory::{fft3_inverse, mod_evaluate_polynomial};
    shares.insert(0, 0);
    let poly = fft3_inverse(&shares, PSS_4_26_3.omega_m, PSS_4_26_3.prime);
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
    assert_eq!(shares.len(), pss.share_count);
}

#[test]
fn test_share_reconstruct() {
    let ref pss = PSS_4_26_3;
    let secrets = vec![5,6,7];
    let shares = pss.share(&secrets);

    use numtheory::positivise;

    // reconstruction must work for all shares
    let indices: Vec<usize> = (0..shares.len()).collect();
    let recovered_secrets = pss.reconstruct(&indices, &shares);
    assert_eq!(positivise(&recovered_secrets, pss.prime), secrets);

    // .. and for only sufficient shares
    let indices: Vec<usize> = (0..pss.reconstruct_limit).collect();
    let recovered_secrets = pss.reconstruct(&indices, &shares[0..pss.reconstruct_limit]);
    assert_eq!(positivise(&recovered_secrets, pss.prime), secrets);
}
