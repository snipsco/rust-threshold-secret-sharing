
pub fn gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let n = a / b;
        let c = a % b;
        let r = gcd(b, c);
        (r.0, r.2, r.1 - r.2 * n)
    }
}

#[test]
fn test_gcd() {
    assert_eq!(gcd(12, 16), (4, -1, 1));
}


pub fn mod_inverse(k: i64, prime: i64) -> i64 {
    let k2 = k % prime;
    let r = if k2 < 0 {
        -gcd(prime, -k2).2
    } else {
        gcd(prime, k2).2
    };
    (prime + r) % prime
}

#[test]
fn test_mod_inverse() {
    assert_eq!(mod_inverse(3, 7), 5);
}


pub fn mod_pow(x: i64, e: u32, prime: i64) -> i64 {
    // TODO optmise by repeated-squaring
    x.pow(e) % prime
}


// TODO don't insist on vectors as input (slices would be better)
pub fn fft2(a_coef: Vec<i64>, omega: i64, prime: i64) -> Vec<i64> {
    if a_coef.len() == 1 {
        a_coef
    } else {
        // split a into b and c
        // TODO avoid copying
        let b_coef: Vec<i64> = a_coef.iter().enumerate().filter_map(|(x,&i)| if x % 2 == 0 { Some(i) } else { None } ).collect();
        let c_coef: Vec<i64> = a_coef.iter().enumerate().filter_map(|(x,&i)| if x % 2 == 1 { Some(i) } else { None } ).collect();

        // recurse
        let b_point = fft2(b_coef, mod_pow(omega, 2, prime), prime);
        let c_point = fft2(c_coef, mod_pow(omega, 2, prime), prime);

        // combine
        let len = a_coef.len();
        let half_len = len >> 1;
        let mut a_point = vec![0; len];  // TODO unsafe { Vec.set_len() } trick
        for i in 0..half_len {
            a_point[i]            = (b_point[i] + mod_pow(omega, i as u32, prime) * c_point[i]) % prime;
            a_point[i + half_len] = (b_point[i] - mod_pow(omega, i as u32, prime) * c_point[i]) % prime;
        }

        // return
        a_point
    }
}

pub fn fft2_inverse(a_point: Vec<i64>, omega: i64, prime: i64) -> Vec<i64> {
    let omega_inv = mod_inverse(omega, prime);
    let len = a_point.len();
    let len_inv = mod_inverse(len as i64, prime);
    let scaled_a_coef = fft2(a_point, omega_inv, prime);
    let a_coef = scaled_a_coef.iter().map(|x| { x * len_inv % prime }).collect();
    a_coef
}

#[test]
fn test_fft2() {
    // field is Z_433 in which 354 is an 8th root of unity
    let prime = 433;
    let omega = 354;

    let a_coef = vec![1,2,3,4,5,6,7,8];
    let a_point = fft2(a_coef, omega, prime);
    assert_eq!(a_point, vec![36, -130, -287, 3, -4, 422, 279, -311])
}

#[test]
fn test_fft2_inverse() {
    // field is Z_433 in which 354 is an 8th root of unity
    let prime = 433;
    let omega = 354;

    let a_point = vec![36, -130, -287, 3, -4, 422, 279, -311];
    let a_coef = fft2_inverse(a_point, omega, prime);
    assert_eq!(a_coef, vec![1, 2, 3, -429, 5, -427, -426, 8])
}


// TODO don't insist on vectors as input (slices would be better)
pub fn fft3(a_coef: Vec<i64>, omega: i64, prime: i64) -> Vec<i64> {
    if a_coef.len() == 1 {
        a_coef
    } else {
        // split a into b and c
        // TODO avoid copying
        let b_coef: Vec<i64> = a_coef.iter().enumerate().filter_map(|(x,&i)| if x % 3 == 0 { Some(i) } else { None } ).collect();
        let c_coef: Vec<i64> = a_coef.iter().enumerate().filter_map(|(x,&i)| if x % 3 == 1 { Some(i) } else { None } ).collect();
        let d_coef: Vec<i64> = a_coef.iter().enumerate().filter_map(|(x,&i)| if x % 3 == 2 { Some(i) } else { None } ).collect();

        // recurse
        let omega_cubed = mod_pow(omega, 3, prime);
        let b_point = fft3(b_coef, omega_cubed, prime);
        let c_point = fft3(c_coef, omega_cubed, prime);
        let d_point = fft3(d_coef, omega_cubed, prime);

        // combine
        let len = a_coef.len();
        let third_len = len / 3;
        let mut a_point = vec![0; len];  // TODO unsafe { Vec.set_len() } trick
        for i in 0..third_len {

            let j = i;
            let x = mod_pow(omega, j as u32, prime);
            let x_squared = (x * x) % prime;
            a_point[j] = (b_point[i] + x * c_point[i] + x_squared * d_point[i]) % prime;

            let j = i + third_len;
            let x = mod_pow(omega, j as u32, prime);
            let x_squared = (x * x) % prime;
            a_point[j] = (b_point[i] + x * c_point[i] + x_squared * d_point[i]) % prime;

            let j = i + third_len + third_len;
            let x = mod_pow(omega, j as u32, prime);
            let x_squared = (x * x) % prime;
            a_point[j] = (b_point[i] + x * c_point[i] + x_squared * d_point[i]) % prime;
        }

        // return
        a_point
    }
}

pub fn fft3_inverse(a_point: Vec<i64>, omega: i64, prime: i64) -> Vec<i64> {
    let omega_inv = mod_inverse(omega, prime);
    let len = a_point.len();
    let len_inv = mod_inverse(len as i64, prime);
    let scaled_a_coef = fft3(a_point, omega_inv, prime);
    let a_coef = scaled_a_coef.iter().map(|x| { x * len_inv % prime }).collect();
    a_coef
}

#[test]
fn test_fft3() {
    // field is Z_433 in which 150 is an 9th root of unity
    let prime = 433;
    let omega = 150;

    let a_coef = vec![1,2,3,4,5,6,7,8,9];
    let a_point = fft3(a_coef, omega, prime);
    assert_eq!(a_point, vec![45, 404, 407, 266, 377, 47, 158, 17, 20])
}

#[test]
fn test_fft3_inverse() {
    // field is Z_433 in which 150 is an 9th root of unity
    let prime = 433;
    let omega = 150;

    let a_point = vec![45, 404, 407, 266, 377, 47, 158, 17, 20];
    let a_coef = fft3_inverse(a_point, omega, prime);
    assert_eq!(a_coef, vec![1,2,3,4,5,6,7,8,9])
}
