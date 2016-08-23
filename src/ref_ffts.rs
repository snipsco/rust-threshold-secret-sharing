use numtheory::{mod_pow, mod_inverse};

#[cfg(test)]
use numtheory::positivise;

pub fn fft2(a_coef: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    fft2_in_place(a_coef, omega, prime)
}

pub fn fft3(a_coef: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    fft3_in_place(a_coef, omega, prime)
}

macro_rules! maybe_reduce {
    ($target:expr, $prime:ident, $max:ident) => {
        if $target > $max {
            $target %= $prime;
        }
    }
}

pub fn fft2_rec(a_coef: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    if a_coef.len() == 1 {
        a_coef.to_vec()
    } else {
        // split A(x) into B(x) and C(x): A(x) = B(x^2) + x C(x^2)
        // TODO avoid copying
        let b_coef: Vec<i64> = a_coef.iter()
            .enumerate()
            .filter_map(|(x, &i)| if x % 2 == 0 { Some(i) } else { None })
            .collect();
        let c_coef: Vec<i64> = a_coef.iter()
            .enumerate()
            .filter_map(|(x, &i)| if x % 2 == 1 { Some(i) } else { None })
            .collect();

        // recurse
        let b_point = fft2_rec(&b_coef, mod_pow(omega, 2, prime), prime);
        let c_point = fft2_rec(&c_coef, mod_pow(omega, 2, prime), prime);

        // combine
        let len = a_coef.len();
        let half_len = len >> 1;
        let mut a_point = vec![0; len];  // TODO trick: unsafe { Vec.set_len() }
        for i in 0..half_len {
            a_point[i] = (b_point[i] + mod_pow(omega, i as u32, prime) * c_point[i]) % prime;
            a_point[i + half_len] = (b_point[i] - mod_pow(omega, i as u32, prime) * c_point[i]) %
                                    prime;
        }

        a_point
    }
}

pub fn fft2_stride(a_coef: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    fft2_rec_stride(a_coef, 0, 1, omega, prime)
}

fn fft2_rec_stride(a_coef: &[i64],
                   offset: usize,
                   stride: usize,
                   omega: i64,
                   prime: i64)
                   -> Vec<i64> {
    if a_coef.len() == stride {
        vec![a_coef[offset]]
    } else {
        // recurse
        let b_point = fft2_rec_stride(&a_coef, offset, 2 * stride, mod_pow(omega, 2, prime), prime);
        let c_point = fft2_rec_stride(&a_coef,
                                      offset + stride,
                                      2 * stride,
                                      mod_pow(omega, 2, prime),
                                      prime);

        // combine
        let len = a_coef.len() / stride;
        let half_len = len >> 1;
        // let mut a_point = vec![0; len];  // TODO trick: unsafe { Vec.set_len() }
        let mut a_point = Vec::with_capacity(len);  // TODO trick: unsafe { Vec.set_len() }
        unsafe {
            a_point.set_len(len);
        }
        for i in 0..half_len {
            a_point[i] = (b_point[i] + mod_pow(omega, i as u32, prime) * c_point[i]) % prime;
            a_point[i + half_len] = (b_point[i] - mod_pow(omega, i as u32, prime) * c_point[i]) %
                                    prime;
        }

        a_point
    }
}

pub fn fft2_in_place_rearrange(data: &mut [u64]) {
    let mut target = 0;
    for pos in 0..data.len() {
        if target > pos {
            data.swap(target, pos)
        }
        let mut mask = data.len() >> 1;
        while target & mask != 0 {
            target &= !mask;
            mask >>= 1;
        }
        target |= mask;
    }
}

pub fn fft2_in_place_compute(data: &mut [u64], omega: u64, prime: u64) {
    let mut depth = 0;
    let max = u64::max_value() / prime / 2;
    let big_omega = mod_pow(omega as i64, data.len() as u32 / 2, prime as i64) as u64 % prime;
    while 1 << depth < data.len() {
        let step = 1 << depth;
        let jump = 2 * step;
        let factor_stride =
            mod_pow(omega as i64, (data.len() / step / 2) as u32, prime as i64) as u64;
        let mut factor = 1;
        for group in 0..step {
            let mut pair = group;
            while pair < data.len() {
                let (x, mut y) = (data[pair], data[pair + step] * factor);
                maybe_reduce!(y, prime, max);

                data[pair] = x + y;
                data[pair + step] = x + big_omega * y;

                maybe_reduce!(data[pair], prime, max);
                maybe_reduce!(data[pair + step], prime, max);

                pair += jump;
            }
            factor = (factor * factor_stride) % prime;
        }
        depth += 1;
    }
}

pub fn fft2_in_place(a_coef: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    let mut data: Vec<u64> =
        a_coef.iter().map(|i| (if *i > 0 { *i } else { i + prime }) as u64).collect();
    fft2_in_place_rearrange(&mut *data);
    fft2_in_place_compute(&mut *data, omega as u64, prime as u64);
    data.iter().map(|coef| (coef % prime as u64) as i64).collect()
}

#[test]
fn test_fft2_in_place_rearrange() {
    let mut input = [0, 1, 2, 3, 4, 5, 6, 7];
    fft2_in_place_rearrange(&mut input);
    assert_eq!(input, [0, 4, 2, 6, 1, 5, 3, 7]);
}

/// Inverse FFT for `fft2`.
pub fn fft2_inverse(a_point: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    let omega_inv = mod_inverse(omega, prime);
    let len = a_point.len();
    let len_inv = mod_inverse(len as i64, prime);
    let scaled_a_coef = fft2(a_point, omega_inv, prime);
    let a_coef = scaled_a_coef.iter().map(|x| x * len_inv % prime).collect();
    a_coef
}

#[test]
fn test_fft2_variants() {
    let prime = 433;
    let omega = (354 * 354) % 433;
    for example in &[vec![0, 0, 0, 0],
                     vec![1, 0, 0, 0],
                     vec![0, 1, 0, 0],
                     vec![0, 0, 1, 0],
                     vec![0, 0, 0, 1]] {
        let rec = positivise(&*fft2_rec(&*example, omega, prime), prime);
        let rec_stride = positivise(&*fft2_stride(&*example, omega, prime), prime);
        let in_place = positivise(&*fft2_in_place(&*example, omega, prime), prime);
        assert_eq!(rec, rec_stride);
        assert_eq!(rec, in_place);
    }
    let prime = 433;
    let omega = 354;
    for example in &[vec![0, 1, 0, 1, 0, 1, 0, 1],
                     vec![0, 1, 0, 0, 0, 0, 0, 0],
                     vec![0, 0, 1, 1, 0, 0, 1, 1],
                     vec![1, 0, 0, 0, 0, 0, 0, 0],
                     vec![1, 1, 1, 1, 1, 1, 1, 1],
                     vec![0, 0, 1, 0, 0, 0, 0, 0]] {
        let rec = positivise(&*fft2_rec(&*example, omega, prime), prime);
        let rec_stride = positivise(&*fft2_stride(&*example, omega, prime), prime);
        let in_place = positivise(&*fft2_in_place(&*example, omega, prime), prime);
        assert_eq!(rec, rec_stride);
        assert_eq!(rec, in_place);
    }
}

#[test]
fn test_fft2() {
    // field is Z_433 in which 354 is an 8th root of unity
    let prime = 433;
    let omega = 354;

    let a_coef = vec![1, 2, 3, 4, 5, 6, 7, 8];
    let a_point = positivise(&*fft2(&a_coef, omega, prime), prime);
    assert_eq!(a_point, vec![36, 303, 146, 3, 429, 422, 279, 122]);
}

#[test]
fn test_fft2_inverse() {
    // field is Z_433 in which 354 is an 8th root of unity
    let prime = 433;
    let omega = 354;

    let a_point = vec![36, -130, -287, 3, -4, 422, 279, -311];
    let a_coef = fft2_inverse(&a_point, omega, prime);
    assert_eq!(positivise(&a_coef, prime), vec![1, 2, 3, 4, 5, 6, 7, 8])
}

#[test]
fn test_fft2_big() {
    let prime = 5038849;
    let omega = 4318906;

    let a_coef: Vec<i64> = (0..256).collect();
    let a_point = fft2(&a_coef, omega, prime);
    let a_coef_back = fft2_inverse(&a_point, omega, prime);

    assert_eq!(positivise(&*a_coef_back, prime), a_coef);
}

#[test]
fn test_trigits_len() {
    assert_eq!(trigits_len(0), 1);
    assert_eq!(trigits_len(1), 1);
    assert_eq!(trigits_len(2), 1);
    assert_eq!(trigits_len(3), 2);
    assert_eq!(trigits_len(8), 2);
    assert_eq!(trigits_len(9), 3);
}

/// Compute recursively the FFT of `a_coef` in the *Zp* field defined by `prime`.
///
/// `omega` must be a principal root of unity. `omega` degree must be equal
/// to the `a_coef` length, and must be a power of 3.
/// The result will contains the same number of elements.
pub fn fft3_rec(a_coef: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    if a_coef.len() == 1 {
        a_coef.to_vec()
    } else {
        // split A(x) into B(x), C(x), and D(x): A(x) = B(x^3) + x C(x^3) + x^2 D(x^3)
        // TODO avoid copying
        let b_coef: Vec<i64> = a_coef.iter()
            .enumerate()
            .filter_map(|(x, &i)| if x % 3 == 0 { Some(i) } else { None })
            .collect();
        let c_coef: Vec<i64> = a_coef.iter()
            .enumerate()
            .filter_map(|(x, &i)| if x % 3 == 1 { Some(i) } else { None })
            .collect();
        let d_coef: Vec<i64> = a_coef.iter()
            .enumerate()
            .filter_map(|(x, &i)| if x % 3 == 2 { Some(i) } else { None })
            .collect();

        // recurse
        let omega_cubed = mod_pow(omega, 3, prime);
        let b_point = fft3_rec(&b_coef, omega_cubed, prime);
        let c_point = fft3_rec(&c_coef, omega_cubed, prime);
        let d_point = fft3_rec(&d_coef, omega_cubed, prime);

        // combine
        let len = a_coef.len();
        let third_len = len / 3;
        let mut a_point = vec![0; len];  // TODO trick: unsafe { Vec.set_len() }
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
        a_point
    }
}

/// Inverse FFT for `fft3`.
pub fn fft3_inverse(a_point: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    let omega_inv = mod_inverse(omega, prime);
    let len = a_point.len();
    let len_inv = mod_inverse(len as i64, prime);
    let scaled_a_coef = fft3(a_point, omega_inv, prime);
    let a_coef = scaled_a_coef.iter().map(|x| x * len_inv % prime).collect();
    a_coef
}

pub fn trigits_len(n: usize) -> usize {
    let mut result = 1;
    let mut value = 3;
    while value < n + 1 {
        result += 1;
        value *= 3;
    }
    result
}

pub fn fft3_in_place_rearrange(data: &mut [u64]) {
    let mut target = 0isize;
    let trigits_len = trigits_len(data.len() - 1);
    let mut trigits: Vec<u8> = ::std::iter::repeat(0).take(trigits_len).collect();
    let powers: Vec<isize> = (0..trigits_len).map(|x| 3isize.pow(x as u32)).rev().collect();
    for pos in 0..data.len() {
        if target as usize > pos {
            data.swap(target as usize, pos)
        }
        for pow in 0..trigits_len {
            if trigits[pow] < 2 {
                trigits[pow] += 1;
                target += powers[pow];
                break;
            } else {
                trigits[pow] = 0;
                target -= 2 * powers[pow];
            }
        }
    }
}

pub fn fft3_in_place_compute(data: &mut [u64], omega: u64, prime: u64) {
    let mut step = 1;
    let max = u64::max_value() / prime / 3;
    let big_omega = mod_pow(omega as i64, (data.len() as u32 / 3), prime as i64) as u64;
    let big_omega_sq = (big_omega * big_omega) % prime;
    while step < data.len() {
        let jump = 3 * step;
        let factor_stride =
            mod_pow(omega as i64, (data.len() / step / 3) as u32, prime as i64) as u64;
        let mut factor = 1;
        for group in 0..step {
            let factor_sq = (factor * factor) % prime;
            let mut pair = group;
            while pair < data.len() {
                let (x, mut y, mut z) =
                    (data[pair], data[pair + step] * factor, data[pair + 2 * step] * factor_sq);

                maybe_reduce!(y, prime, max);
                maybe_reduce!(z, prime, max);

                data[pair] = x + y + z;
                data[pair + step] = x + big_omega * y + big_omega_sq * z;
                data[pair + 2 * step] = x + big_omega_sq * y + big_omega * z;

                maybe_reduce!(data[pair], prime, max);
                maybe_reduce!(data[pair + step], prime, max);
                maybe_reduce!(data[pair + 2 * step], prime, max);

                pair += jump;
            }
            factor = (factor * factor_stride) % prime;
        }
        step = jump;
    }
}

pub fn fft3_in_place(a_coef: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    let mut data: Vec<u64> =
        a_coef.iter().map(|i| (if *i > 0 { *i } else { i + prime }) as u64).collect();
    fft3_in_place_rearrange(&mut *data);
    fft3_in_place_compute(&mut *data, omega as u64, prime as u64);
    data.iter().map(|coef| (coef % prime as u64) as i64).collect()
}

#[test]
fn test_fft3_in_place_rearrange() {
    let mut input = [0, 1, 2, 3, 4, 5, 6, 7, 8];
    fft3_in_place_rearrange(&mut input);
    assert_eq!(input, [0, 3, 6, 1, 4, 7, 2, 5, 8]);
}

#[test]
fn test_fft3() {
    // field is Z_433 in which 150 is an 9th root of unity
    let prime = 433;
    let omega = 150;

    let a_coef = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
    let a_point = positivise(&*fft3_in_place(&a_coef, omega, prime), prime);
    assert_eq!(a_point, vec![45, 404, 407, 266, 377, 47, 158, 17, 20])
}


#[test]
fn test_fft3_variants() {
    let prime = 433;
    let omega = 198;
    for example in &[vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1], vec![1, 1, 1]] {
        let a_point_ref = positivise(&*fft3(&*example, omega, prime), prime);
        let a_in_place = positivise(&*fft3_in_place(&*example, omega, prime), prime);
        assert_eq!(a_in_place, a_point_ref);
    }

    let prime = 433;
    let omega = 27;
    for example in &[vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
                     vec![1, 0, 0, 0, 0, 0, 0, 0, 0],
                     vec![1, 0, 0, 0, 0, 0, 0, 0, 0],
                     vec![0, 0, 0, 1, 0, 0, 0, 0, 0],
                     vec![1, 0, 0, 1, 0, 0, 1, 0, 0],
                     vec![1, 1, 1, 1, 1, 1, 1, 1, 1],
                     vec![1, 2, 3, 4, 5, 6, 7, 8, 9]] {
        let a_point_ref = positivise(&*fft3(&*example, omega, prime), prime);
        let a_in_place = positivise(&*fft3_in_place(&*example, omega, prime), prime);
        assert_eq!(a_in_place, a_point_ref);
    }

    let prime = 5038849;
    let omega = 1814687;
    let example: Vec<i64> = (0..19683).collect();
    let a_point_ref = positivise(&*fft3(&*example, omega, prime), prime);
    let a_in_place = positivise(&*fft3_in_place(&*example, omega, prime), prime);
    assert_eq!(a_in_place, a_point_ref);
}

#[test]
fn test_fft3_inverse() {
    // field is Z_433 in which 150 is an 9th root of unity
    let prime = 433;
    let omega = 150;

    let a_point = vec![45, 404, 407, 266, 377, 47, 158, 17, 20];
    let a_coef = fft3_inverse(&a_point, omega, prime);
    assert_eq!(a_coef, vec![1, 2, 3, 4, 5, 6, 7, 8, 9])
}

#[test]
fn test_fft3_big() {
    let prime = 5038849;
    let omega = 1814687;

    let a_coef: Vec<i64> = (0..19683).collect();
    let a_point = fft3(&a_coef, omega, prime);
    let a_coef_back = fft3_inverse(&a_point, omega, prime);

    assert_eq!(positivise(&*a_coef_back, prime), a_coef);
}
