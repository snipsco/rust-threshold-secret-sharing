// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

pub trait Field {
    type U: Copy;

    fn from(&self, a: u64) -> Self::U;
    fn back(&self, a: Self::U) -> u64;

    fn zero(&self) -> Self::U {
        self.from(0)
    }

    fn one(&self) -> Self::U {
        self.from(1)
    }

    fn optimize(&self, a: Self::U) -> Self::U {
        a
    }

    fn add(&self, a: Self::U, b: Self::U) -> Self::U;
    fn sub(&self, a: Self::U, b: Self::U) -> Self::U;
    fn mul(&self, a: Self::U, b: Self::U) -> Self::U;

    fn inv(&self, a: Self::U) -> Self::U;

    fn qpow(&self, mut x: Self::U, mut e: u32) -> Self::U {
        if e == 0 {
            return self.one();
        }
        let mut acc = self.one();
        while e > 0 {
            if e % 2 == 0 {
                // even
                // no-op
            } else {
                // odd
                acc = self.mul(acc, x);
            }
            x = self.mul(x, x);  // waste one of these by having it here but code is simpler (tiny bit)
            e = e >> 1;
        }
        acc
    }

    fn from_v(&self, values: &[u64]) -> Vec<Self::U> {
        values.iter().map(|&a| self.from(a)).collect()
    }

    fn back_v(&self, values: &[Self::U]) -> Vec<u64> {
        values.iter().map(|&a| self.back(a)).collect()
    }
}

pub trait ZpField: Field {
    fn new(prime: u64) -> Self;
}

pub fn fft2_in_place_rearrange<F: Field>(_zp: &F, data: &mut [F::U]) {
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

pub fn fft2_in_place_compute<F: Field>(zp: &F, data: &mut [F::U], omega: F::U) {
    let mut depth = 0usize;
    let big_omega = zp.optimize(zp.qpow(omega, data.len() as u32 / 2));
    while 1usize << depth < data.len() {
        let step = 1usize << depth;
        let jump = 2 * step;
        let factor_stride = zp.optimize(zp.qpow(omega, (data.len() / step / 2) as u32));
        let mut factor = zp.one();
        for group in 0usize..step {
            let mut pair = group;
            while pair < data.len() {
                let (x, y) = (data[pair], zp.mul(data[pair + step], factor));

                data[pair] = zp.add(x, y);
                data[pair + step] = zp.add(x, zp.mul(big_omega, y));

                pair += jump;
            }
            factor = zp.optimize(zp.mul(factor, factor_stride));
        }
        depth += 1;
    }
}

pub fn fft2<F: Field>(zp: &F, a_coef: &[F::U], omega: F::U) -> Vec<F::U> {
    let mut data: Vec<F::U> = a_coef.to_vec();
    fft2_in_place_rearrange(zp, &mut *data);
    fft2_in_place_compute(zp, &mut *data, omega);
    data
}

pub fn fft2_inverse<F: Field>(zp: &F, a_point: &[F::U], omega: F::U) -> Vec<F::U> {
    let omega_inv = zp.inv(omega);
    let len = a_point.len();
    let len_inv = zp.inv(zp.from(len as u64));
    let scaled_a_coef = fft2(zp, a_point, omega_inv);
    let a_coef = scaled_a_coef.iter().map(|&x| zp.mul(x, len_inv)).collect();
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

pub fn fft3_in_place_rearrange<F:Field>(_zp:&F, data: &mut [F::U]) {
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

pub fn fft3_in_place_compute<F:Field>(zp: &F, data: &mut [F::U], omega: F::U) {
    let mut step = 1;
    let big_omega = zp.optimize(zp.qpow(omega, (data.len() as u32 / 3)));
    let big_omega_sq = zp.optimize(zp.mul(big_omega, big_omega));
    while step < data.len() {
        let jump = 3 * step;
        let factor_stride = zp.optimize(zp.qpow(omega, (data.len() / step / 3) as u32));
        let mut factor = zp.one();
        for group in 0usize..step {
            let factor_sq = zp.mul(factor, factor);
            let mut pair = group;
            while pair < data.len() {
                let (x, y, z) =
                    (data[pair], zp.mul(data[pair + step], factor), zp.mul(data[pair + 2 * step],factor_sq));

                data[pair] = zp.add(zp.add(x, y), z);
                data[pair + step] = zp.add(zp.add(x, zp.mul(big_omega,y)), zp.mul(big_omega_sq,z));
                data[pair + 2 * step] = zp.add(zp.add(x, zp.mul(big_omega_sq, y)), zp.mul(big_omega, z));

                pair += jump;
            }
            factor = zp.optimize(zp.mul(factor, factor_stride));
        }
        step = jump;
    }
}

pub fn fft3<F:Field>(zp: &F, a_coef: &[F::U], omega: F::U) -> Vec<F::U> {
    let mut data: Vec<F::U> = a_coef.to_vec();
    fft3_in_place_rearrange(zp, &mut *data);
    fft3_in_place_compute(zp, &mut *data, omega);
    data
}

pub fn fft3_inverse<F:Field>(zp: &F, a_point: &[F::U], omega: F::U) -> Vec<F::U> {
    let omega_inv = zp.inv(omega);
    let len_inv = zp.inv(zp.from(a_point.len() as u64));
    let scaled_a_coef = fft3(zp, a_point, omega_inv);
    let a_coef = scaled_a_coef.iter().map(|&x| zp.mul(x, len_inv)).collect();
    a_coef
}

macro_rules! all_fields_test {
    ($field:ty) => {
        #[test] fn test_convert() { super::test_convert::<$field>(); }
        #[test] fn test_add() { super::test_add::<$field>(); }
        #[test] fn test_sub() { super::test_sub::<$field>(); }
        #[test] fn test_mul() { super::test_mul::<$field>(); }
        #[test] fn test_qpow() { super::test_qpow::<$field>(); }
        #[test ]fn test_fft2() { super::test_fft2::<$field>(); }
        #[test ]fn test_fft2_inverse() { super::test_fft2_inverse::<$field>(); }
        #[test ]fn test_fft2_big() { super::test_fft2_big::<$field>(); }
        #[test ]fn test_fft3() { super::test_fft3::<$field>(); }
        #[test ]fn test_fft3_inverse() { super::test_fft3_inverse::<$field>(); }
        #[test ]fn test_fft3_big() { super::test_fft3_big::<$field>(); }
    }
}

pub mod naive_zp;
pub mod lazy_zp;
pub mod montgomery;

#[cfg(test)]
fn test_convert<F: ZpField>() {
    let zp = F::new(17);
    for i in 0u64..20 {
        assert_eq!(zp.back(zp.from(i)), i % 17);
    }
}

#[cfg(test)]
fn test_add<F: ZpField>() {
    let zp = F::new(17);
    assert_eq!(zp.back(zp.add(zp.from(8), zp.from(2))), 10);
    assert_eq!(zp.back(zp.add(zp.from(8), zp.from(13))), 4);
}

#[cfg(test)]
fn test_sub<F: ZpField>() {
    let zp = F::new(17);
    assert_eq!(zp.back(zp.sub(zp.from(8), zp.from(2))), 6);
    assert_eq!(zp.back(zp.sub(zp.from(8), zp.from(13))), 12);
}

#[cfg(test)]
fn test_mul<F: ZpField>() {
    let zp = F::new(17);
    assert_eq!(zp.back(zp.mul(zp.from(8), zp.from(2))), 16);
    assert_eq!(zp.back(zp.mul(zp.from(8), zp.from(5))), 6);
}

#[cfg(test)]
fn test_qpow<F: ZpField>() {
    let zp = F::new(17);
    assert_eq!(zp.back(zp.qpow(zp.from(2), 0)), 1);
    assert_eq!(zp.back(zp.qpow(zp.from(2), 3)), 8);
    assert_eq!(zp.back(zp.qpow(zp.from(2), 6)), 13);
}

#[cfg(test)]
fn test_fft2<F: ZpField>() {
    // field is Z_433 in which 354 is an 8th root of unity
    let zp = F::new(433);
    let omega = zp.from(354);

    let a_coef = zp.from_v(&[1, 2, 3, 4, 5, 6, 7, 8]);
    let a_point = fft2(&zp, &*a_coef, omega);
    assert_eq!(zp.back_v(&a_point), [36, 303, 146, 3, 429, 422, 279, 122]);
}

#[cfg(test)]
fn test_fft2_inverse<F: ZpField>() {
    // field is Z_433 in which 354 is an 8th root of unity
    let zp = F::new(433);
    let omega = zp.from(354);

    let a_point = zp.from_v(&[36, 303, 146, 3, 429, 422, 279, 122]);
    let a_coef = fft2_inverse(&zp, &a_point, omega);
    assert_eq!(zp.back_v(&a_coef), [1, 2, 3, 4, 5, 6, 7, 8])
}

#[cfg(test)]
fn test_fft2_big<F: ZpField>() {
    let zp = F::new(5038849);
    let omega = zp.from(4318906);

    let a_coef: Vec<_> = (0..256).map(|a| zp.from(a)).collect();
    let a_point = fft2(&zp, &a_coef, omega);
    let a_coef_back = fft2_inverse(&zp, &a_point, omega);

    assert_eq!(zp.back_v(&a_coef_back), zp.back_v(&a_coef));
}

#[cfg(test)]
fn test_fft3<F:ZpField>() {
    // field is Z_433 in which 150 is an 9th root of unity
    let zp = F::new(433);
    let omega = zp.from(150);

    let a_coef = zp.from_v(&[1, 2, 3, 4, 5, 6, 7, 8, 9]);
    let a_point = fft3(&zp, &a_coef, omega);
    assert_eq!(zp.back_v(&a_point), [45, 404, 407, 266, 377, 47, 158, 17, 20]);
}

#[cfg(test)]
fn test_fft3_inverse<F:ZpField>() {
    // field is Z_433 in which 150 is an 9th root of unity
    let zp = F::new(433);
    let omega = zp.from(150);

    let a_point = zp.from_v(&[45, 404, 407, 266, 377, 47, 158, 17, 20]);
    let a_coef = fft3_inverse(&zp, &a_point, omega);
    assert_eq!(zp.back_v(&a_coef), [1, 2, 3, 4, 5, 6, 7, 8, 9])
}

#[cfg(test)]
fn test_fft3_big<F:ZpField>() {
    let zp = F::new(5038849);
    let omega = zp.from(1814687);

    let a_coef: Vec<_> = (0..19683).map(|a| zp.from(a)).collect();
    let a_point = fft3(&zp, &a_coef, omega);
    let a_coef_back = fft3_inverse(&zp, &a_point, omega);

    assert_eq!(zp.back_v(&a_coef_back), zp.back_v(&a_coef));
}
