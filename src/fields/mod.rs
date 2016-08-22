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
    let big_omega = zp.qpow(omega, data.len() as u32 / 2);
    while 1usize << depth < data.len() {
        let step = 1usize << depth;
        let jump = 2 * step;
        let factor_stride = zp.qpow(omega, (data.len() / step / 2) as u32);
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

pub fn fft2_in_place<F: Field>(zp: &F, a_coef: &[F::U], omega: F::U) -> Vec<F::U> {
    let mut data: Vec<F::U> = a_coef.to_vec();
    fft2_in_place_rearrange(zp, &mut *data);
    fft2_in_place_compute(zp, &mut *data, omega);
    data
}

pub fn fft2<F: Field>(zp: &F, a_coef: &[F::U], omega: F::U) -> Vec<F::U> {
    fft2_in_place(zp, a_coef, omega)
}

/// Inverse FFT for `fft2`.
pub fn fft2_inverse<F: Field>(zp: &F, a_point: &[F::U], omega: F::U) -> Vec<F::U> {
    let omega_inv = zp.inv(omega);
    let len = a_point.len();
    let len_inv = zp.inv(zp.from(len as u64));
    let scaled_a_coef = fft2_in_place(zp, a_point, omega_inv);
    let a_coef = scaled_a_coef.iter().map(|&x| zp.mul(x, len_inv)).collect();
    a_coef
}


macro_rules! all_fields_test {
    ($field:ty) => {
        #[test] fn test_qpow() { super::test_qpow::<$field>(); }
        #[test ]fn test_fft2() { super::test_fft2::<$field>(); }
        #[test ]fn test_fft2_inverse() { super::test_fft2_inverse::<$field>(); }
        #[test ]fn test_fft2_big() { super::test_fft2_big::<$field>(); }
    }
}

pub mod naive_zp;
pub mod lazy_zp;

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
