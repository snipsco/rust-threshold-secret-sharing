// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

pub mod ffts;

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

macro_rules! all_fields_test {
    ($field:ty) => {
        #[test] fn test_convert() { super::test_convert::<$field>(); }
        #[test] fn test_add() { super::test_add::<$field>(); }
        #[test] fn test_sub() { super::test_sub::<$field>(); }
        #[test] fn test_mul() { super::test_mul::<$field>(); }
        #[test] fn test_qpow() { super::test_qpow::<$field>(); }
        #[test ]fn test_fft2() { super::ffts::test_fft2::<$field>(); }
        #[test ]fn test_fft2_inverse() { super::ffts::test_fft2_inverse::<$field>(); }
        #[test ]fn test_fft2_big() { super::ffts::test_fft2_big::<$field>(); }
        #[test ]fn test_fft3() { super::ffts::test_fft3::<$field>(); }
        #[test ]fn test_fft3_inverse() { super::ffts::test_fft3_inverse::<$field>(); }
        #[test ]fn test_fft3_big() { super::ffts::test_fft3_big::<$field>(); }
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
    assert_eq!(zp.back(zp.sub(zp.from(8), zp.from(13))), (17 + 8 - 13) % 17);
}

#[cfg(test)]
fn test_mul<F: ZpField>() {
    let zp = F::new(17);
    assert_eq!(zp.back(zp.mul(zp.from(8), zp.from(2))), (8 * 2) % 17);
    assert_eq!(zp.back(zp.mul(zp.from(8), zp.from(5))), (8 * 5) % 17);
}

#[cfg(test)]
fn test_qpow<F: ZpField>() {
    let zp = F::new(17);
    assert_eq!(zp.back(zp.qpow(zp.from(2), 0)), 1);
    assert_eq!(zp.back(zp.qpow(zp.from(2), 3)), 8);
    assert_eq!(zp.back(zp.qpow(zp.from(2), 6)), 13);
}

