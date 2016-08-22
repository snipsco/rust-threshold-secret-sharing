// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

use super::{Field, ZpField};

#[derive(Copy,Clone,Debug)]
pub struct ZprimeU64(pub u64);

pub struct ZprimeField64 {
    pub n: u64, // the prime
    pub n_quote: u64,
    pub r_inv: u64, // r = 2^32
}

impl ZprimeField64 {
    pub fn new(prime: u64) -> ZprimeField64 {
        let r = 1u64 << 32;
        let tmp = ::numtheory::mod_inverse(r as i64, prime as i64);
        let r_inv = if tmp < 0 {
            (tmp + prime as i64) as u64
        } else {
            tmp as u64
        };
        let tmp = ::numtheory::mod_inverse(prime as i64, r as i64);
        let n_quote = if tmp > 0 {
            (r as i64 - tmp) as u64
        } else {
            (r as i64-tmp) as u64
        };
        ZprimeField64 {
            n: prime,
            r_inv: r_inv,
            n_quote: n_quote,
        }
    }

    pub fn mul_r(&self, a: u64) -> u64 {
        a << 32
    }

    pub fn mod_r(&self, a: u64) -> u64 {
        a & (0xFFFF_FFFF as u64)
    }

    pub fn div_r(&self, a: u64) -> u64 {
        a >> 32
    }

    fn redc(&self, a: u64) -> ZprimeU64 {
        let m = self.mod_r(self.mod_r(a).wrapping_mul(self.n_quote));
        let t = self.div_r(a + m * self.n);
        ZprimeU64(if t >= self.n { t - self.n } else { t })
    }
}

impl Field for ZprimeField64 {
    type U = ZprimeU64;

    fn from(&self, a: u64) -> Self::U {
        ZprimeU64(self.mul_r(a) % self.n)
    }

    fn back(&self, a: Self::U) -> u64 {
        a.0 * self.r_inv % self.n
    }

    fn add(&self, a: Self::U, b: Self::U) -> Self::U {
        ZprimeU64(self.mod_r(a.0.wrapping_add(b.0)))
    }

    fn sub(&self, a: Self::U, b: Self::U) -> Self::U {
        unimplemented!();
        //        ZprimeU64(a.0.wrapping_sub(b.0) & 0x0000_0000_FFFF_FFFFu64)
    }

    fn mul(&self, a: Self::U, b: Self::U) -> Self::U {
        self.redc(a.0.wrapping_mul(b.0))
    }

    fn inv(&self, a: Self::U) -> Self::U {
        let a = self.back(a);
        self.from(::numtheory::mod_inverse(a as i64, self.n as i64) as u64)
    }
}

impl ZpField for ZprimeField64 {
    fn new(prime: u64) -> ZprimeField64 {
        ZprimeField64::new(prime)
    }
}

all_fields_test!(ZprimeField64);
