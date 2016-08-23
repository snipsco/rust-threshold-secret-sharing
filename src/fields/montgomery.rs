// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

use super::{Field, ZpField};

#[derive(Copy,Clone,Debug)]
pub struct ZprimeU32(pub u32);

pub struct ZprimeField32 {
    pub n: u32, // the prime
    pub n_quote: u32,
    pub r_inv: u32, // r = 2^32
}

impl ZprimeField32 {
    pub fn new(prime: u32) -> ZprimeField32 {
        let r = 1u64 << 32;
        let tmp = ::numtheory::mod_inverse(r as i64, prime as i64);
        let r_inv = if tmp < 0 {
            (tmp + prime as i64) as u32
        } else {
            tmp as u32
        };
        let tmp = ::numtheory::mod_inverse(prime as i64, r as i64);
        let n_quote = if tmp > 0 {
            (r as i64 - tmp) as u32
        } else {
            (r as i64 - tmp) as u32
        };
        ZprimeField32 {
            n: prime,
            r_inv: r_inv,
            n_quote: n_quote,
        }
    }

    fn redc(&self, a: u64) -> ZprimeU32 {
        let m: u64 = (a as u32).wrapping_mul(self.n_quote) as u64;
        let t: u32 = ((a + m * (self.n as u64)) >> 32) as u32;
        ZprimeU32((if t >= (self.n) { t - (self.n) } else { t }))
    }
}

impl Field for ZprimeField32 {
    type U = ZprimeU32;

    fn from(&self, a: u64) -> Self::U {
        ZprimeU32(((a << 32) % self.n as u64) as u32)
    }

    fn back(&self, a: Self::U) -> u64 {
        a.0 as u64 * self.r_inv as u64 % self.n as u64
    }

    fn add(&self, a: Self::U, b: Self::U) -> Self::U {
        let sum = a.0 as u64 + b.0 as u64;
        if sum > self.n as u64 {
            ZprimeU32((sum - self.n as u64) as u32)
        } else {
            ZprimeU32(sum as u32)
        }
    }

    fn sub(&self, a: Self::U, b: Self::U) -> Self::U {
        if a.0 > b.0 {
            ZprimeU32(a.0 - b.0)
        } else {
            ZprimeU32((a.0 as u64 + self.n as u64 - b.0 as u64) as u32)
        }
    }

    fn mul(&self, a: Self::U, b: Self::U) -> Self::U {
        self.redc((a.0 as u64).wrapping_mul(b.0 as u64))
    }

    fn inv(&self, a: Self::U) -> Self::U {
        let a = self.back(a);
        self.from(::numtheory::mod_inverse(a as i64, self.n as i64) as u64)
    }
}

impl ZpField for ZprimeField32 {
    fn new(prime: u64) -> ZprimeField32 {
        ZprimeField32::new(prime as u32)
    }
}

all_fields_test!(ZprimeField32);
