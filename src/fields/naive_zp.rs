// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

use super::{Field, ZpField};

#[derive(Copy,Clone,Debug)]
pub struct ZprimeU64(pub i64);

pub struct ZprimeField64(pub i64);

impl Field for ZprimeField64 {
    type U = ZprimeU64;

    fn from(&self, a: u64) -> Self::U {
        ZprimeU64(a as i64 % self.0)
    }

    fn back(&self, a: Self::U) -> u64 {
        a.0 as u64
    }

    fn add(&self, a: Self::U, b: Self::U) -> Self::U {
        ZprimeU64((a.0 + b.0) % self.0)
    }

    fn sub(&self, a: Self::U, b: Self::U) -> Self::U {
        let tmp = a.0 - b.0;
        if tmp > 0 {
            ZprimeU64(tmp)
        } else {
            ZprimeU64(tmp + self.0)
        }
    }

    fn mul(&self, a: Self::U, b: Self::U) -> Self::U {
        ZprimeU64((a.0 * b.0) % self.0)
    }

    fn inv(&self, a: Self::U) -> Self::U {
        let tmp = ::numtheory::mod_inverse((a.0 % self.0) as i64, self.0 as i64);
        self.from((tmp % self.0 as i64) as u64)
    }
}

impl ZpField for ZprimeField64 {
    fn new(prime: u64) -> ZprimeField64 {
        ZprimeField64(prime as i64)
    }
}

all_fields_test!(ZprimeField64);
