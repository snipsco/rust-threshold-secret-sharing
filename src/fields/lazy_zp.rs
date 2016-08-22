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
    pub prime: u64,
}

impl ZprimeField64 {
    pub fn new(prime: u64) -> ZprimeField64 {
        ZprimeField64 { prime: prime }
    }
}

macro_rules! maybe_reduce {
    ($target:expr, $prime:ident, $max:ident) => {
        if $target > $max {
            $target %= $prime;
        }
    }
}

impl Field for ZprimeField64 {
    type U = ZprimeU64;

    fn from(&self, a: u64) -> Self::U {
        if a > u32::max_value() as u64 {
            ZprimeU64(a % self.prime)
        } else {
            ZprimeU64(a)
        }
    }

    fn back(&self, a: Self::U) -> u64 {
        a.0 % self.prime
    }

    fn optimize(&self, a: Self::U) -> Self::U {
        ZprimeU64(a.0 % self.prime)
    }

    fn add(&self, a: Self::U, b: Self::U) -> Self::U {
        self.from(a.0 + b.0)
    }

    fn sub(&self, a: Self::U, b: Self::U) -> Self::U {
        if b.0 > a.0 {
            self.from(self.prime + a.0 - (b.0 % self.prime))
        } else {
            self.from(a.0 - b.0)
        }
    }

    fn mul(&self, a: Self::U, b: Self::U) -> Self::U {
        self.from(a.0 * b.0)
    }

    fn inv(&self, a: Self::U) -> Self::U {
        let tmp = ::numtheory::mod_inverse((a.0 % self.prime) as i64, self.prime as i64);
        self.from((tmp % self.prime as i64) as u64)
    }
}

impl ZpField for ZprimeField64 {
    fn new(prime: u64) -> ZprimeField64 {
        ZprimeField64::new(prime)
    }
}

all_fields_test!(ZprimeField64);
