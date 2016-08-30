// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

//! # Threshold Secret Sharing
//! Pure-Rust library for [secret sharing](https://en.wikipedia.org/wiki/Secret_sharing),
//! offering efficient share generation and reconstruction for both
//! traditional Shamir sharing and packet sharing. For now, secrets and shares
//! are fixed as prime field elements represented by `i64` values.

extern crate rand;

pub mod numtheory;

pub mod shamir;
pub mod packed;
pub mod fields;
