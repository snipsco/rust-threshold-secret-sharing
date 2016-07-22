
pub type Secret = u64;
pub type Share = (u64, u64);

extern crate rand;

mod numtheory;
mod shamir;
mod packed;

pub use shamir::*;
pub use packed::*;
