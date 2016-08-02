extern crate rand;

pub type Secret = u64;
pub type Share = (u64, u64);

pub mod numtheory;

pub mod shamir;
pub mod packed;
