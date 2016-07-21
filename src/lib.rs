
pub type Secret = u64;
pub type Share = (u64, u64);

extern crate rand;

mod numtheory;
mod single;
mod multi;

pub use single::*;
pub use multi::*;
