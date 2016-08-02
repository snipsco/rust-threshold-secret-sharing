# Threshold Secret Sharing

Pure Rust library for [secret sharing](https://en.wikipedia.org/wiki/Secret_sharing),
offering efficient share generation and reconstruction for both traditional Shamir sharing and packet sharing.


# Installation

## Cargo
```toml
[dependencies]
threshold_secret_sharing = "0.1"
```

## Manually run from source code
```bash
git clone https://github.com/snipsco/rust-threshold-secret-sharing
cd rust-threshold-secret-sharing/
cargo build
```

# Examples

## Basic use
```rust
extern crate threshold_secret_sharing as tss;
```
TODO

## Shamir sharing
TODO

## Packed sharing
In this example we're going to pack 3 secrets into each share, using a reconstruction threshold of 4 out of a total of 8 shares.
```rust
let ref pss = tss::PSS_4_8_3;
let secrets = vec![1, 2, 3];
let shares = pss.share(&secrets);
```

## Homomorphic properties
In this example we're going to pack 3 secrets into each share, using a reconstruction threshold of 4 out of a total of 8 shares.
```rust
let ref pss = tss::PSS_4_8_3;

let secrets_1 = vec![1, 2, 3];
let shares_1 = pss.share(&secrets_1);

let secrets_2 = vec![4, 5, 6];
let shares_2 = pss.share(&secrets_2);

// sum the shares pointwise
let shares_sum: Vec<_> = shares_1.iter().zip(&shares_2)
  .map(|(a, b)| (a + b) % pss.prime)
  .collect();

let indices_sum: Vec<usize> = (0..shares_sum.len()).collect();
let secrets_sum = pss.reconstruct(indices_sum, shares_sum);
assert_eq!(secrets_sum, vec![5, 7, 9]);
```
