extern crate threshold_secret_sharing as tss;

fn main() {

    let ref tss = tss::shamir::ShamirSecretSharing {
        threshold: 5,
        parts: 20,
        prime: 41  // any large enough prime will do
    };

    let secret = 5;
    let all_shares = tss.share(secret);

    let reconstruct_share_count = 10;
    assert!(reconstruct_share_count >= tss.threshold);

    let indices: Vec<usize> = (0..reconstruct_share_count).collect();
    let shares: &[i64] = &all_shares[0..reconstruct_share_count];
    let recovered_secret = tss.reconstruct(&indices, shares);

    println!("The recovered secret is {}", recovered_secret);
    assert_eq!(recovered_secret, secret);

}
