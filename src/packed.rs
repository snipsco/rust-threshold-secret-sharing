
use super::{Secret, Share};
use numtheory::*;


#[derive(Debug)]
pub struct PackedSecretSharing {
    threshold: usize,
    parts: usize,
    prime: i64,
    omega_n: i64,
    omega_m: i64
}


impl PackedSecretSharing {

    pub fn new() -> PackedSecretSharing {
        
    }

    // pub fn new(parts: usize, threshold: usize, prime: u64) -> PackedSecretSharing {
    //     let tss = ShamirSecretSharing::new(parts, threshold, prime);
    //     PackedSecretSharing {
    //         tss: tss
    //     }
    // }
    //
    // fn share(&self, secrets: &Vec<Secret>) -> Vec<Vec<Share>> {
    //     secrets.iter().map(|&secret| {
    //         self.tss.share(secret)
    //     }).collect()
    // }
    //
    // fn reconstruct(&self, shares: &[Vec<Share>]) -> Vec<Secret> {
    //     shares.iter().map(|shares| {
    //         self.tss.reconstruct(&*shares)
    //     }).collect()
    // }

}


// #[test]
// fn multi() {
//     let tss = PackedSecretSharing::new(20, 5, 41);
//
//     let secrets = vec![1, 2, 3];
//     let shares = tss.share(&secrets);
//     let recon_secrets = tss.reconstruct(&*shares);
//     assert_eq!(secrets, recon_secrets);
// }
