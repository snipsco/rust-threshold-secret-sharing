
use super::{Secret, Share};
use single::*;


#[derive(Debug)]
pub struct MultiThresholdSecretSharing {
    tss: SingleThresholdSecretSharing
}


impl MultiThresholdSecretSharing {

    pub fn new(parts: usize, threshold: usize, prime: u64) -> MultiThresholdSecretSharing {
        let tss = SingleThresholdSecretSharing::new(parts, threshold, prime);
        MultiThresholdSecretSharing {
            tss: tss
        }
    }

    fn share(&self, secrets: &Vec<Secret>) -> Vec<Vec<Share>> {
        secrets.iter().map(|&secret| {
            self.tss.share(secret)
        }).collect()
    }

    fn reconstruct(&self, shares: &[Vec<Share>]) -> Vec<Secret> {
        shares.iter().map(|shares| {
            self.tss.reconstruct(&*shares)
        }).collect()
    }

}


#[test]
fn multi() {
    let tss = MultiThresholdSecretSharing::new(20, 5, 41);

    let secrets = vec![1, 2, 3];
    let shares = tss.share(&secrets);
    let recon_secrets = tss.reconstruct(&*shares);
    assert_eq!(secrets, recon_secrets);
}
