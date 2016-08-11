#![feature(test)]

extern crate test;
extern crate threshold_secret_sharing as tss;

mod packed {

    use test::Bencher;
    use tss::packed::*;

    #[bench]
    fn bench_large_share(b: &mut Bencher) {
        let ref pss = PSS_155_19682_100;
        let secrets = vec![5 ; pss.secret_count];
        b.iter(|| {
            let _shares = pss.share(&secrets);
        });
    }

    #[bench]
    fn bench_large_reconstruct(b: &mut Bencher) {
        let ref pss = PSS_155_19682_100;
        let secrets = vec![5 ; pss.secret_count];
        let all_shares = pss.share(&secrets);

        // reconstruct using minimum number of shares required
        let indices: Vec<usize> = (0..pss.reconstruct_limit).collect();
        let shares = &all_shares[0..pss.reconstruct_limit];

        b.iter(|| {
            let _recovered_secrets = pss.reconstruct(&indices, &shares);
        });
    }

}
