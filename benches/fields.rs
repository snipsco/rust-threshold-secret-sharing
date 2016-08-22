#![feature(test)]

extern crate test;
extern crate threshold_secret_sharing as tss;

use tss::numtheory;
use tss::fields::ZpField;

#[bench]
fn bench_ref_fft2(b: &mut test::Bencher) {
    let prime = 5038849;
    let omega = 4318906;
    assert_eq!(numtheory::mod_pow(omega, 256, prime), 1);

    let a_coef: Vec<i64> = (0..256).collect();
    let _a_point = b.iter(|| numtheory::fft2_in_place(&a_coef, omega, prime));
}

#[bench]
fn bench_ref_fft2_naive_zp(b: &mut test::Bencher) {
    bench_ref_fft2_zp::<tss::fields::naive_zp::ZprimeField64>(b)
}

#[bench]
fn bench_ref_fft2_lazy_zp(b: &mut test::Bencher) {
    bench_ref_fft2_zp::<tss::fields::lazy_zp::ZprimeField64>(b)
}

fn bench_ref_fft2_zp<F: ZpField>(b: &mut test::Bencher) {
    let zp = F::new(5038849);
    let omega = zp.from(4318906);
    assert_eq!(zp.back(zp.qpow(omega, 256)), 1);

    let a_coef:Vec<_> = (0..256).map(|i| zp.from(i)).collect();
    let _a_point = b.iter(|| tss::fields::fft2(&zp, &a_coef, omega));
}
