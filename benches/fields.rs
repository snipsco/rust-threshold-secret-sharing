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
fn bench_ref_fft3(b: &mut test::Bencher) {
    let prime = 5038849;
    let omega = 1814687;

    let a_coef: Vec<i64> = (0..19683).collect();
    assert_eq!(numtheory::mod_pow(omega, a_coef.len() as u32, prime), 1);

    let _a_point = b.iter(|| numtheory::fft3_in_place(&a_coef, omega, prime));
}

macro_rules! all_fields_test {
    ($field:ty) => {
        #[bench] fn test_fft2(b:&mut ::test::Bencher) { super::bench_fft2::<$field>(b); }
        #[bench] fn test_fft3(b:&mut ::test::Bencher) { super::bench_fft3::<$field>(b); }
    }
}

mod naive_zp {
    all_fields_test!(::tss::fields::naive_zp::ZprimeField64);
}

mod lazy_zp {
    all_fields_test!(::tss::fields::lazy_zp::ZprimeField64);
}

fn bench_fft2<F: ZpField>(b: &mut test::Bencher) {
    let zp = F::new(5038849);
    let omega = zp.from(4318906);
    assert_eq!(zp.back(zp.qpow(omega, 256)), 1);

    let a_coef: Vec<_> = (0..256).map(|i| zp.from(i)).collect();
    let _a_point = b.iter(|| tss::fields::fft2(&zp, &a_coef, omega));
}

fn bench_fft3<F: ZpField>(b: &mut test::Bencher) {
    let zp = F::new(5038849);
    let omega = zp.from(1814687);
    assert_eq!(zp.back(zp.qpow(omega, 19683)), 1);

    let a_coef: Vec<_> = (0..19683).map(|a| zp.from(a)).collect();
    let _a_point = b.iter(|| tss::fields::fft3(&zp, &a_coef, omega));
}
