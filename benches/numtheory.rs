#![feature(test)]

extern crate test;
extern crate threshold_secret_sharing as tss;

use tss::numtheory;

#[bench]
fn bench_fft2_ref(b: &mut test::Bencher) {
    let prime = 5038849;
    let omega = 4318906;
    assert_eq!(numtheory::mod_pow(omega, 256, prime), 1);

    let a_coef: Vec<i64> = (0..256).collect();
    let _a_point = b.iter(|| numtheory::fft2_ref(&a_coef, omega, prime));
}

#[bench]
fn bench_fft2_stride(b: &mut test::Bencher) {
    let prime = 5038849;
    let omega = 4318906;
    assert_eq!(numtheory::mod_pow(omega, 256, prime), 1);

    let a_coef: Vec<i64> = (0..256).collect();
    let _a_point = b.iter(|| numtheory::fft2_stride(&a_coef, omega, prime));
}

#[bench]
fn bench_fft2_in_place(b: &mut test::Bencher) {
    let prime = 5038849;
    let omega = 4318906;
    assert_eq!(numtheory::mod_pow(omega, 256, prime), 1);

    let a_coef: Vec<i64> = (0..256).collect();
    let _a_point = b.iter(|| numtheory::fft2_in_place(&a_coef, omega, prime));
}

#[bench]
fn bench_fft3_ref(b: &mut test::Bencher) {
    let prime = 5038849;
    let omega = 1814687;

    let a_coef: Vec<i64> = (0..19683).collect();
    assert_eq!(numtheory::mod_pow(omega, a_coef.len() as u32, prime), 1);

    let _a_point = b.iter(|| numtheory::fft3_ref(&a_coef, omega, prime));
}

#[bench]
fn bench_fft3_in_place(b: &mut test::Bencher) {
    let prime = 5038849;
    let omega = 1814687;

    let a_coef: Vec<i64> = (0..19683).collect();
    assert_eq!(numtheory::mod_pow(omega, a_coef.len() as u32, prime), 1);

    let _a_point = b.iter(|| numtheory::fft3_in_place(&a_coef, omega, prime));
}
