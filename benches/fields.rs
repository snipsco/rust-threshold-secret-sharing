// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.
#[macro_use]
extern crate bencher;
extern crate threshold_secret_sharing as tss;

use bencher::Bencher;

use tss::numtheory;
use tss::fields::ZpField;

fn ref_fft2(b: &mut Bencher) {
    let prime = 5038849;
    let omega = 4318906;
    assert_eq!(numtheory::mod_pow(omega, 256, prime), 1);

    let a_coef: Vec<i64> = (0..256).collect();
    let _a_point = b.iter(|| numtheory::fft2_in_place(&a_coef, omega, prime));
}

fn ref_fft3(b: &mut Bencher) {
    let prime = 5038849;
    let omega = 1814687;

    let a_coef: Vec<i64> = (0..19683).collect();
    assert_eq!(numtheory::mod_pow(omega, a_coef.len() as u32, prime), 1);

    let _a_point = b.iter(|| numtheory::fft3_in_place(&a_coef, omega, prime));
}

macro_rules! all_fields_test {
    ($group:ident, $field:ty) => {
        mod $group {
            use bencher::Bencher;
            pub fn test_fft2(b:&mut Bencher) { super::bench_fft2::<$field>(b); }
            pub fn test_fft3(b:&mut Bencher) { super::bench_fft3::<$field>(b); }
        }
        benchmark_group!($group, $group::test_fft2, $group::test_fft3);
    }
}

fn bench_fft2<F: ZpField>(b: &mut Bencher) {
    let zp = F::new(5038849);
    let omega = zp.from(4318906);
    assert_eq!(zp.back(zp.qpow(omega, 256)), 1);

    let a_coef: Vec<_> = (0..256).map(|i| zp.from(i)).collect();
    let _a_point = b.iter(|| tss::fields::fft2(&zp, &a_coef, omega));
}

fn bench_fft3<F: ZpField>(b: &mut Bencher) {
    let zp = F::new(5038849);
    let omega = zp.from(1814687);
    assert_eq!(zp.back(zp.qpow(omega, 19683)), 1);

    let a_coef: Vec<_> = (0..19683).map(|a| zp.from(a)).collect();
    let _a_point = b.iter(|| tss::fields::fft3(&zp, &a_coef, omega));
}

all_fields_test!(naive_zp, ::tss::fields::naive_zp::ZprimeField64);
all_fields_test!(lazy_zp, ::tss::fields::lazy_zp::ZprimeField64);
all_fields_test!(montgomery, ::tss::fields::montgomery::ZprimeField32);
benchmark_group!(reference, ref_fft2, ref_fft3);

benchmark_main!(reference, naive_zp, lazy_zp, montgomery);
