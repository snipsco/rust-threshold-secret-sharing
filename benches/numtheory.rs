#[macro_use]
extern crate bencher;
extern crate threshold_secret_sharing as tss;

use bencher::Bencher;

use tss::numtheory;

fn fft2_div_rec(b: &mut Bencher) {
    let prime = 5038849;
    let omega = 4318906;
    assert_eq!(numtheory::mod_pow(omega, 256, prime), 1);

    let a_coef: Vec<i64> = (0..256).collect();
    let _a_point = b.iter(|| tss::ref_ffts::fft2_rec(&a_coef, omega, prime));
}

fn fft2_div_rec_stride(b: &mut Bencher) {
    let prime = 5038849;
    let omega = 4318906;
    assert_eq!(numtheory::mod_pow(omega, 256, prime), 1);

    let a_coef: Vec<i64> = (0..256).collect();
    let _a_point = b.iter(|| tss::ref_ffts::fft2_stride(&a_coef, omega, prime));
}

fn fft2_div_in_place_eager(b: &mut Bencher) {
    let prime = 5038849;
    let omega = 4318906;
    assert_eq!(numtheory::mod_pow(omega, 256, prime), 1);

    let a_coef: Vec<i64> = (0..256).collect();
    let _a_point = b.iter(|| tss::ref_ffts::fft2_in_place_eager(&a_coef, omega, prime));
}

fn fft2_div_in_place(b: &mut Bencher) {
    let prime = 5038849;
    let omega = 4318906;
    assert_eq!(numtheory::mod_pow(omega, 256, prime), 1);

    let a_coef: Vec<i64> = (0..256).collect();
    let _a_point = b.iter(|| tss::ref_ffts::fft2_in_place(&a_coef, omega, prime));
}

fn fft3_div_rec(b: &mut Bencher) {
    let prime = 5038849;
    let omega = 1814687;

    let a_coef: Vec<i64> = (0..19683).collect();
    assert_eq!(numtheory::mod_pow(omega, a_coef.len() as u32, prime), 1);

    let _a_point = b.iter(|| tss::ref_ffts::fft3_rec(&a_coef, omega, prime));
}

fn fft3_div_in_place(b: &mut Bencher) {
    let prime = 5038849;
    let omega = 1814687;

    let a_coef: Vec<i64> = (0..19683).collect();
    assert_eq!(numtheory::mod_pow(omega, a_coef.len() as u32, prime), 1);

    let _a_point = b.iter(|| tss::ref_ffts::fft3_in_place(&a_coef, omega, prime));
}

fn fft3_div_in_place_eager(b: &mut Bencher) {
    let prime = 5038849;
    let omega = 1814687;

    let a_coef: Vec<i64> = (0..19683).collect();
    assert_eq!(numtheory::mod_pow(omega, a_coef.len() as u32, prime), 1);

    let _a_point = b.iter(|| tss::ref_ffts::fft3_in_place_eager(&a_coef, omega, prime));
}


benchmark_group!(numtheory,
                 fft2_div_rec,
                 fft2_div_rec_stride,
                 fft2_div_in_place_eager,
                 fft2_div_in_place,
                 fft3_div_rec,
                 fft3_div_in_place,
                 fft3_div_in_place_eager
                 );
benchmark_main!(numtheory);
