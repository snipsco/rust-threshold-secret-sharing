// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

pub fn gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let n = a / b;
        let c = a % b;
        let r = gcd(b, c);
        (r.0, r.2, r.1 - r.2 * n)
    }
}

#[test]
fn test_gcd() {
    assert_eq!(gcd(12, 16), (4, -1, 1));
}


pub fn mod_inverse(k: i64, prime: i64) -> i64 {
    let k2 = k % prime;
    let r = if k2 < 0 {
        -gcd(prime, -k2).2
    } else {
        gcd(prime, k2).2
    };
    (prime + r) % prime
}

#[test]
fn test_mod_inverse() {
    assert_eq!(mod_inverse(3, 7), 5);
}


pub fn mod_pow(mut x: i64, mut e: u32, prime: i64) -> i64 {
    let mut acc = 1;
    while e > 0 {
        if e % 2 == 0 {
            // even
            // no-op
        }
        else {
            // odd
            acc = (acc * x) % prime;
        }
        x = (x * x) % prime;  // waste one of these by having it here but code is simpler (tiny bit)
        e = e >> 1;
    }
    acc
}

#[test]
fn test_mod_pow() {
    assert_eq!(mod_pow(2, 0, 17), 1);
    assert_eq!(mod_pow(2, 3, 17), 8);
    assert_eq!(mod_pow(2, 6, 17), 13);

    assert_eq!(mod_pow(-3, 0, 17), 1);
    assert_eq!(mod_pow(-3, 1, 17), -3);
    assert_eq!(mod_pow(-3, 15, 17), -6);
}


pub fn fft2(a_coef: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    if a_coef.len() == 1 {
        a_coef.to_vec()
    } else {
        // split A(x) into B(x) and C(x): A(x) = B(x^2) + x C(x^2)
        // TODO avoid copying
        let b_coef: Vec<i64> = a_coef.iter().enumerate().filter_map(|(x,&i)| if x % 2 == 0 { Some(i) } else { None } ).collect();
        let c_coef: Vec<i64> = a_coef.iter().enumerate().filter_map(|(x,&i)| if x % 2 == 1 { Some(i) } else { None } ).collect();

        // recurse
        let b_point = fft2(&b_coef, mod_pow(omega, 2, prime), prime);
        let c_point = fft2(&c_coef, mod_pow(omega, 2, prime), prime);

        // combine
        let len = a_coef.len();
        let half_len = len >> 1;
        let mut a_point = vec![0; len];  // TODO trick: unsafe { Vec.set_len() }
        for i in 0..half_len {
            a_point[i]            = (b_point[i] + mod_pow(omega, i as u32, prime) * c_point[i]) % prime;
            a_point[i + half_len] = (b_point[i] - mod_pow(omega, i as u32, prime) * c_point[i]) % prime;
        }

        // return
        a_point
    }
}

pub fn fft2_inverse(a_point: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    let omega_inv = mod_inverse(omega, prime);
    let len = a_point.len();
    let len_inv = mod_inverse(len as i64, prime);
    let scaled_a_coef = fft2(a_point, omega_inv, prime);
    let a_coef = scaled_a_coef.iter().map(|x| { x * len_inv % prime }).collect();
    a_coef
}

#[test]
fn test_fft2() {
    // field is Z_433 in which 354 is an 8th root of unity
    let prime = 433;
    let omega = 354;

    let a_coef = vec![1,2,3,4,5,6,7,8];
    let a_point = fft2(&a_coef, omega, prime);
    assert_eq!(a_point, vec![36, -130, -287, 3, -4, 422, 279, -311])
}

#[test]
fn test_fft2_inverse() {
    // field is Z_433 in which 354 is an 8th root of unity
    let prime = 433;
    let omega = 354;

    let a_point = vec![36, -130, -287, 3, -4, 422, 279, -311];
    let a_coef = fft2_inverse(&a_point, omega, prime);
    assert_eq!(positivise(&a_coef, prime), vec![1,2,3,4,5,6,7,8])
}


pub fn fft3(a_coef: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    if a_coef.len() == 1 {
        a_coef.to_vec()
    } else {
        // split A(x) into B(x), C(x), and D(x): A(x) = B(x^3) + x C(x^3) + x^2 D(x^3)
        // TODO avoid copying
        let b_coef: Vec<i64> = a_coef.iter().enumerate().filter_map(|(x,&i)| if x % 3 == 0 { Some(i) } else { None } ).collect();
        let c_coef: Vec<i64> = a_coef.iter().enumerate().filter_map(|(x,&i)| if x % 3 == 1 { Some(i) } else { None } ).collect();
        let d_coef: Vec<i64> = a_coef.iter().enumerate().filter_map(|(x,&i)| if x % 3 == 2 { Some(i) } else { None } ).collect();

        // recurse
        let omega_cubed = mod_pow(omega, 3, prime);
        let b_point = fft3(&b_coef, omega_cubed, prime);
        let c_point = fft3(&c_coef, omega_cubed, prime);
        let d_point = fft3(&d_coef, omega_cubed, prime);

        // combine
        let len = a_coef.len();
        let third_len = len / 3;
        let mut a_point = vec![0; len];  // TODO trick: unsafe { Vec.set_len() }
        for i in 0..third_len {

            let j = i;
            let x = mod_pow(omega, j as u32, prime);
            let x_squared = (x * x) % prime;
            a_point[j] = (b_point[i] + x * c_point[i] + x_squared * d_point[i]) % prime;

            let j = i + third_len;
            let x = mod_pow(omega, j as u32, prime);
            let x_squared = (x * x) % prime;
            a_point[j] = (b_point[i] + x * c_point[i] + x_squared * d_point[i]) % prime;

            let j = i + third_len + third_len;
            let x = mod_pow(omega, j as u32, prime);
            let x_squared = (x * x) % prime;
            a_point[j] = (b_point[i] + x * c_point[i] + x_squared * d_point[i]) % prime;
        }

        // return
        a_point
    }
}

pub fn fft3_inverse(a_point: &[i64], omega: i64, prime: i64) -> Vec<i64> {
    let omega_inv = mod_inverse(omega, prime);
    let len = a_point.len();
    let len_inv = mod_inverse(len as i64, prime);
    let scaled_a_coef = fft3(a_point, omega_inv, prime);
    let a_coef = scaled_a_coef.iter().map(|x| { x * len_inv % prime }).collect();
    a_coef
}

#[test]
fn test_fft3() {
    // field is Z_433 in which 150 is an 9th root of unity
    let prime = 433;
    let omega = 150;

    let a_coef = vec![1,2,3,4,5,6,7,8,9];
    let a_point = fft3(&a_coef, omega, prime);
    assert_eq!(a_point, vec![45, 404, 407, 266, 377, 47, 158, 17, 20])
}

#[test]
fn test_fft3_inverse() {
    // field is Z_433 in which 150 is an 9th root of unity
    let prime = 433;
    let omega = 150;

    let a_point = vec![45, 404, 407, 266, 377, 47, 158, 17, 20];
    let a_coef = fft3_inverse(&a_point, omega, prime);
    assert_eq!(a_coef, vec![1,2,3,4,5,6,7,8,9])
}


pub fn lagrange_interpolation_at_zero(points: &[i64], values: &[i64], prime: i64) -> i64 {
    assert_eq!(points.len(), values.len());
    // Lagrange interpolation for point 0
    let mut acc = 0i64;
    for i in 0..values.len() {
        let xi = points[i];
        let yi = values[i];
        let mut num = 1i64;
        let mut denum = 1i64;
        for j in 0..values.len() {
            if j != i {
                let xj = points[j];
                num = (num * xj) % prime;
                denum = (denum * (xj - xi)) % prime;
            }
        }
        acc = (acc + yi * num * mod_inverse(denum, prime)) % prime;
    }
    acc
}


pub struct NewtonPolynomial<'a> {
    points: &'a[i64],
    coefficients: Vec<i64>
}

pub fn newton_interpolation_general<'a>(points: &'a[i64], values: &[i64], prime: i64) -> NewtonPolynomial<'a> {
    let coefficients = compute_newton_coefficients(points, values, prime);
    NewtonPolynomial {
        points: points,
        coefficients: coefficients
    }
}

#[test]
fn test_newton_interpolation_general() {
    let prime = 17;

    let poly = [1,2,3,4];
    let points = vec![5, 6, 7, 8, 9];
    let values: Vec<i64> = points.iter().map(|&point| mod_evaluate_polynomial(&poly, point, prime)).collect();
    assert_eq!(values, vec![8, 16, 4, 13, 16]);

    let recovered_poly = newton_interpolation_general(&points, &values, prime);
    let recovered_values: Vec<i64> = points.iter().map(|&point| newton_evaluate(&recovered_poly, point, prime)).collect();
    assert_eq!(recovered_values, values);

    assert_eq!(newton_evaluate(&recovered_poly, 10, prime), 3);
    assert_eq!(newton_evaluate(&recovered_poly, 11, prime), -2);
    assert_eq!(newton_evaluate(&recovered_poly, 12, prime), 8);
}

pub fn newton_evaluate(poly: &NewtonPolynomial, point: i64, prime: i64) -> i64 {
    // compute Newton points
    let mut newton_points = vec![1];
    for i in 0..poly.points.len()-1 {
        let diff = (point - poly.points[i]) % prime;
        let product = (newton_points[i] * diff) % prime;
        newton_points.push(product);
    }
    let ref newton_coefs = poly.coefficients;
    // sum up
    newton_coefs.iter()
        .zip(newton_points)
        .map(|(coef, point)| (coef * point) % prime)
        .fold(0, |a, b| (a + b) % prime)
}

fn compute_newton_coefficients(points: &[i64], values: &[i64], prime: i64) -> Vec<i64> {
    assert_eq!(points.len(), values.len());

    let mut store: Vec<(usize, usize, i64)> = values.iter().enumerate().map(|(index, &value)| (index, index, value)).collect();

    for j in 1..store.len() {
        for i in (j..store.len()).rev() {
            let index_lower = store[i-1].0;
            let index_upper = store[i].1;

            let point_lower = points[index_lower];
            let point_upper = points[index_upper];
            let point_diff = (point_upper - point_lower) % prime;
            let point_diff_inverse = mod_inverse(point_diff, prime);

            let coef_lower = store[i-1].2;
            let coef_upper = store[i].2;
            let coef_diff = (coef_upper - coef_lower) % prime;

            let fraction = (coef_diff * point_diff_inverse) % prime;

            store[i] = (index_lower, index_upper, fraction);
        }
    }

    store.iter().map(|&(_, _, v)| v).collect()
}

#[test]
fn test_compute_newton_coefficients() {
    let points = vec![5, 6, 7, 8, 9];
    let values = vec![8, 16, 4, 13, 16];
    let prime = 17;

    let coefficients = compute_newton_coefficients(&points, &values, prime);
    assert_eq!(coefficients, vec![8,8,-10,4,0]);
}


pub fn positivise(values: &[i64], prime: i64) -> Vec<i64> {
    values.iter()
        .map(|&value| if value < 0 { value + prime } else { value })
        .collect()
}

// deprecated
// fn mod_evaluate_polynomial_naive(coefficients: &[i64], point: i64, prime: i64) -> i64 {
//     // evaluate naively
//     coefficients.iter()
//        .enumerate()
//        .map(|(deg, coef)| (coef * mod_pow(point, deg as u32, prime)) % prime)
//        .fold(0, |a, b| (a + b) % prime)
// }
//
// #[test]
// fn test_mod_evaluate_polynomial_naive() {
//     let poly = vec![1,2,3,4,5,6];
//     let point = 5;
//     let prime = 17;
//     assert_eq!(mod_evaluate_polynomial_naive(&poly, point, prime), 4);
// }

pub fn mod_evaluate_polynomial(coefficients: &[i64], point: i64, prime: i64) -> i64 {
    // evaluate using Horner's rule
    //  - to combine with fold we consider the coefficients in reverse order
    let mut reversed_coefficients = coefficients.iter().rev();
    // manually split due to fold insisting on an initial value
    let head = *reversed_coefficients.next().unwrap();
    let tail = reversed_coefficients;
    tail.fold(head, |partial, coef| (partial * point + coef) % prime)
}

#[test]
fn test_mod_evaluate_polynomial() {
    let poly = vec![1,2,3,4,5,6];
    let point = 5;
    let prime = 17;
    assert_eq!(mod_evaluate_polynomial(&poly, point, prime), 4);
}
