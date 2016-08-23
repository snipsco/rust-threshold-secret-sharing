use super::Field;

fn fft2_in_place_rearrange<F: Field>(_zp: &F, data: &mut [F::U]) {
    let mut target = 0;
    for pos in 0..data.len() {
        if target > pos {
            data.swap(target, pos)
        }
        let mut mask = data.len() >> 1;
        while target & mask != 0 {
            target &= !mask;
            mask >>= 1;
        }
        target |= mask;
    }
}

fn fft2_in_place_compute<F: Field>(zp: &F, data: &mut [F::U], omega: F::U) {
    let mut depth = 0usize;
    let big_omega = zp.optimize(zp.qpow(omega, data.len() as u32 / 2));
    while 1usize << depth < data.len() {
        let step = 1usize << depth;
        let jump = 2 * step;
        let factor_stride = zp.optimize(zp.qpow(omega, (data.len() / step / 2) as u32));
        let mut factor = zp.one();
        for group in 0usize..step {
            let mut pair = group;
            while pair < data.len() {
                let (x, y) = (data[pair], zp.mul(data[pair + step], factor));

                data[pair] = zp.add(x, y);
                data[pair + step] = zp.add(x, zp.mul(big_omega, y));

                pair += jump;
            }
            factor = zp.optimize(zp.mul(factor, factor_stride));
        }
        depth += 1;
    }
}

pub fn fft2<F: Field>(zp: &F, data: &mut [F::U], omega: F::U) {
    fft2_in_place_rearrange(zp, &mut *data);
    fft2_in_place_compute(zp, &mut *data, omega);
}

pub fn fft2_inverse<F: Field>(zp: &F, data: &mut [F::U], omega: F::U) {
    let omega_inv = zp.inv(omega);
    let len = data.len();
    let len_inv = zp.inv(zp.from(len as u64));
    fft2(zp, data, omega_inv);
    for mut x in data {
        *x = zp.mul(*x, len_inv);
    }
}

fn trigits_len(n: usize) -> usize {
    let mut result = 1;
    let mut value = 3;
    while value < n + 1 {
        result += 1;
        value *= 3;
    }
    result
}

fn fft3_in_place_rearrange<F: Field>(_zp: &F, data: &mut [F::U]) {
    let mut target = 0isize;
    let trigits_len = trigits_len(data.len() - 1);
    let mut trigits: Vec<u8> = ::std::iter::repeat(0).take(trigits_len).collect();
    let powers: Vec<isize> = (0..trigits_len).map(|x| 3isize.pow(x as u32)).rev().collect();
    for pos in 0..data.len() {
        if target as usize > pos {
            data.swap(target as usize, pos)
        }
        for pow in 0..trigits_len {
            if trigits[pow] < 2 {
                trigits[pow] += 1;
                target += powers[pow];
                break;
            } else {
                trigits[pow] = 0;
                target -= 2 * powers[pow];
            }
        }
    }
}

fn fft3_in_place_compute<F: Field>(zp: &F, data: &mut [F::U], omega: F::U) {
    let mut step = 1;
    let big_omega = zp.optimize(zp.qpow(omega, (data.len() as u32 / 3)));
    let big_omega_sq = zp.optimize(zp.mul(big_omega, big_omega));
    while step < data.len() {
        let jump = 3 * step;
        let factor_stride = zp.optimize(zp.qpow(omega, (data.len() / step / 3) as u32));
        let mut factor = zp.one();
        for group in 0usize..step {
            let factor_sq = zp.mul(factor, factor);
            let mut pair = group;
            while pair < data.len() {
                let (x, y, z) = (data[pair],
                                 zp.mul(data[pair + step], factor),
                                 zp.mul(data[pair + 2 * step], factor_sq));

                data[pair] = zp.add(zp.add(x, y), z);
                data[pair + step] =
                    zp.add(zp.add(x, zp.mul(big_omega, y)), zp.mul(big_omega_sq, z));
                data[pair + 2 * step] =
                    zp.add(zp.add(x, zp.mul(big_omega_sq, y)), zp.mul(big_omega, z));

                pair += jump;
            }
            factor = zp.optimize(zp.mul(factor, factor_stride));
        }
        step = jump;
    }
}

pub fn fft3<F: Field>(zp: &F, data: &mut [F::U], omega: F::U) {
    fft3_in_place_rearrange(zp, &mut *data);
    fft3_in_place_compute(zp, &mut *data, omega);
}

pub fn fft3_inverse<F: Field>(zp: &F, data: &mut [F::U], omega: F::U) {
    let omega_inv = zp.inv(omega);
    let len_inv = zp.inv(zp.from(data.len() as u64));
    fft3(zp, data, omega_inv);
    for mut x in data {
        *x = zp.mul(*x, len_inv);
    }
}

#[cfg(test)]
use super::ZpField;

#[cfg(test)]
pub fn test_fft2<F: ZpField>() {
    // field is Z_433 in which 354 is an 8th root of unity
    let zp = F::new(433);
    let omega = zp.from(354);

    let mut data = zp.from_v(&[1, 2, 3, 4, 5, 6, 7, 8]);
    fft2(&zp, &mut data, omega);
    assert_eq!(zp.back_v(&data), [36, 303, 146, 3, 429, 422, 279, 122]);
}

#[cfg(test)]
pub fn test_fft2_inverse<F: ZpField>() {
    // field is Z_433 in which 354 is an 8th root of unity
    let zp = F::new(433);
    let omega = zp.from(354);

    let mut data = zp.from_v(&[36, 303, 146, 3, 429, 422, 279, 122]);
    fft2_inverse(&zp, &mut *data, omega);
    assert_eq!(zp.back_v(&data), [1, 2, 3, 4, 5, 6, 7, 8])
}

#[cfg(test)]
pub fn test_fft2_big<F: ZpField>() {
    let zp = F::new(5038849);
    let omega = zp.from(4318906);

    let mut data: Vec<_> = (0..256).map(|a| zp.from(a)).collect();
    fft2(&zp, &mut *data, omega);
    fft2_inverse(&zp, &mut data, omega);

    assert_eq!(zp.back_v(&data), (0..256).collect::<Vec<_>>());
}

#[cfg(test)]
pub fn test_fft3<F: ZpField>() {
    // field is Z_433 in which 150 is an 9th root of unity
    let zp = F::new(433);
    let omega = zp.from(150);

    let mut data = zp.from_v(&[1, 2, 3, 4, 5, 6, 7, 8, 9]);
    fft3(&zp, &mut data, omega);
    assert_eq!(zp.back_v(&data), [45, 404, 407, 266, 377, 47, 158, 17, 20]);
}

#[cfg(test)]
pub fn test_fft3_inverse<F: ZpField>() {
    // field is Z_433 in which 150 is an 9th root of unity
    let zp = F::new(433);
    let omega = zp.from(150);

    let mut data = zp.from_v(&[45, 404, 407, 266, 377, 47, 158, 17, 20]);
    fft3_inverse(&zp, &mut *data, omega);
    assert_eq!(zp.back_v(&data), [1, 2, 3, 4, 5, 6, 7, 8, 9])
}

#[cfg(test)]
pub fn test_fft3_big<F: ZpField>() {
    let zp = F::new(5038849);
    let omega = zp.from(1814687);

    let mut data: Vec<_> = (0..19683).map(|a| zp.from(a)).collect();
    fft3(&zp, &mut data, omega);
    fft3_inverse(&zp, &mut data, omega);

    assert_eq!(zp.back_v(&data), (0..19683).collect::<Vec<_>>());
}
