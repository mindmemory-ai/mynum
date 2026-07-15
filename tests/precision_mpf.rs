//! Mpf precision tests — arithmetic operations validated against MPFR via rug.
//! Tests cover: add, sub, mul, div, sqrt at 64/128/256-bit precision.

use mynum::mpf::Mpf;
use rug::Float;

fn random_f64(min: f64, max: f64) -> f64 {
    use std::collections::hash_map::RandomState;
    use std::hash::{BuildHasher, Hasher};
    let mut h = RandomState::new().build_hasher();
    h.write_u8(0);
    let r = (h.finish() as f64) / (u64::MAX as f64);
    min + r * (max - min)
}

/// Test a binary Mpf operation against MPFR.
fn test_binary_op(
    name: &str,
    precisions: &[usize],
    mynum_op: fn(&Mpf, &Mpf) -> Mpf,
    rug_op: fn(Float, &Float) -> Float,
    tolerance: f64,
) {
    for &prec in precisions {
        for _ in 0..30 {
            let a_val = random_f64(-100.0, 100.0);
            let b_val = if name == "div" {
                random_f64(0.1, 100.0)
            } else {
                random_f64(-100.0, 100.0)
            };

            let a_m = Mpf::from_f64(a_val, prec);
            let b_m = Mpf::from_f64(b_val, prec);
            let result = mynum_op(&a_m, &b_m);

            let a_r = Float::with_val(prec as u32, a_val);
            let b_r = Float::with_val(prec as u32, b_val);
            let expected = rug_op(a_r, &b_r);

            let result_f = result.to_f64().unwrap();
            let expected_f = expected.to_f64();
            let diff = (result_f - expected_f).abs();
            let mag = expected_f.abs().max(1.0);
            let rel_err = diff / mag;

            assert!(
                rel_err < tolerance,
                "{} error at {} bits: mynum={}, rug={}, rel_err={}",
                name,
                prec,
                result_f,
                expected_f,
                rel_err
            );
        }
    }
}

fn rug_add(a: Float, b: &Float) -> Float {
    a + b
}
fn rug_sub(a: Float, b: &Float) -> Float {
    a - b
}
fn rug_mul(a: Float, b: &Float) -> Float {
    a * b
}
fn rug_div(a: Float, b: &Float) -> Float {
    a / b
}

#[test]
fn test_mpf_add_precision() {
    test_binary_op("add", &[64, 128, 256], |a, b| a.add(b), rug_add, 1e-14);
}

#[test]
fn test_mpf_sub_precision() {
    test_binary_op("sub", &[64, 128, 256], |a, b| a.sub(b), rug_sub, 1e-14);
}

#[test]
fn test_mpf_mul_precision() {
    test_binary_op("mul", &[64, 128, 256], |a, b| a.mul(b), rug_mul, 1e-13);
}

#[test]
fn test_mpf_div_precision() {
    test_binary_op(
        "div",
        &[64, 128, 256],
        |a, b| a.div(b).unwrap(),
        rug_div,
        1e-12,
    );
}

#[test]
fn test_mpf_sqrt_precision() {
    for &prec in &[64, 128, 256] {
        for _ in 0..30 {
            let val = random_f64(0.0, 1000.0);
            let a = Mpf::from_f64(val, prec);
            let result = a.sqrt().unwrap();
            let expected = Float::with_val(prec as u32, val).sqrt();
            let diff = (result.to_f64().unwrap() - expected.to_f64()).abs();
            let rel_err = diff / expected.to_f64().abs().max(1.0);
            assert!(
                rel_err < 1e-13,
                "sqrt error at {} bits: rel_err={}",
                prec,
                rel_err
            );
        }
    }
}

#[test]
fn test_mpf_add_sub_identity() {
    for &prec in &[64, 128, 256] {
        for _ in 0..20 {
            let a = random_f64(-50.0, 50.0);
            let b = random_f64(-50.0, 50.0);
            let am = Mpf::from_f64(a, prec);
            let bm = Mpf::from_f64(b, prec);
            // (a + b) - b == a
            let result = am.add(&bm).sub(&bm);
            let diff = (result.to_f64().unwrap() - a).abs();
            assert!(
                diff < 1e-12,
                "identity failed at {} bits: diff={}",
                prec,
                diff
            );
        }
    }
}

#[test]
fn test_mpf_mul_div_identity() {
    for &prec in &[64, 128, 256] {
        for _ in 0..20 {
            let a = random_f64(-50.0, 50.0);
            let b = random_f64(0.1, 50.0);
            let am = Mpf::from_f64(a, prec);
            let bm = Mpf::from_f64(b, prec);
            // (a * b) / b == a
            let result = am.mul(&bm).div(&bm).unwrap();
            let diff = (result.to_f64().unwrap() - a).abs();
            assert!(
                diff < 1e-11,
                "mul/div identity failed at {} bits: diff={}",
                prec,
                diff
            );
        }
    }
}

#[test]
fn test_mpf_sqrt_identity() {
    for _ in 0..20 {
        let val = random_f64(0.0, 100.0);
        let a = Mpf::from_f64(val, 128);
        let s = a.sqrt().unwrap();
        let back = s.mul(&s);
        let diff = (back.to_f64().unwrap() - val).abs();
        assert!(
            diff < 1e-12,
            "sqrt identity failed: sqrt({})^2 = {}",
            val,
            back.to_f64().unwrap()
        );
    }
}
