//! Mpz precision tests — multiplication (all 6 backends), division, pow.
//! Uses rug (GMP) as the oracle for correctness validation.

use mynum::mpz::{Mpz, MpzMultiplicationConfig, MultiplicationBackend};
use rug::ops::Pow;
use rug::Integer;
use std::str::FromStr;

fn random_hex(bytes: usize) -> String {
    use std::collections::hash_map::RandomState;
    use std::hash::{BuildHasher, Hasher};
    let mut s = String::with_capacity(bytes * 2);
    let mut hasher = RandomState::new().build_hasher();
    for _ in 0..bytes {
        hasher.write_u8(0);
        let b = (hasher.finish() & 0xFF) as u8;
        s.push_str(&format!("{:02x}", b));
    }
    s.trim_start_matches('0').to_string()
}

fn test_mul_vs_rug(a_hex: &str, b_hex: &str) {
    let a_rug = Integer::from_str_radix(a_hex, 16).unwrap();
    let b_rug = Integer::from_str_radix(b_hex, 16).unwrap();
    let expected = (a_rug * b_rug).to_string_radix(10);

    let a = Mpz::from_str(a_hex, 16).unwrap();
    let b = Mpz::from_str(b_hex, 16).unwrap();
    let result = a.mul(&b);
    assert_eq!(
        result.to_string(10),
        expected,
        "Mul mismatch: {} * {}",
        a_hex,
        b_hex
    );
}

#[test]
fn test_mul_schoolbook() {
    MpzMultiplicationConfig::set_global_backend(MultiplicationBackend::Schoolbook);
    for _ in 0..10 {
        test_mul_vs_rug(&random_hex(64), &random_hex(64));
    }
    for _ in 0..5 {
        test_mul_vs_rug(&random_hex(128), &random_hex(128));
    }
}

#[test]
fn test_mul_karatsuba() {
    MpzMultiplicationConfig::set_global_backend(MultiplicationBackend::Karatsuba);
    for _ in 0..10 {
        test_mul_vs_rug(&random_hex(256), &random_hex(256));
    }
    for _ in 0..5 {
        test_mul_vs_rug(&random_hex(512), &random_hex(256));
    }
    test_mul_vs_rug(&random_hex(800), &random_hex(200));
}

#[test]
fn test_mul_toom3() {
    MpzMultiplicationConfig::set_global_backend(MultiplicationBackend::ToomCook3);
    for _ in 0..10 {
        test_mul_vs_rug(&random_hex(2048), &random_hex(2048));
    }
    for _ in 0..5 {
        test_mul_vs_rug(&random_hex(4000), &random_hex(2400));
    }
    test_mul_vs_rug(&random_hex(4096), &random_hex(4096));
}

#[test]
#[cfg_attr(tarpaulin, ignore)]
fn test_mul_toom4() {
    MpzMultiplicationConfig::set_global_backend(MultiplicationBackend::ToomCook4);
    for _ in 0..10 {
        test_mul_vs_rug(&random_hex(8192), &random_hex(8192));
    }
    for _ in 0..3 {
        test_mul_vs_rug(&random_hex(16000), &random_hex(12000));
    }
}

#[test]
#[cfg_attr(tarpaulin, ignore)]
fn test_mul_fft() {
    MpzMultiplicationConfig::set_global_backend(MultiplicationBackend::FFT);
    for _ in 0..10 {
        test_mul_vs_rug(&random_hex(32768), &random_hex(32768));
    }
    for _ in 0..3 {
        test_mul_vs_rug(&random_hex(64000), &random_hex(32000));
    }
}

#[test]
#[cfg_attr(tarpaulin, ignore)]
fn test_mul_ntt() {
    MpzMultiplicationConfig::set_global_backend(MultiplicationBackend::NTT);
    for _ in 0..10 {
        test_mul_vs_rug(&random_hex(131072), &random_hex(131072));
    }
    for _ in 0..2 {
        test_mul_vs_rug(&random_hex(200000), &random_hex(100000));
    }
}

#[test]
#[cfg_attr(tarpaulin, ignore)]
fn test_mul_adaptive() {
    MpzMultiplicationConfig::set_global_backend(MultiplicationBackend::Adaptive);
    for limbs in &[50, 200, 1000, 4000, 16000, 30000] {
        for _ in 0..5 {
            test_mul_vs_rug(&random_hex(*limbs * 8), &random_hex(*limbs * 8));
        }
    }
}

#[test]
fn test_mul_edge_cases() {
    let zero = Mpz::new();
    let one = Mpz::from_i64(1);
    let big = Mpz::from_str("ffffffffffffffffffffffffffffffffffffffff", 16).unwrap();
    assert_eq!(zero.mul(&big).to_string(10), "0");
    assert_eq!(one.mul(&big), big);
    let big_sq = big.mul(&big);
    assert!(big_sq.to_string(10).len() > 60, "large product too small");
    // Identity: a * 1 = a
    let a = Mpz::from_str("123456789012345678901234567890", 10).unwrap();
    assert_eq!(a.mul(&one), a);
    assert_eq!(a.mul(&zero).to_string(10), "0");
}

#[test]
fn test_mul_commutative() {
    for _ in 0..20 {
        let a_hex = random_hex(128);
        let b_hex = random_hex(128);
        let a = Mpz::from_str(&a_hex, 16).unwrap();
        let b = Mpz::from_str(&b_hex, 16).unwrap();
        assert_eq!(a.mul(&b), b.mul(&a), "Multiplication not commutative");
    }
}

// ── Division tests ──

#[test]
fn test_div_precision() {
    for _ in 0..30 {
        let a_hex = random_hex(256);
        let mut b_hex = random_hex(128);
        if b_hex.is_empty() || b_hex == "0" {
            b_hex = "1".to_string();
        }
        let a_rug = Integer::from_str_radix(&a_hex, 16).unwrap();
        let b_rug = Integer::from_str_radix(&b_hex, 16).unwrap();
        let (q_rug, r_rug) = a_rug.div_rem(b_rug.clone());
        let a = Mpz::from_str(&a_hex, 16).unwrap();
        let b = Mpz::from_str(&b_hex, 16).unwrap();
        let (q, r) = a.div_rem(&b).unwrap();
        assert_eq!(
            q.to_string(10),
            q_rug.to_string_radix(10),
            "div quotient mismatch"
        );
        assert_eq!(
            r.to_string(10),
            r_rug.to_string_radix(10),
            "div remainder mismatch"
        );
    }
}

#[test]
fn test_div_large() {
    for _ in 0..10 {
        let a_hex = random_hex(4096);
        let mut b_hex = random_hex(2048);
        if b_hex.is_empty() {
            b_hex = "1".to_string();
        }
        let a_rug = Integer::from_str_radix(&a_hex, 16).unwrap();
        let b_rug = Integer::from_str_radix(&b_hex, 16).unwrap();
        let (q_rug, r_rug) = a_rug.div_rem(b_rug.clone());
        let a = Mpz::from_str(&a_hex, 16).unwrap();
        let b = Mpz::from_str(&b_hex, 16).unwrap();
        let (q, r) = a.div_rem(&b).unwrap();
        assert_eq!(q.to_string(10), q_rug.to_string_radix(10));
        assert_eq!(r.to_string(10), r_rug.to_string_radix(10));
    }
}

#[test]
fn test_div_edge() {
    let a = Mpz::from_i64(42);
    let one = Mpz::from_i64(1);
    let (q, r) = a.div_rem(&one).unwrap();
    assert_eq!(q, a);
    assert_eq!(r.to_string(10), "0");
    // a / a = 1, a % a = 0
    let (q, r) = a.div_rem(&a).unwrap();
    assert_eq!(q, one);
    assert_eq!(r.to_string(10), "0");
}

// ── Pow tests ──

#[test]
fn test_pow_u32_precision() {
    for _ in 0..20 {
        let base_hex = random_hex(16);
        let base = Mpz::from_str(&base_hex, 16).unwrap();
        let exp: u32 = 3;
        let rug_base = Integer::from_str_radix(&base_hex, 16).unwrap();
        let expected = rug_base.pow(exp).to_string_radix(10);
        let result = base.pow_u32(exp);
        assert_eq!(result.to_string(10), expected);
    }
}

#[test]
fn test_pow_big_vs_rug() {
    let cases = [
        ("2", "10"),
        ("3", "20"),
        ("5", "15"),
        ("7", "12"),
        ("10", "8"),
    ];
    for (base_str, exp_str) in &cases {
        let base = Mpz::from_str(base_str, 10).unwrap();
        let exp = Mpz::from_str(exp_str, 10).unwrap();
        let result = base.pow(&exp);
        let rug_base = Integer::from_str(base_str).unwrap();
        let rug_exp = u32::from_str(exp_str).unwrap();
        let expected = rug_base.pow(rug_exp).to_string_radix(10);
        assert_eq!(result.to_string(10), expected);
    }
}

#[test]
fn test_pow_zero_one() {
    let zero = Mpz::new();
    let one = Mpz::from_i64(1);
    let two = Mpz::from_i64(2);
    let ten = Mpz::from_i64(10);
    assert_eq!(two.pow(&zero), one);
    assert_eq!(two.pow(&one), two);
    assert_eq!(zero.pow(&ten).to_string(10), "0");
    assert_eq!(one.pow(&ten), one);
}

#[test]
fn test_pow_large() {
    let two = Mpz::from_i64(2);
    let hundred = Mpz::from_i64(100);
    let result = two.pow(&hundred);
    let rug_two = Integer::from(2u32);
    let expected = rug_two.pow(100).to_string_radix(10);
    assert_eq!(result.to_string(10), expected);
    // Verify a known value: 2^100
    assert!(result.to_string(10).len() > 30);
}

// ── GCD / LCM tests ──

#[test]
fn test_gcd_precision() {
    for _ in 0..20 {
        let a_hex = random_hex(128);
        let b_hex = random_hex(128);
        let a = Mpz::from_str(&a_hex, 16).unwrap();
        let b = Mpz::from_str(&b_hex, 16).unwrap();
        let a_rug = Integer::from_str_radix(&a_hex, 16).unwrap();
        let b_rug = Integer::from_str_radix(&b_hex, 16).unwrap();
        let g = a.gcd(&b);
        let expected = a_rug.gcd(&b_rug).to_string_radix(10);
        assert_eq!(g.to_string(10), expected);
    }
}
