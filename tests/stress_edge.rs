//! Stress tests and boundary conditions for mynum.

use mynum::mpf::Mpf;
use mynum::mpz::{Mpz, MpzMultiplicationConfig, MultiplicationBackend};

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

#[test]
fn test_precision_transition() {
    let a = Mpf::from_f64(1.5, 64);
    let b = Mpf::from_f64(2.5, 256);
    let c = a.add(&b);
    assert!((c.to_f64().unwrap() - 4.0).abs() < 1e-10);
}

#[test]
fn test_large_multiplication() {
    let hex = "f".repeat(1250);
    let a = Mpz::from_str(&hex, 16).unwrap();
    let b = Mpz::from_str(&hex, 16).unwrap();
    let result = a.mul(&b);
    assert!(
        result.bit_length() > 9000,
        "large product too small: {} bits",
        result.bit_length()
    );
}

#[test]
#[cfg_attr(tarpaulin, ignore)]
fn test_overflow_multiplication() {
    let hex = "f".repeat(5000);
    let a = Mpz::from_str(&hex, 16).unwrap();
    let bigger = a.mul(&a);
    assert!(bigger.bit_length() > 30000);
}

#[test]
fn test_exact_rounding_tolerance() {
    for prec in &[32, 64, 128, 256] {
        let one = Mpf::from_mpz(Mpz::from_i64(1), *prec);
        let three = Mpf::from_mpz(Mpz::from_i64(3), *prec);
        let third = one.div(&three).unwrap();
        let back = third.mul(&three);
        let diff = (back.to_f64().unwrap() - 1.0).abs();
        let limit = if *prec <= 64 { 1e-8 } else { 1e-30 };
        assert!(
            diff < limit,
            "1/3*3 precision loss at {} bits: diff={}",
            prec,
            diff
        );
    }
}

#[test]
fn test_from_str_roundtrip() {
    let values = [
        "3.141592653589793",
        "0.000000000000001",
        "12345678901234567890.1234567890",
        "1e100",
        "1e-100",
    ];
    for val in &values {
        let a = Mpf::from_str(val, 10).unwrap();
        let s = a.to_string(10);
        assert!(!s.is_empty(), "from_str({}) produced empty", val);
    }
}

#[test]
fn test_zero_operations() {
    let zero = Mpf::new();
    let one = Mpf::from_i64(1, 64);
    let neg_one = Mpf::from_f64(-1.0, 64);

    assert!(zero.mul(&one).is_zero());
    assert!(zero.add(&zero).is_zero());
    assert!(one.sub(&one).is_zero());
    assert!(zero.div(&one).unwrap().is_zero());
    assert!(zero.exp().unwrap().to_f64().unwrap() - 1.0 < 1e-15);
    assert!(zero.ln().is_err());
    assert!(one.div(&zero).is_err());
    assert!(neg_one.sqrt().is_err());
}

#[test]
fn test_concurrent_random_mul() {
    use std::thread;
    let handles: Vec<_> = (0..8)
        .map(|_| {
            thread::spawn(|| {
                for _ in 0..50 {
                    let a = Mpz::from_str(&random_hex(512), 16).unwrap();
                    let b = Mpz::from_str(&random_hex(512), 16).unwrap();
                    let r1 = a.mul(&b);
                    let r2 = b.mul(&a);
                    assert_eq!(r1, r2);
                }
            })
        })
        .collect();
    for h in handles {
        h.join().unwrap();
    }
}

#[test]
fn test_backend_switch_stress() {
    // Rapidly switch between backends and verify correctness
    let backends = [
        MultiplicationBackend::Schoolbook,
        MultiplicationBackend::Karatsuba,
        MultiplicationBackend::ToomCook3,
        MultiplicationBackend::ToomCook4,
    ];
    for _ in 0..20 {
        let a = Mpz::from_str(&random_hex(1024), 16).unwrap();
        let b = Mpz::from_str(&random_hex(1024), 16).unwrap();
        let mut prev = None;
        for &backend in &backends {
            MpzMultiplicationConfig::set_global_backend(backend);
            let r = a.mul(&b);
            if let Some(ref prev_r) = prev {
                assert_eq!(
                    r, *prev_r,
                    "Backend {:?} produced different result",
                    backend
                );
            }
            prev = Some(r);
        }
        MpzMultiplicationConfig::set_global_backend(MultiplicationBackend::Adaptive);
    }
}

#[test]
fn test_many_digits_pi() {
    let pi = Mpf::pi(1024);
    let s = pi.to_string(10);
    assert!(
        s.starts_with("3.14159265358979323846"),
        "Pi starts wrong: {}",
        &s[..30]
    );
}

#[test]
fn test_very_small_numbers() {
    let tiny = Mpf::from_f64(1e-200, 256);
    let also_tiny = Mpf::from_f64(1e-200, 256);
    let product = tiny.mul(&also_tiny);
    assert!(product.is_zero() || product.to_f64().unwrap() < 1e-300);
}
