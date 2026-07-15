//! Mpf elementary function precision tests — validated against MPFR via rug.
//! Uses deterministic test values for reproducible regression detection.

use mynum::mpf::Mpf;
use rug::Float;

fn check_value(name: &str, prec: usize, x: f64, result: f64, expected: f64, tol: f64) {
    let diff = (result - expected).abs();
    let mag = expected.abs().max(1e-10);
    let rel = diff / mag;
    assert!(
        rel < tol,
        "{} error at {} bits x={}: mynum={}, rug={}, rel_err={}",
        name,
        prec,
        x,
        result,
        expected,
        rel
    );
}

fn test_at_point(
    name: &str,
    x: f64,
    prec: usize,
    mynum_fn: fn(&Mpf) -> mynum::Result<Mpf>,
    rug_fn: fn(Float) -> Float,
    tol: f64,
) {
    let a = Mpf::from_f64(x, prec);
    let result = mynum_fn(&a).unwrap();
    let r = Float::with_val(prec as u32, x);
    let expected = rug_fn(r);
    check_value(
        name,
        prec,
        x,
        result.to_f64().unwrap(),
        expected.to_f64(),
        tol,
    );
}

macro_rules! test_points {
    ($name:ident, $fn:expr, $rug:expr, $tol:expr, [$(($x:expr, $desc:expr)),+]) => {
        #[test]
        fn $name() {
            for &prec in &[64, 128] {
                let epsilon = if prec <= 64 { $tol } else { $tol * 5e-2 };
                $(
                    test_at_point(stringify!($name), $x, prec, $fn, $rug, epsilon);
                )+
            }
        }
    };
}

test_points!(
    test_sin_vs_mpfr,
    Mpf::sin,
    |f| f.sin(),
    1e-7,
    [
        (0.0, "zero"),
        (1.0, "one"),
        (-1.0, "neg_one"),
        (0.5, "half"),
        (2.0, "two"),
        (3.0, "three")
    ]
);

test_points!(
    test_cos_vs_mpfr,
    Mpf::cos,
    |f| f.cos(),
    1e-7,
    [
        (0.0, "zero"),
        (1.0, "one"),
        (-1.0, "neg_one"),
        (0.5, "half"),
        (2.0, "two"),
        (3.0, "three")
    ]
);

test_points!(
    test_tan_vs_mpfr,
    Mpf::tan,
    |f| f.tan(),
    1e-4,
    [
        (0.0, "zero"),
        (0.5, "half"),
        (1.0, "one"),
        (-0.5, "neg_half")
    ]
);

test_points!(
    test_exp_vs_mpfr,
    Mpf::exp,
    |f| f.exp(),
    1e-7,
    [
        (0.0, "zero"),
        (1.0, "one"),
        (-1.0, "neg_one"),
        (2.0, "two"),
        (10.0, "ten"),
        (-5.0, "neg_five")
    ]
);

test_points!(
    test_ln_vs_mpfr,
    Mpf::ln,
    |f| f.ln(),
    1e-7,
    [
        (1.0, "one"),
        (2.0, "two"),
        (10.0, "ten"),
        (0.5, "half"),
        (100.0, "hundred")
    ]
);

test_points!(
    test_asin_vs_mpfr,
    Mpf::asin,
    |f| f.asin(),
    1e-4,
    [
        (0.0, "zero"),
        (0.5, "half"),
        (-0.5, "neg_half"),
        (0.9, "near_one")
    ]
);

test_points!(
    test_acos_vs_mpfr,
    Mpf::acos,
    |f| f.acos(),
    1e-4,
    [
        (0.0, "zero"),
        (0.5, "half"),
        (-0.5, "neg_half"),
        (0.9, "near_one")
    ]
);

test_points!(
    test_atan_vs_mpfr,
    Mpf::atan,
    |f| f.atan(),
    1e-4,
    [
        (0.0, "zero"),
        (1.0, "one"),
        (-1.0, "neg_one"),
        (10.0, "ten")
    ]
);

// ── Edge cases ──

#[test]
fn test_sin_zero() {
    let x = Mpf::new();
    let s = x.sin().unwrap();
    assert!(
        s.is_zero(),
        "sin(0) should be 0, got {}",
        s.to_f64().unwrap()
    );
}

#[test]
fn test_cos_zero() {
    let x = Mpf::new();
    let c = x.cos().unwrap();
    assert!(
        (c.to_f64().unwrap() - 1.0).abs() < 1e-15,
        "cos(0) should be 1"
    );
}

#[test]
fn test_exp_zero() {
    let x = Mpf::new();
    let e = x.exp().unwrap();
    assert!(
        (e.to_f64().unwrap() - 1.0).abs() < 1e-15,
        "exp(0) should be 1"
    );
}

#[test]
fn test_ln_one() {
    let x = Mpf::from_i64(1, 64);
    let l = x.ln().unwrap();
    assert!(l.is_zero(), "ln(1) should be 0");
}

#[test]
fn test_ln_negative() {
    let x = Mpf::from_f64(-1.0, 64);
    assert!(x.ln().is_err(), "ln(-1) should fail");
}

#[test]
fn test_sqrt_zero() {
    let x = Mpf::new();
    assert!(x.sqrt().unwrap().is_zero());
}

#[test]
fn test_sqrt_negative() {
    let x = Mpf::from_f64(-1.0, 64);
    assert!(x.sqrt().is_err(), "sqrt(-1) should fail");
}

#[test]
fn test_asin_domain() {
    let x = Mpf::from_f64(2.0, 64);
    assert!(x.asin().is_err(), "asin(2) should fail");
    let x = Mpf::from_f64(-2.0, 64);
    assert!(x.asin().is_err(), "asin(-2) should fail");
}

#[test]
fn test_tan_near_pi_half() {
    let pi_half = Mpf::pi(128).div(&Mpf::from_i64(2, 128)).unwrap();
    let t = pi_half.tan().unwrap();
    let abs_t = t.abs();
    assert!(
        abs_t.to_f64().unwrap() > 1e8,
        "tan(pi/2) should be large, got {}",
        t.to_f64().unwrap()
    );
}

#[test]
fn test_precision_stability() {
    let x = Mpf::from_f64(1.5, 128);
    let mut val = x.clone();
    for _ in 0..100 {
        val = val.sin().unwrap();
        val = val.asin().unwrap();
    }
    let diff = (val.to_f64().unwrap() - 1.5).abs();
    assert!(
        diff < 0.5,
        "Precision drift after 100 sin/asin cycles: diff={}",
        diff
    );
}
