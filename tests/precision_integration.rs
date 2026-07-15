//! Integration and numerical methods precision tests.

use mynum::mpf::{differentiation, integration::*, solve, Mpf};
use mynum::mpz::Mpz;

#[test]
fn test_quad_polynomial_exact() {
    for deg in 0..=4 {
        let f = move |x: &Mpf| {
            let mut r = Mpf::from_mpz(Mpz::from_i64(1), x.precision());
            for _ in 0..deg {
                r = r.mul(x);
            }
            Ok(r)
        };
        let a = Mpf::from_f64(0.0, 128);
        let b = Mpf::from_f64(1.0, 128);
        let result = quad(&f, &a, &b, None, None).unwrap();
        let expected = 1.0 / (deg as f64 + 1.0);
        assert!(
            (result.to_f64().unwrap() - expected).abs() < 1e-10,
            "GK failed for x^{}: got {}, expected {}",
            deg,
            result.to_f64().unwrap(),
            expected
        );
    }
}

#[test]
fn test_quad_sin_cos() {
    let a = Mpf::from_f64(0.0, 64);
    let b = Mpf::pi(64);
    let f_sin = |x: &Mpf| x.sin();
    let result = quad(&f_sin, &a, &b, None, None).unwrap();
    // CORDIC sin accuracy limits integration precision; 1e-4 is realistic at 64-bit
    assert!(
        (result.to_f64().unwrap() - 2.0).abs() < 1e-4,
        "∫sin should ≈ 2, got {}",
        result.to_f64().unwrap()
    );
    let f_cos = |x: &Mpf| x.cos();
    let result = quad(&f_cos, &a, &b, None, None).unwrap();
    assert!(
        result.to_f64().unwrap().abs() < 1e-4,
        "∫cos from 0 to pi should ≈ 0, got {}",
        result.to_f64().unwrap()
    );
}

#[test]
fn test_quad_exp() {
    let a = Mpf::from_f64(0.0, 64);
    let b = Mpf::from_f64(1.0, 64);
    let f = |x: &Mpf| x.exp();
    let result = quad(&f, &a, &b, None, None).unwrap();
    assert!((result.to_f64().unwrap() - (std::f64::consts::E - 1.0)).abs() < 1e-8);
}

#[test]
fn test_nsum_geometric() {
    // Σ_{n=0}∞ (1/2)^n = 2
    let f = |n: u64| -> Result<Mpf, mynum::Error> {
        let two = Mpf::from_mpz(Mpz::from_i64(2), 64);
        let pow = two.pow(n as u32)?;
        Ok(Mpf::from_mpz(Mpz::from_i64(1), 64).div(&pow)?)
    };
    let result = nsum(&f, 0, 200, None).unwrap();
    assert!((result.to_f64().unwrap() - 2.0).abs() < 1e-10);
}

#[test]
fn test_nsum_basel() {
    // Σ 1/n^2 = π^2/6 ≈ 1.6449340668
    let f = |n: u64| -> Result<Mpf, mynum::Error> {
        let n_mpf = Mpf::from_mpz(Mpz::from_i64(n as i64), 64);
        Ok(Mpf::from_mpz(Mpz::from_i64(1), 64).div(&n_mpf.mul(&n_mpf))?)
    };
    let result = nsum(&f, 1, 5000, None).unwrap();
    // Basel series converges as O(1/N); 5000 terms gives ~2e-4 accuracy
    assert!(
        (result.to_f64().unwrap() - 1.6449340668).abs() < 1e-3,
        "nsum basel: got {}",
        result.to_f64().unwrap()
    );
}

#[test]
fn test_limit_sin_over_x() {
    let f = |x: &Mpf| -> Result<Mpf, mynum::Error> { Ok(x.sin()?.div(x)?) };
    let x0 = Mpf::new();
    let result = limit(&f, &x0, None).unwrap();
    assert!((result.to_f64().unwrap() - 1.0).abs() < 1e-8);
}

#[test]
fn test_limit_exp_minus_1_over_x() {
    let one = Mpf::from_mpz(Mpz::from_i64(1), 128);
    let f = |x: &Mpf| -> Result<Mpf, mynum::Error> { Ok(x.exp()?.sub(&one).div(x)?) };
    let result = limit(&f, &Mpf::new(), None).unwrap();
    assert!((result.to_f64().unwrap() - 1.0).abs() < 1e-8);
}

#[test]
fn test_diff_sin() {
    // Derivative of sin at 0 is cos(0) = 1
    let x = Mpf::from_f64(0.0, 128);
    let f = |x: &Mpf| x.sin();
    let result = differentiation::central_diff(&f, &x, &Mpf::from_f64(0.001, 128)).unwrap();
    assert!((result.to_f64().unwrap() - 1.0).abs() < 1e-6);
}

#[test]
fn test_solve_bisection() {
    // Find root of x^2 - 2 = 0 → x = √2
    let f = |x: &Mpf| -> Result<Mpf, mynum::Error> {
        Ok(x.mul(x)
            .sub(&Mpf::from_mpz(Mpz::from_i64(2), x.precision())))
    };
    let a = Mpf::from_f64(1.0, 64);
    let b = Mpf::from_f64(2.0, 64);
    let root = solve::bisection(&f, &a, &b, &Mpf::from_f64(1e-10, 64), 50).unwrap();
    assert!((root.to_f64().unwrap() - 1.41421356).abs() < 1e-6);
}
