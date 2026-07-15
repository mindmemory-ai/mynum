//! Special function precision tests — Gamma, Bessel, Zeta, erf vs known reference values.
//! Reference values from DLMF (NIST Digital Library of Mathematical Functions).

use mynum::mpf::special::*;
use mynum::mpf::Mpf;

#[test]
fn test_gamma_known_values() {
    let cases: [(f64, f64); 8] = [
        (0.5, 1.7724538509055160273),
        (1.0, 1.0),
        (2.0, 1.0),
        (3.0, 2.0),
        (5.0, 24.0),
        (10.0, 362880.0),
        (0.25, 3.6256099082219083119),
        (1.5, 0.8862269254527580137),
    ];
    for (x, expected) in &cases {
        let val = gamma(&Mpf::from_f64(*x, 128)).unwrap();
        let v = val.to_f64().unwrap();
        assert!(
            (v - expected).abs() < 1e-12,
            "Gamma({}) = {}, expected {}",
            x,
            v,
            expected
        );
    }
}

#[test]
fn test_gamma_reflection() {
    for x in &[0.1, 0.25, 0.33, 0.5, 0.67, 0.75, 0.9] {
        let gx = gamma(&Mpf::from_f64(*x, 128)).unwrap();
        let g1x = gamma(&Mpf::from_f64(1.0 - x, 128)).unwrap();
        let product = gx.mul(&g1x).to_f64().unwrap();
        let expected = std::f64::consts::PI / (std::f64::consts::PI * x).sin();
        assert!(
            (product - expected).abs() < 1e-11,
            "Gamma reflection failed at x={}: got {}",
            x,
            product
        );
    }
}

#[test]
fn test_gamma_poles() {
    assert!(gamma(&Mpf::new()).is_err());
    assert!(gamma(&Mpf::from_i64(-1, 64)).is_err());
    assert!(gamma(&Mpf::from_i64(-2, 64)).is_err());
}

#[test]
fn test_bessel_jn_vs_dlmf() {
    let cases = [
        (0, 1.0, 0.7651976865579666),
        (1, 1.0, 0.4400505857449335),
        (2, 1.0, 0.1149034849319005),
        (0, 5.0, -0.1775967713143383),
        (3, 2.0, 0.12894324947440205),
        (5, 10.0, 0.23406152818679366),
        (10, 1.0, 2.630615123687453e-17),
    ];
    for (n, x, expected) in &cases {
        let val = bessel_jn(*n, &Mpf::from_f64(*x, 128)).unwrap();
        let v = val.to_f64().unwrap();
        let tol = if *expected == 0.0 { 1e-15 } else { 1e-12 };
        assert!(
            (v - expected).abs() < tol,
            "J_{}({}) = {}, expected {}",
            n,
            x,
            v,
            expected
        );
    }
}

#[test]
fn test_bessel_yn_vs_dlmf() {
    let cases = [
        (0, 1.0, 0.08825696421567698),
        (1, 1.0, -0.7812128213002887),
        (0, 5.0, -0.30851762524903377),
        (2, 0.5, -5.441370837174266),
    ];
    for (n, x, expected) in &cases {
        let val = bessel_yn(*n, &Mpf::from_f64(*x, 128)).unwrap();
        let v = val.to_f64().unwrap();
        assert!(
            (v - expected).abs() < 1e-10,
            "Y_{}({}) = {}, expected {}",
            n,
            x,
            v,
            expected
        );
    }
}

#[test]
fn test_bessel_k0_k1() {
    // K_0(1) ≈ 0.42102443824070834
    let x = Mpf::from_f64(1.0, 128);
    let k0 = bessel_k0(&x).unwrap();
    assert!((k0.to_f64().unwrap() - 0.4210244382).abs() < 1e-8);
    // K_1(1) ≈ 0.6019072301972346
    let k1 = bessel_k1(&x).unwrap();
    assert!((k1.to_f64().unwrap() - 0.6019072302).abs() < 1e-8);
}

#[test]
fn test_bessel_i0_i1() {
    let x = Mpf::from_f64(1.0, 128);
    let i0 = bessel_i0(&x).unwrap();
    assert!((i0.to_f64().unwrap() - 1.2660658780).abs() < 1e-8);
    let i1 = bessel_i1(&x).unwrap();
    assert!((i1.to_f64().unwrap() - 0.5651591040).abs() < 1e-8);
}

#[test]
fn test_digamma_known() {
    let cases = [
        (1.0, -0.5772156649015329),
        (2.0, 0.42278433509846713),
        (10.0, 2.2517525890667211),
        (0.5, -1.9635100260214235),
    ];
    for (x, expected) in &cases {
        let val = digamma(&Mpf::from_f64(*x, 128)).unwrap();
        assert!(
            (val.to_f64().unwrap() - expected).abs() < 1e-12,
            "digamma({})",
            x
        );
    }
}

#[test]
fn test_zeta_known() {
    let cases = [
        (2.0, 1.6449340668482264),
        (3.0, 1.2020569031595942),
        (4.0, 1.0823232337111382),
    ];
    for (s, expected) in &cases {
        let val = zeta(&Mpf::from_f64(*s, 128)).unwrap();
        assert!(
            (val.to_f64().unwrap() - expected).abs() < 1e-12,
            "zeta({})",
            s
        );
    }
}

#[test]
fn test_erf_known() {
    let cases = [
        (0.0, 0.0),
        (0.5, 0.5204998778130465),
        (1.0, 0.8427007929497149),
        (2.0, 0.9953222650189527),
        (-1.0, -0.8427007929497149),
    ];
    for (x, expected) in &cases {
        let val = erf(&Mpf::from_f64(*x, 128)).unwrap();
        assert!(
            (val.to_f64().unwrap() - expected).abs() < 1e-12,
            "erf({})",
            x
        );
    }
}

#[test]
fn test_elliptic_k_known() {
    // K(0) = π/2, K(0.5) ≈ 1.8540746773013719
    let m_zero = Mpf::new();
    let k0 = elliptic_k(&m_zero).unwrap();
    let pi_half = std::f64::consts::PI / 2.0;
    assert!((k0.to_f64().unwrap() - pi_half).abs() < 1e-10);
    let m_half = Mpf::from_f64(0.5, 128);
    let k_half = elliptic_k(&m_half).unwrap();
    assert!((k_half.to_f64().unwrap() - 1.8540746773).abs() < 1e-9);
}

#[test]
fn test_hypergeometric_2f1() {
    // ₂F₁(1,1;2;0.5) = 2*ln(2) ≈ 1.38629436111989
    let a = Mpf::from_i64(1, 128);
    let b = Mpf::from_i64(1, 128);
    let c = Mpf::from_i64(2, 128);
    let z = Mpf::from_f64(0.5, 128);
    let result = hypergeometric_2f1(&a, &b, &c, &z).unwrap();
    assert!((result.to_f64().unwrap() - 1.3862943611).abs() < 1e-8);
}

#[test]
fn test_lambert_w() {
    let one = Mpf::from_i64(1, 64);
    let w = lambert_w(&one, 0).unwrap();
    assert!((w.to_f64().unwrap() - 0.5671432904).abs() < 1e-8);
    let neg = Mpf::from_f64(-0.2, 64);
    let w_neg = lambert_w(&neg, -1).unwrap();
    assert!((w_neg.to_f64().unwrap() - (-2.5426413578)).abs() < 1e-6);
}

#[test]
fn test_expint_ei_known() {
    let cases = [
        (-1.0, -0.21938393439552),
        (1.0, 1.89511781635594),
        (0.5, 0.4542199048631736),
    ];
    for (x, expected) in &cases {
        let val = expint_ei(&Mpf::from_f64(*x, 128)).unwrap();
        assert!(
            (val.to_f64().unwrap() - expected).abs() < 1e-10,
            "Ei({}) = {}, expected {}",
            x,
            val.to_f64().unwrap(),
            expected
        );
    }
}

#[test]
fn test_expint_ei_large() {
    let x = Mpf::from_f64(40.0, 256);
    let val = expint_ei(&x).unwrap();
    let v = val.to_f64().unwrap();
    assert!(v > 1e14, "Ei(40) should be large, got {}", v);
}

#[test]
fn test_expint_ei_zero() {
    assert!(expint_ei(&Mpf::new()).is_err());
}

#[test]
fn test_erfi_known() {
    let cases = [
        (0.5, 0.614952094696511),
        (1.0, 1.650425758797543),
        (2.0, 18.56480241457555),
    ];
    for (x, expected) in &cases {
        let val = erfi(&Mpf::from_f64(*x, 128)).unwrap();
        assert!(
            (val.to_f64().unwrap() - expected).abs() < 1e-10,
            "erfi({}) = {}, expected {}",
            x,
            val.to_f64().unwrap(),
            expected
        );
    }
}

#[test]
fn test_erfi_symmetry() {
    let x = Mpf::from_f64(1.5, 128);
    let v1 = erfi(&x).unwrap();
    let v2 = erfi(&x.neg()).unwrap();
    assert!((v1.to_f64().unwrap() + v2.to_f64().unwrap()).abs() < 1e-12);
}

#[test]
fn test_erfi_zero() {
    assert!(erfi(&Mpf::new()).unwrap().is_zero());
}

#[test]
fn test_sinint_si_known() {
    let cases = [
        (0.0, 0.0),
        (1.0, 0.9460830703671830),
        (std::f64::consts::PI, 1.851937051982466),
        (10.0, 1.658347594218874),
    ];
    for (x, expected) in &cases {
        let val = sinint_si(&Mpf::from_f64(*x, 128)).unwrap();
        let tol = if *x == 0.0 { 1e-15 } else { 1e-10 };
        assert!(
            (val.to_f64().unwrap() - expected).abs() < tol,
            "Si({}) = {}, expected {}",
            x,
            val.to_f64().unwrap(),
            expected
        );
    }
}

#[test]
fn test_cosint_ci_known() {
    let cases = [
        (1.0, 0.33740392290096813),
        (std::f64::consts::PI, 0.07366791204645162),
        (5.0, -0.19002974965664395),
    ];
    for (x, expected) in &cases {
        let val = cosint_ci(&Mpf::from_f64(*x, 128)).unwrap();
        assert!(
            (val.to_f64().unwrap() - expected).abs() < 1e-10,
            "Ci({}) = {}, expected {}",
            x,
            val.to_f64().unwrap(),
            expected
        );
    }
}

#[test]
fn test_cosint_ci_domain() {
    assert!(cosint_ci(&Mpf::new()).is_err()); // x ≤ 0
}

#[test]
fn test_sinint_si_negative() {
    // Si(-x) = -Si(x) (odd function)
    let x = Mpf::from_f64(2.5, 128);
    let pos = sinint_si(&x).unwrap();
    let neg = sinint_si(&x.neg()).unwrap();
    assert!(
        (pos.to_f64().unwrap() + neg.to_f64().unwrap()).abs() < 1e-12,
        "Si({}) and Si({}) should be opposites",
        x.to_f64().unwrap(),
        (-x.to_f64().unwrap())
    );
}

// ── Fresnel integrals ──

#[test]
fn test_fresnel_s_known() {
    let cases = [
        (0.5, 0.06473243285999928),
        (1.0, 0.43825914739035477),
        (5.0, 0.4991913819171169),
    ];
    for (x, expected) in &cases {
        let val = fresnel_s(&Mpf::from_f64(*x, 128)).unwrap();
        assert!(
            (val.to_f64().unwrap() - expected).abs() < 1e-10,
            "S({}) = {}, expected {}",
            x,
            val.to_f64().unwrap(),
            expected
        );
    }
}

#[test]
fn test_fresnel_c_known() {
    let cases = [
        (0.5, 0.4923442258714464),
        (1.0, 0.7798934003768228),
        (5.0, 0.5636311887040122),
    ];
    for (x, expected) in &cases {
        let val = fresnel_c(&Mpf::from_f64(*x, 128)).unwrap();
        assert!(
            (val.to_f64().unwrap() - expected).abs() < 1e-10,
            "C({}) = {}, expected {}",
            x,
            val.to_f64().unwrap(),
            expected
        );
    }
}

#[test]
fn test_fresnel_odd() {
    let x = Mpf::from_f64(1.0, 128);
    assert!(
        (fresnel_s(&x.neg()).unwrap().to_f64().unwrap() + fresnel_s(&x).unwrap().to_f64().unwrap())
            .abs()
            < 1e-12,
        "S(x) should be odd"
    );
    assert!(
        (fresnel_c(&x.neg()).unwrap().to_f64().unwrap() + fresnel_c(&x).unwrap().to_f64().unwrap())
            .abs()
            < 1e-12,
        "C(x) should be odd"
    );
}

#[test]
fn test_fresnel_zero() {
    let zero = Mpf::new();
    assert!(fresnel_s(&zero).unwrap().is_zero());
    assert!(fresnel_c(&zero).unwrap().is_zero());
}

// ── Struve H ──

#[test]
fn test_struve_h_known() {
    // H_0(1) ≈ 0.568656627048, H_1(1) ≈ 0.1984573362, H_0(10) ≈ 0.118743
    let x1 = Mpf::from_f64(1.0, 128);
    assert!((struve_h(0, &x1).unwrap().to_f64().unwrap() - 0.568656627048).abs() < 1e-8);
    assert!((struve_h(1, &x1).unwrap().to_f64().unwrap() - 0.1984573362).abs() < 1e-8);
    let x10 = Mpf::from_f64(10.0, 128);
    assert!((struve_h(0, &x10).unwrap().to_f64().unwrap() - 0.118743).abs() < 1e-5);
}

// ── Lommel s ──

#[test]
fn test_lommel_s_known() {
    // s_{μ,ν}(x) = x^{μ+1}/((μ+1)²-ν²)·₁F₂(1; (μ-ν+3)/2, (μ+ν+3)/2; -x²/4)
    // s_{1,0}(1) = 1/4 * ₁F₂(1; 2, 2; -0.25) ≈ 0.2348023
    let mu = Mpf::from_i64(1, 128);
    let x = Mpf::from_f64(1.0, 128);
    let val = lommel_s(&mu, 0, &x).unwrap();
    assert!((val.to_f64().unwrap() - 0.2348023).abs() < 1e-5);
}

// ── Whittaker M ──

#[test]
fn test_whittaker_m_known() {
    // M_{κ,μ}(x) = e^{-x/2}·x^{μ+1/2}·₁F₁(1/2+μ-κ; 1+2μ; x)
    // M_{0,0}(2) = e^{-1}·√2·₁F₁(1/2;1;2) ≈ 1.79049
    let kappa = Mpf::new();
    let mu = Mpf::new();
    let x = Mpf::from_f64(2.0, 128);
    let val = whittaker_m(&kappa, &mu, &x).unwrap();
    assert!((val.to_f64().unwrap() - 1.79049).abs() < 1e-4);
}

#[test]
fn test_whittaker_m_domain() {
    // μ = -0.5: 1+2μ = 0, should error
    let kappa = Mpf::new();
    let mu = Mpf::from_f64(-0.5, 64);
    let x = Mpf::from_f64(1.0, 64);
    assert!(whittaker_m(&kappa, &mu, &x).is_err());

    // μ = -1: 1+2μ = -1, should error
    let mu2 = Mpf::from_i64(-1, 64);
    assert!(whittaker_m(&kappa, &mu2, &x).is_err());
}

// ── Hyperbolic sine/cosine integrals Shi(x), Chi(x) ──

#[test]
fn test_sinhint_shi_known() {
    let cases = [
        (0.5, 0.5069967498196672),
        (1.0, 1.0572508753757285),
        (2.0, 2.501567433354976),
    ];
    for (x, expected) in &cases {
        let val = sinhint_shi(&Mpf::from_f64(*x, 128)).unwrap();
        let actual = val.to_f64().unwrap();
        assert!((actual - expected).abs() < 1e-10);
    }
}

#[test]
fn test_coshint_chi_known() {
    let cases = [(1.0, 0.8378669409802082), (2.0, 2.452666922646914)];
    for (x, expected) in &cases {
        let val = coshint_chi(&Mpf::from_f64(*x, 128)).unwrap();
        let actual = val.to_f64().unwrap();
        assert!((actual - expected).abs() < 1e-10);
    }
}

#[test]
fn test_coshint_chi_domain() {
    assert!(coshint_chi(&Mpf::new()).is_err());
}

// ── Real-order Bessel functions ──

#[test]
fn test_bessel_j_real_half() {
    // J_{0.5}(1) = √(2/π)·sin(1) ≈ 0.6713967071
    let nu = Mpf::from_f64(0.5, 128);
    let x = Mpf::from_f64(1.0, 128);
    let val = bessel_j_real(&nu, &x).unwrap();
    assert!((val.to_f64().unwrap() - 0.6713967071).abs() < 1e-8);
}

#[test]
fn test_bessel_j_real_match_integer() {
    // J_{2.0} should match bessel_jn(2, x)
    let nu = Mpf::from_f64(2.0, 128);
    let x = Mpf::from_f64(3.0, 128);
    let val = bessel_j_real(&nu, &x).unwrap();
    let expected = bessel_jn(2, &x).unwrap();
    assert!((val.to_f64().unwrap() - expected.to_f64().unwrap()).abs() < 1e-8);
}

#[test]
fn test_bessel_j_real_negative_int() {
    // J_{-1}(x) = -J_1(x)
    let nu = Mpf::from_i64(-1, 128);
    let x = Mpf::from_f64(1.0, 128);
    let val = bessel_j_real(&nu, &x).unwrap();
    let expected = bessel_jn(1, &x).unwrap().neg();
    assert!((val.to_f64().unwrap() - expected.to_f64().unwrap()).abs() < 1e-12);
}

#[test]
fn test_bessel_y_real_half() {
    // Y_{0.5}(1) ≈ -0.4310988689
    let nu = Mpf::from_f64(0.5, 128);
    let x = Mpf::from_f64(1.0, 128);
    let val = bessel_y_real(&nu, &x).unwrap();
    assert!((val.to_f64().unwrap() - (-0.4310988689)).abs() < 1e-8);
}

#[test]
fn test_bessel_i_real_half() {
    // I_{0.5}(1) = √(2/π)·sinh(1) ≈ 0.9376748882
    let nu = Mpf::from_f64(0.5, 128);
    let x = Mpf::from_f64(1.0, 128);
    let val = bessel_i_real(&nu, &x).unwrap();
    assert!((val.to_f64().unwrap() - 0.9376748882).abs() < 1e-7);
}

#[test]
fn test_bessel_k_real_half() {
    // K_{0.5}(x) = √(π/(2x))·e^{-x}
    // K_{0.5}(1) = √(π/2)·e^{-1} ≈ 0.4610685073
    let nu = Mpf::from_f64(0.5, 128);
    let x = Mpf::from_f64(1.0, 128);
    let val = bessel_k_real(&nu, &x).unwrap();
    assert!((val.to_f64().unwrap() - 0.4610685073).abs() < 1e-7);
}

// ── Elliptic integrals ──

#[test]
fn test_elliptic_pi_known() {
    // Π(0.5, 0.5) ≈ 2.7012877621 (verified via Gauss-Legendre & trapezoidal)
    let n = Mpf::from_f64(0.5, 64);
    let m = Mpf::from_f64(0.5, 64);
    let val = elliptic_pi(&n, &m).unwrap();
    assert!((val.to_f64().unwrap() - 2.70128776).abs() < 1e-5);
}

#[test]
fn test_elliptic_pi_zero_n() {
    // Π(0, m) = K(m)
    let n = Mpf::new(); // n = 0
    let m = Mpf::from_f64(0.5, 128);
    let pi_val = elliptic_pi(&n, &m).unwrap();
    let k_val = elliptic_k(&m).unwrap();
    assert!((pi_val.to_f64().unwrap() - k_val.to_f64().unwrap()).abs() < 1e-8);
}

// ── Jacobi elliptic functions ──

#[test]
fn test_jacobi_sn_known() {
    // sn(0.5, 0.5) ≈ 0.4707504736555765
    let u = Mpf::from_f64(0.5, 64);
    let m = Mpf::from_f64(0.5, 64);
    let val = jacobi_sn(&u, &m).unwrap();
    assert!((val.to_f64().unwrap() - 0.4707504736).abs() < 1e-5);
}

#[test]
fn test_jacobi_sn_special_cases() {
    // sn(u, 0) = sin(u)
    let u = Mpf::from_f64(0.8, 64);
    let m_zero = Mpf::new();
    let sn = jacobi_sn(&u, &m_zero).unwrap();
    assert!((sn.to_f64().unwrap() - 0.8_f64.sin()).abs() < 1e-9);
    // sn(0, m) = 0
    let zero = Mpf::new();
    let m = Mpf::from_f64(0.3, 64);
    assert!(jacobi_sn(&zero, &m).unwrap().is_zero());
}

#[test]
fn test_jacobi_cn_identity() {
    // sn² + cn² = 1
    let u = Mpf::from_f64(0.7, 64);
    let m = Mpf::from_f64(0.3, 64);
    let s = jacobi_sn(&u, &m).unwrap();
    let c = jacobi_cn(&u, &m).unwrap();
    let sum = s.mul(&s).add(&c.mul(&c));
    assert!((sum.to_f64().unwrap() - 1.0).abs() < 1e-6);
}

#[test]
fn test_jacobi_dn_identity() {
    // m·sn² + dn² = 1
    let m = Mpf::from_f64(0.3, 64);
    let u = Mpf::from_f64(0.7, 64);
    let s = jacobi_sn(&u, &m).unwrap();
    let d = jacobi_dn(&u, &m).unwrap();
    let sum = m.mul(&s.mul(&s)).add(&d.mul(&d));
    assert!((sum.to_f64().unwrap() - 1.0).abs() < 1e-6);
}

#[test]
fn test_jacobi_cn_special_cases() {
    // cn(u, 0) = cos(u)
    let u = Mpf::from_f64(0.8, 64);
    let m_zero = Mpf::new();
    let cn = jacobi_cn(&u, &m_zero).unwrap();
    assert!((cn.to_f64().unwrap() - 0.8_f64.cos()).abs() < 1e-9);
    // cn(0, m) = 1
    let zero = Mpf::new();
    let m = Mpf::from_f64(0.3, 64);
    assert!((jacobi_cn(&zero, &m).unwrap().to_f64().unwrap() - 1.0).abs() < 1e-12);
}

#[test]
fn test_jacobi_dn_special_cases() {
    // dn(u, 0) = 1
    let u = Mpf::from_f64(0.8, 64);
    let m_zero = Mpf::new();
    let dn = jacobi_dn(&u, &m_zero).unwrap();
    assert!((dn.to_f64().unwrap() - 1.0).abs() < 1e-12);
    // dn(0, m) = 1
    let zero = Mpf::new();
    let m = Mpf::from_f64(0.3, 64);
    assert!((jacobi_dn(&zero, &m).unwrap().to_f64().unwrap() - 1.0).abs() < 1e-12);
}
