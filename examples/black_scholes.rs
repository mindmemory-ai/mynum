//! Black-Scholes 期权定价
//!
//! 使用高精度 Mpf 计算欧式看涨期权价格。

use mynum::mpf::special::erf;
use mynum::mpf::Mpf;

fn normal_cdf(x: &Mpf) -> Mpf {
    let prec = x.precision();
    let half = Mpf::from_f64(0.5, prec);
    let sqrt2 = Mpf::from_f64(std::f64::consts::SQRT_2, prec);
    let erf_arg = x.div(&sqrt2).unwrap();
    let erf_val = erf(&erf_arg).unwrap();
    half.add(&half.mul(&erf_val))
}

fn black_scholes_call(s: &Mpf, k: &Mpf, t: &Mpf, r: &Mpf, sigma: &Mpf) -> Mpf {
    let prec = s.precision();
    let sqrt_t = Mpf::sqrt(t).unwrap();
    let sigma_sqrt_t = sigma.mul(&sqrt_t);
    let half_sigma2 = sigma.mul(sigma).div(&Mpf::from_f64(2.0, prec)).unwrap();

    let d1 = s
        .div(k)
        .unwrap()
        .ln()
        .unwrap()
        .add(&r.add(&half_sigma2).mul(t))
        .div(&sigma_sqrt_t)
        .unwrap();
    let d2 = d1.sub(&sigma_sqrt_t);

    let discount = r.mul(t).neg().exp().unwrap();
    let k_discounted = k.mul(&discount);

    s.mul(&normal_cdf(&d1))
        .sub(&k_discounted.mul(&normal_cdf(&d2)))
}

fn main() {
    let prec = 64;
    let s = Mpf::from_f64(100.0, prec);
    let k = Mpf::from_f64(105.0, prec);
    let t = Mpf::from_f64(1.0, prec);
    let r = Mpf::from_f64(0.05, prec);
    let sigma = Mpf::from_f64(0.2, prec);

    let call = black_scholes_call(&s, &k, &t, &r, &sigma);
    println!("Black-Scholes 欧式看涨期权: {:.6}", call);
    // 预期: S=100, K=105, T=1, r=5%, σ=20% → C ≈ 8.02
}
