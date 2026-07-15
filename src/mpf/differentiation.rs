//! Mpf 数值微分模块

use crate::error::Result;
use crate::mpf::Mpf;
use crate::mpz::Mpz;

/// 前向差分: f'(x) ≈ (f(x+h) - f(x)) / h
/// 一阶精度 O(h)
pub fn forward_diff<F>(f: &F, x: &Mpf, h: &Mpf) -> Result<Mpf>
where
    F: Fn(&Mpf) -> Result<Mpf>,
{
    let fxh = f(&x.add(h))?;
    let fx = f(x)?;
    fxh.sub(&fx).div(h)
}

/// 中心差分: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
/// 二阶精度 O(h^2)
pub fn central_diff<F>(f: &F, x: &Mpf, h: &Mpf) -> Result<Mpf>
where
    F: Fn(&Mpf) -> Result<Mpf>,
{
    let fxh = f(&x.add(h))?;
    let fx_mh = f(&x.sub(h))?;
    let two_h = h.mul(&Mpf::from_mpz(Mpz::from_i64(2), h.precision()));
    fxh.sub(&fx_mh).div(&two_h)
}

/// 使用 Richardson 外推提高导数精度
///
/// 从步长 h 开始，反复减半步长并外推以消除误差项。
/// `levels` 控制外推深度（通常 3-5 足够）。
pub fn richardson_extrapolation<F>(f: &F, x: &Mpf, h: &Mpf, levels: usize) -> Result<Mpf>
where
    F: Fn(&Mpf) -> Result<Mpf>,
{
    if levels == 0 {
        return central_diff(f, x, h);
    }

    let two = Mpf::from_mpz(Mpz::from_i64(2), h.precision());
    let mut table: Vec<Vec<Mpf>> = Vec::with_capacity(levels + 1);

    for i in 0..=levels {
        let mut pow = Mpf::from_mpz(Mpz::from_i64(1), h.precision());
        for _ in 0..i {
            pow = pow.mul(&two);
        }
        let hi = h.div(&pow)?;
        table.push(vec![central_diff(f, x, &hi)?]);
    }

    for j in 1..=levels {
        let mut four_pow = Mpf::from_mpz(Mpz::from_i64(4), h.precision());
        for _ in 1..j {
            four_pow = four_pow.mul(&Mpf::from_mpz(Mpz::from_i64(4), h.precision()));
        }
        let one = Mpf::from_mpz(Mpz::from_i64(1), h.precision());

        for i in 0..=levels - j {
            let num = table[i + 1][j - 1].mul(&four_pow).sub(&table[i][j - 1]);
            let denom = four_pow.sub(&one);
            table[i].push(num.div(&denom)?);
        }
    }

    Ok(table[0][levels].clone())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_forward_diff_exp() {
        let f = |x: &Mpf| x.exp();
        let x = Mpf::new();
        let h = Mpf::from_f64(1e-5, 64);
        let result = forward_diff(&f, &x, &h).unwrap();
        assert!((result.to_f64().unwrap() - 1.0).abs() < 1e-4);
    }

    #[test]
    fn test_central_diff_sin() {
        let f = |x: &Mpf| x.sin();
        let x = Mpf::new();
        let h = Mpf::from_f64(1e-5, 64);
        let result = central_diff(&f, &x, &h).unwrap();
        assert!((result.to_f64().unwrap() - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_richardson_accuracy() {
        let f = |x: &Mpf| x.exp();
        let x = Mpf::new();
        let h = Mpf::from_f64(0.1, 64);
        let central = central_diff(&f, &x, &h).unwrap();
        let rich = richardson_extrapolation(&f, &x, &h, 3).unwrap();
        let one = Mpf::from_i64(1, 64);
        let central_err = central.sub(&one).abs();
        let rich_err = rich.sub(&one).abs();
        assert!(rich_err.to_f64().unwrap() < central_err.to_f64().unwrap());
    }
}
