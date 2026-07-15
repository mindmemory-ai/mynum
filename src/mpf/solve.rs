//! Mpf 方程求解模块
//!
//! 提供单变量实值方程的求根算法。

use crate::error::{Error, Result};
use crate::mpf::Mpf;

/// 二分法求根
///
/// 要求 f(a) 和 f(b) 异号。在 [a, b] 中寻找 f(x) = 0 的根。
pub fn bisection<F>(f: &F, a: &Mpf, b: &Mpf, tol: &Mpf, max_iter: usize) -> Result<Mpf>
where
    F: Fn(&Mpf) -> Result<Mpf>,
{
    let fa = f(a)?;
    let fb = f(b)?;

    if fa.is_zero() {
        return Ok(a.clone());
    }
    if fb.is_zero() {
        return Ok(b.clone());
    }
    if fa.is_positive() == fb.is_positive() {
        return Err(Error::domain(
            "bisection requires f(a) and f(b) to have opposite signs",
        ));
    }

    let half = Mpf::from_f64(0.5, a.precision());
    let mut left = a.clone();
    let mut right = b.clone();
    let mut f_left = fa;

    for _ in 0..max_iter {
        let mid = left.add(&right).mul(&half);
        let f_mid = f(&mid)?;

        if f_mid.abs().cmp(tol) == core::cmp::Ordering::Less {
            return Ok(mid);
        }

        if f_left.is_positive() == f_mid.is_positive() {
            left = mid;
            f_left = f_mid;
        } else {
            right = mid;
        }
    }

    Err(Error::convergence(format!(
        "bisection did not converge after {} iterations",
        max_iter
    )))
}

/// 牛顿法求根
///
/// x_{n+1} = x_n - f(x_n) / f'(x_n)
pub fn newton<F, DF>(f: &F, df: &DF, x0: &Mpf, tol: &Mpf, max_iter: usize) -> Result<Mpf>
where
    F: Fn(&Mpf) -> Result<Mpf>,
    DF: Fn(&Mpf) -> Result<Mpf>,
{
    let mut x = x0.clone();

    for _ in 0..max_iter {
        let fx = f(&x)?;
        if fx.abs().cmp(tol) == core::cmp::Ordering::Less {
            return Ok(x);
        }
        let dfx = df(&x)?;
        if dfx.is_zero() {
            return Err(Error::domain("Newton method: derivative is zero"));
        }
        x = x.sub(&fx.div(&dfx)?);
    }

    Err(Error::convergence(format!(
        "Newton method did not converge after {} iterations",
        max_iter
    )))
}

/// 割线法求根
///
/// x_{n+1} = x_n - f(x_n) * (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))
/// 不需要导数。
pub fn secant<F>(f: &F, x0: &Mpf, x1: &Mpf, tol: &Mpf, max_iter: usize) -> Result<Mpf>
where
    F: Fn(&Mpf) -> Result<Mpf>,
{
    let mut x_prev = x0.clone();
    let mut f_prev = f(&x_prev)?;
    let mut x_curr = x1.clone();

    if f_prev.abs().cmp(tol) == core::cmp::Ordering::Less {
        return Ok(x_prev);
    }

    for _ in 0..max_iter {
        let f_curr = f(&x_curr)?;
        if f_curr.abs().cmp(tol) == core::cmp::Ordering::Less {
            return Ok(x_curr);
        }

        let denom = f_curr.sub(&f_prev);
        if denom.is_zero() {
            return Err(Error::domain("secant method: division by zero"));
        }

        let x_next = x_curr.sub(&f_curr.mul(&x_curr.sub(&x_prev)).div(&denom)?);
        x_prev = x_curr;
        f_prev = f_curr;
        x_curr = x_next;
    }

    Err(Error::convergence(format!(
        "secant method did not converge after {} iterations",
        max_iter
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mpf::Mpf;

    #[test]
    fn test_bisection_quadratic() {
        let f = |x: &Mpf| Ok(x.mul(x).sub(&Mpf::from_i64(2, 64)));
        let a = Mpf::from_i64(1, 64);
        let b = Mpf::from_i64(2, 64);
        let tol = Mpf::from_f64(1e-10, 64);
        let result = bisection(&f, &a, &b, &tol, 100).unwrap();
        assert!((result.to_f64().unwrap() - 1.41421356).abs() < 1e-8);
    }

    #[test]
    fn test_newton_sqrt2() {
        let f = |x: &Mpf| Ok(x.mul(x).sub(&Mpf::from_i64(2, 64)));
        let df = |x: &Mpf| Ok(x.mul(&Mpf::from_i64(2, 64)));
        let x0 = Mpf::from_f64(1.5, 64);
        let tol = Mpf::from_f64(1e-10, 64);
        let result = newton(&f, &df, &x0, &tol, 50).unwrap();
        assert!((result.to_f64().unwrap() - 1.41421356).abs() < 1e-8);
    }

    #[test]
    fn test_secant_cubic() {
        let f = |x: &Mpf| {
            let x3 = x.mul(x).mul(x);
            Ok(x3.sub(x).sub(&Mpf::from_i64(2, 64)))
        };
        let x0 = Mpf::from_i64(1, 64);
        let x1 = Mpf::from_i64(2, 64);
        let tol = Mpf::from_f64(1e-10, 64);
        let result = secant(&f, &x0, &x1, &tol, 50).unwrap();
        assert!((result.to_f64().unwrap() - 1.5213797).abs() < 1e-6);
    }
}
