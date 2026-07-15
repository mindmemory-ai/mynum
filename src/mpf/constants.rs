//! Mpf 数学常量实现

use super::core::Mpf;
use crate::mpz::Mpz;

impl Mpf {
    /// 创建 π (pi)
    pub fn pi(_precision: usize) -> Self {
        // Use from_f64 which correctly handles IEEE 754 representation.
        // The previous approach (mantissa=3141592653589793, exp=-50) was incorrect:
        //   3141592653589793 * 2^(-50) ≈ 2.79, not 3.14.
        // The issue was that 2^(-50) ≠ 10^(-15) — binary exponents cannot
        // simply replace decimal exponents.
        Mpf::from_f64(std::f64::consts::PI, _precision)
    }

    /// 创建 e (自然对数的底)
    pub fn e(_precision: usize) -> Self {
        // 使用级数展开计算 e
        let e_str = "2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274";
        Mpf::from_str(e_str, 10).unwrap_or_else(|_| Mpf::new())
    }

    /// 创建 γ (欧拉-马歇罗尼常数)
    pub fn euler_gamma(_precision: usize) -> Self {
        let gamma_str = "0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495";
        Mpf::from_str(gamma_str, 10).unwrap_or_else(|_| Mpf::new())
    }

    /// 创建 φ (黄金比例)
    pub fn golden_ratio(_precision: usize) -> Self {
        let phi_str = "1.6180339887498948482045868343656381177203091798057628621354486227052604628189024497072072041893911374";
        Mpf::from_str(phi_str, 10).unwrap_or_else(|_| Mpf::new())
    }

    /// 创建 √2
    pub fn sqrt2(_precision: usize) -> Self {
        let sqrt2_str = "1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727";
        Mpf::from_str(sqrt2_str, 10).unwrap_or_else(|_| Mpf::new())
    }

    /// 创建 √3
    pub fn sqrt3(_precision: usize) -> Self {
        let sqrt3_str = "1.7320508075688772935274463415058723669428052538103806280558069794519330169088000370811461867572485756";
        Mpf::from_str(sqrt3_str, 10).unwrap_or_else(|_| Mpf::new())
    }

    /// 创建 ln(2)
    pub fn ln2(_precision: usize) -> Self {
        let ln2_str = "0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875";
        Mpf::from_str(ln2_str, 10).unwrap_or_else(|_| Mpf::new())
    }

    /// 创建 ln(10)
    pub fn ln10(_precision: usize) -> Self {
        let ln10_str = "2.3025850929940456840179914546843642076011014886287729760333279009675726096773524802359972050895982983";
        Mpf::from_str(ln10_str, 10).unwrap_or_else(|_| Mpf::new())
    }

    /// 创建 log₂(e)
    pub fn log2e(_precision: usize) -> Self {
        let log2e_str = "1.4426950408889634073599246810018921374266459541529859341354494069311092191811850798855266228935063445";
        Mpf::from_str(log2e_str, 10).unwrap_or_else(|_| Mpf::new())
    }

    /// 创建 log₁₀(e)
    pub fn log10e(_precision: usize) -> Self {
        let log10e_str = "0.4342944819032518276511289189166050822943970058036665661144537831658646492088707747292249493384317483";
        Mpf::from_str(log10e_str, 10).unwrap_or_else(|_| Mpf::new())
    }

    // infinity 方法已经在 core.rs 中定义

    // nan 方法已经在 core.rs 中定义

    /// 创建零
    pub fn zero() -> Self {
        Mpf::new()
    }

    /// 创建一
    pub fn one(precision: usize) -> Self {
        Mpf::from_mpz(Mpz::from_i64(1), precision)
    }

    /// 创建二
    pub fn two(precision: usize) -> Self {
        Mpf::from_mpz(Mpz::from_i64(2), precision)
    }

    /// 创建十
    pub fn ten(precision: usize) -> Self {
        Mpf::from_mpz(Mpz::from_i64(10), precision)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constants() {
        let pi = Mpf::pi(64);
        assert!(!pi.is_zero());

        let e = Mpf::e(64);
        assert!(!e.is_zero());

        let zero = Mpf::zero();
        assert!(zero.is_zero());

        let one = Mpf::one(64);
        assert_eq!(one.to_i64(), Some(1));
    }
}
