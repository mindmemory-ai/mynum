//! 复数指数函数模块
//!
//! 包含复数的指数和对数函数，基于欧拉公式和复数的极坐标表示

use crate::complex::core::Complex;
use crate::error::{Error, Result};
use crate::mpf::Mpf;
use crate::mpz::Mpz;

impl Complex {
    /// 复数的指数函数 e^z
    ///
    /// 使用欧拉公式：e^(a+bi) = e^a * (cos(b) + i*sin(b))
    ///
    /// # 参数
    /// - `self`: 复数 z = a + bi
    ///
    /// # 返回
    /// - `Result<Complex>`: e^z 的结果
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(1, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.exp().unwrap();
    /// // result ≈ e^1 ≈ 2.718...
    /// ```
    pub fn exp(&self) -> Result<Self> {
        let a = self.real();
        let b = self.imaginary();

        // 计算 e^a
        let exp_a = a.exp()?;

        // 如果 b = 0，则结果就是 e^a
        if b.is_zero() {
            return Ok(Self::from_real(exp_a));
        }

        // 计算 cos(b) 和 sin(b)
        let cos_b = b.cos()?;
        let sin_b = b.sin()?;

        // 使用欧拉公式：e^(a+bi) = e^a * (cos(b) + i*sin(b))
        let real_part = exp_a.mul(&cos_b);
        let imag_part = exp_a.mul(&sin_b);

        Ok(Self::from_real_imag(real_part, imag_part))
    }

    /// 复数的自然对数 ln(z)
    ///
    /// 使用公式：ln(z) = ln|z| + i*arg(z)
    /// 其中 |z| 是模长，arg(z) 是幅角
    ///
    /// # 参数
    /// - `self`: 复数 z
    ///
    /// # 返回
    /// - `Result<Complex>`: ln(z) 的结果
    ///
    /// # 注意
    /// - 当 z = 0 时，ln(z) 无定义
    /// - 结果的主值在 (-π, π] 范围内
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(1, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.ln().unwrap();
    /// // result ≈ 0
    /// ```
    pub fn ln(&self) -> Result<Self> {
        if self.is_zero() {
            return Err(Error::DomainError(
                "Natural logarithm of zero is undefined".into(),
            ));
        }

        // 计算模长 |z|
        let magnitude = self.magnitude()?;

        // 计算幅角 arg(z)
        let argument = self.argument()?;

        // ln(z) = ln|z| + i*arg(z)
        let ln_magnitude = magnitude.ln()?;
        let imag_part = Mpf::from_mpz(Mpz::from_i64(1), self.precision()).mul(&argument);

        Ok(Self::from_real_imag(ln_magnitude, imag_part))
    }

    /// 复数的对数函数 log_b(z)
    ///
    /// 使用换底公式：log_b(z) = ln(z) / ln(b)
    ///
    /// # 参数
    /// - `self`: 复数 z
    /// - `base`: 底数 b
    ///
    /// # 返回
    /// - `Result<Complex>`: log_b(z) 的结果
    ///
    /// # 注意
    /// - 当 z = 0 或 b = 0 或 b = 1 时，函数无定义
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(100, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let base = Complex::from_real_imag(
    ///     Mpf::from_i64(10, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.log(&base).unwrap();
    /// // result ≈ 2 (因为 10^2 = 100)
    /// ```
    pub fn log(&self, base: &Complex) -> Result<Self> {
        if self.is_zero() {
            return Err(Error::DomainError("Logarithm of zero is undefined".into()));
        }

        if base.is_zero() || base.is_one() {
            return Err(Error::DomainError("Invalid logarithm base".into()));
        }

        // 计算 ln(z)
        let ln_z = self.ln()?;

        // 计算 ln(base)
        let ln_base = base.ln()?;

        // 使用复数除法：ln(z) / ln(base)
        ln_z.div(&ln_base)
    }

    /// 复数的幂函数 z^w
    ///
    /// 使用公式：z^w = e^(w * ln(z))
    ///
    /// # 参数
    /// - `self`: 底数 z
    /// - `exponent`: 指数 w
    ///
    /// # 返回
    /// - `Result<Complex>`: z^w 的结果
    ///
    /// # 注意
    /// - 当 z = 0 且 w 的实部 ≤ 0 时，函数无定义
    /// - 结果可能不是唯一的（多值函数）
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(2, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let w = Complex::from_real_imag(
    ///     Mpf::from_i64(3, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.pow_complex(&w).unwrap();
    /// // result ≈ 8 (因为 2^3 = 8)
    /// ```
    pub fn pow_complex(&self, exponent: &Complex) -> Result<Self> {
        if self.is_zero() {
            if exponent.real().is_zero() || exponent.real().is_negative() {
                return Err(Error::DomainError(
                    "Zero raised to non-positive power is undefined".into(),
                ));
            }
            return Ok(Self::new());
        }

        // 计算 ln(z)
        let ln_z = self.ln()?;

        // 计算 w * ln(z)
        let product = exponent.mul(&ln_z)?;

        // 计算 e^(w * ln(z))
        product.exp()
    }

    /// 复数的平方根 sqrt(z)
    ///
    /// 使用公式：sqrt(z) = sqrt(|z|) * (cos(θ/2) + i*sin(θ/2))
    /// 其中 |z| 是模长，θ 是幅角
    ///
    /// # 参数
    /// - `self`: 复数 z
    ///
    /// # 返回
    /// - `Result<Complex>`: sqrt(z) 的结果
    ///
    /// # 注意
    /// - 当 z = 0 时，sqrt(0) = 0
    /// - 结果的主值在 (-π/2, π/2] 范围内
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(4, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.sqrt().unwrap();
    /// // result ≈ 2
    /// ```
    pub fn sqrt(&self) -> Result<Self> {
        if self.is_zero() {
            return Ok(Self::new());
        }

        // 计算模长 |z|
        let magnitude = self.magnitude()?;

        // 计算幅角 θ
        let argument = self.argument()?;

        // 计算 sqrt(|z|)
        let sqrt_magnitude = magnitude.sqrt()?;

        // 计算 θ/2
        let half = Mpf::from_mpz(Mpz::from_i64(2), self.precision());
        let half_argument = argument.div(&half)?;

        // 计算 cos(θ/2) 和 sin(θ/2)
        let cos_half = half_argument.cos()?;
        let sin_half = half_argument.sin()?;

        // 使用公式：sqrt(z) = sqrt(|z|) * (cos(θ/2) + i*sin(θ/2))
        let real_part = sqrt_magnitude.mul(&cos_half);
        let imag_part = sqrt_magnitude.mul(&sin_half);

        Ok(Self::from_real_imag(real_part, imag_part))
    }

    /// 复数的立方根 cbrt(z)
    ///
    /// 使用公式：cbrt(z) = z^(1/3) = e^(ln(z)/3)
    ///
    /// # 参数
    /// - `self`: 复数 z
    ///
    /// # 返回
    /// - `Result<Complex>`: cbrt(z) 的结果
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(8, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.cbrt().unwrap();
    /// // result ≈ 2
    /// ```
    pub fn cbrt(&self) -> Result<Self> {
        if self.is_zero() {
            return Ok(Self::new());
        }

        // 计算 ln(z)
        let ln_z = self.ln()?;

        // 计算 ln(z)/3
        let third = Mpf::from_mpz(Mpz::from_i64(3), self.precision());
        let ln_z_third = ln_z.div(&Self::from_real(third))?;

        // 计算 e^(ln(z)/3)
        ln_z_third.exp()
    }

    /// 复数的n次方根 nth_root(z, n)
    ///
    /// 使用公式：nth_root(z, n) = z^(1/n) = e^(ln(z)/n)
    ///
    /// # 参数
    /// - `self`: 复数 z
    /// - `n`: 根数（正整数）
    ///
    /// # 返回
    /// - `Result<Complex>`: z的n次方根
    ///
    /// # 注意
    /// - 当 z = 0 且 n ≤ 0 时，函数无定义
    /// - 当 n = 0 时，函数无定义
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(16, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.nth_root(4).unwrap();
    /// // result ≈ 2 (因为 2^4 = 16)
    /// ```
    pub fn nth_root(&self, n: u32) -> Result<Self> {
        if n == 0 {
            return Err(Error::DomainError("Root order cannot be zero".into()));
        }

        if self.is_zero() {
            return Ok(Self::new());
        }

        // 计算 ln(z)
        let ln_z = self.ln()?;

        // 计算 ln(z)/n
        let n_mpf = Mpf::from_mpz(Mpz::from_u64(n as u64), self.precision());
        let ln_z_nth = ln_z.div(&Self::from_real(n_mpf))?;

        // 计算 e^(ln(z)/n)
        ln_z_nth.exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exp_real() {
        // 测试实数指数：e^0 = 1
        let z = Complex::from_i64(0, 0, 64);
        let result = z.exp().unwrap();

        // 验证结果接近1
        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!((real_part - 1.0).abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_exp_imaginary() {
        // 测试纯虚数指数：e^(i*0) = 1
        let z_pi = Complex::from_i64(0, 0, 64);
        let result = z_pi.exp().unwrap();

        // 验证结果接近1
        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!((real_part - 1.0).abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_ln_real() {
        // 测试实数对数：ln(1) = 0
        let z = Complex::from_i64(1, 0, 64);
        let result = z.ln().unwrap();

        // 验证结果接近0
        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!(real_part.abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_ln_zero() {
        // 测试ln(0)应该返回错误
        let z = Complex::new();
        assert!(z.ln().is_err());
    }

    #[test]
    fn test_sqrt_real() {
        // 测试实数平方根：sqrt(1) = 1
        let z = Complex::from_i64(1, 0, 64);
        let result = z.sqrt().unwrap();

        // 验证结果接近1
        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!((real_part - 1.0).abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_sqrt_negative() {
        // 测试负数平方根：sqrt(-1) = i
        let z = Complex::from_i64(-1, 0, 64);
        let result = z.sqrt().unwrap();

        // 验证结果：sqrt(-1) 应该是一个复数，不是实数
        assert!(!result.is_real());

        // 验证结果不是零
        assert!(!result.is_zero());

        // 验证结果的模长应该接近1
        let magnitude = result.magnitude().unwrap();
        let mag_real = magnitude.to_f64().unwrap();
        assert!((mag_real - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_pow_real() {
        // 测试实数幂：1^1 = 1
        let z = Complex::from_i64(1, 0, 64);
        let w = Complex::from_i64(1, 0, 64);
        let result = z.pow_complex(&w).unwrap();

        // 验证结果接近1
        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!((real_part - 1.0).abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_log_base_10() {
        // 测试以10为底的对数：log_10(100) = 2
        let z = Complex::from_i64(100, 0, 64);
        let base = Complex::from_i64(10, 0, 64);
        let result = z.log(&base).unwrap();

        // 验证结果接近2
        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!((real_part - 2.0).abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_nth_root() {
        // 测试1次方根：1^(1/1) = 1
        let z = Complex::from_i64(1, 0, 64);
        let result = z.nth_root(1).unwrap();

        // 验证结果接近1
        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!((real_part - 1.0).abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_cbrt() {
        // 测试立方根：1^(1/3) = 1
        let z = Complex::from_i64(1, 0, 64);
        let result = z.cbrt().unwrap();

        // 验证结果接近1
        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!((real_part - 1.0).abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }
}
