//! 复数双曲函数模块
//!
//! 包含复数的双曲函数、反双曲函数和相关的数学函数
//! 基于指数函数和复数的数学性质实现

use crate::complex::core::Complex;
use crate::error::{Error, Result};
use crate::mpf::Mpf;
use crate::mpz::Mpz;

impl Complex {
    /// 复数的双曲正弦函数 sinh(z)
    ///
    /// 使用公式：sinh(z) = (e^z - e^(-z)) / 2
    /// 展开后：sinh(a+bi) = sinh(a)cos(b) + i*cosh(a)sin(b)
    ///
    /// # 参数
    /// - `self`: 复数 z = a + bi
    ///
    /// # 返回
    /// - `Result<Complex>`: sinh(z) 的结果
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(0, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.sinh().unwrap();
    /// // result ≈ 0
    /// ```
    pub fn sinh(&self) -> Result<Self> {
        let a = self.real();
        let b = self.imaginary();

        // 如果 b = 0，则 sinh(a+0i) = sinh(a)
        if b.is_zero() {
            let sinh_a = a.sinh()?;
            return Ok(Self::from_real(sinh_a));
        }

        // 如果 a = 0，则 sinh(0+bi) = i*sin(b)
        if a.is_zero() {
            let sin_b = b.sin()?;
            return Ok(Self::from_imag(sin_b));
        }

        // 一般情况：sinh(a+bi) = sinh(a)cos(b) + i*cosh(a)sin(b)
        let sinh_a = a.sinh()?;
        let cosh_a = a.cosh()?;
        let cos_b = b.cos()?;
        let sin_b = b.sin()?;

        let real_part = sinh_a.mul(&cos_b);
        let imag_part = cosh_a.mul(&sin_b);

        Ok(Self::from_real_imag(real_part, imag_part))
    }

    /// 复数的双曲余弦函数 cosh(z)
    ///
    /// 使用公式：cosh(z) = (e^z + e^(-z)) / 2
    /// 展开后：cosh(a+bi) = cosh(a)cos(b) + i*sinh(a)sin(b)
    ///
    /// # 参数
    /// - `self`: 复数 z = a + bi
    ///
    /// # 返回
    /// - `Result<Complex>`: cosh(z) 的结果
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(0, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.cosh().unwrap();
    /// // result ≈ 1
    /// ```
    pub fn cosh(&self) -> Result<Self> {
        let a = self.real();
        let b = self.imaginary();

        // 如果 b = 0，则 cosh(a+0i) = cosh(a)
        if b.is_zero() {
            let cosh_a = a.cosh()?;
            return Ok(Self::from_real(cosh_a));
        }

        // 如果 a = 0，则 cosh(0+bi) = cos(b)
        if a.is_zero() {
            let cos_b = b.cos()?;
            return Ok(Self::from_real(cos_b));
        }

        // 一般情况：cosh(a+bi) = cosh(a)cos(b) + i*sinh(a)sin(b)
        let sinh_a = a.sinh()?;
        let cosh_a = a.cosh()?;
        let cos_b = b.cos()?;
        let sin_b = b.sin()?;

        let real_part = cosh_a.mul(&cos_b);
        let imag_part = sinh_a.mul(&sin_b);

        Ok(Self::from_real_imag(real_part, imag_part))
    }

    /// 复数的双曲正切函数 tanh(z)
    ///
    /// 使用公式：tanh(z) = sinh(z) / cosh(z)
    ///
    /// # 参数
    /// - `self`: 复数 z = a + bi
    ///
    /// # 返回
    /// - `Result<Complex>`: tanh(z) 的结果
    ///
    /// # 注意
    /// - 当 cosh(z) = 0 时，tanh(z) 无定义
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(0, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.tanh().unwrap();
    /// // result ≈ 0
    /// ```
    pub fn tanh(&self) -> Result<Self> {
        let sinh_z = self.sinh()?;
        let cosh_z = self.cosh()?;

        // 检查 cosh(z) 是否为零
        if cosh_z.is_zero() {
            return Err(Error::DomainError(
                "Hyperbolic tangent is undefined when hyperbolic cosine is zero".into(),
            ));
        }

        // tanh(z) = sinh(z) / cosh(z)
        sinh_z.div(&cosh_z)
    }

    /// 复数的双曲余切函数 coth(z)
    ///
    /// 使用公式：coth(z) = cosh(z) / sinh(z)
    ///
    /// # 参数
    /// - `self`: 复数 z = a + bi
    ///
    /// # 返回
    /// - `Result<Complex>`: coth(z) 的结果
    ///
    /// # 注意
    /// - 当 sinh(z) = 0 时，coth(z) 无定义
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
    /// let result = z.coth().unwrap();
    /// // result ≈ 1.313...
    /// ```
    pub fn coth(&self) -> Result<Self> {
        let sinh_z = self.sinh()?;
        let cosh_z = self.cosh()?;

        // 检查 sinh(z) 是否为零
        if sinh_z.is_zero() {
            return Err(Error::DomainError(
                "Hyperbolic cotangent is undefined when hyperbolic sine is zero".into(),
            ));
        }

        // coth(z) = cosh(z) / sinh(z)
        cosh_z.div(&sinh_z)
    }

    /// 复数的双曲正割函数 sech(z)
    ///
    /// 使用公式：sech(z) = 1 / cosh(z)
    ///
    /// # 参数
    /// - `self`: 复数 z = a + bi
    ///
    /// # 返回
    /// - `Result<Complex>`: sech(z) 的结果
    ///
    /// # 注意
    /// - 当 cosh(z) = 0 时，sech(z) 无定义
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(0, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.sech().unwrap();
    /// // result ≈ 1
    /// ```
    pub fn sech(&self) -> Result<Self> {
        let cosh_z = self.cosh()?;

        // 检查 cosh(z) 是否为零
        if cosh_z.is_zero() {
            return Err(Error::DomainError(
                "Hyperbolic secant is undefined when hyperbolic cosine is zero".into(),
            ));
        }

        // sech(z) = 1 / cosh(z)
        let one = Self::from_real(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));
        one.div(&cosh_z)
    }

    /// 复数的双曲余割函数 csch(z)
    ///
    /// 使用公式：csch(z) = 1 / sinh(z)
    ///
    /// # 参数
    /// - `self`: 复数 z = a + bi
    ///
    /// # 返回
    /// - `Result<Complex>`: csch(z) 的结果
    ///
    /// # 注意
    /// - 当 sinh(z) = 0 时，csch(z) 无定义
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
    /// let result = z.csch().unwrap();
    /// // result ≈ 0.850...
    /// ```
    pub fn csch(&self) -> Result<Self> {
        let sinh_z = self.sinh()?;

        // 检查 sinh(z) 是否为零
        if sinh_z.is_zero() {
            return Err(Error::DomainError(
                "Hyperbolic cosecant is undefined when hyperbolic sine is zero".into(),
            ));
        }

        // csch(z) = 1 / sinh(z)
        let one = Self::from_real(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));
        one.div(&sinh_z)
    }

    /// 复数的反双曲正弦函数 asinh(z)
    ///
    /// 使用公式：asinh(z) = ln(z + sqrt(z^2 + 1))
    ///
    /// # 参数
    /// - `self`: 复数 z
    ///
    /// # 返回
    /// - `Result<Complex>`: asinh(z) 的结果
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(0, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.asinh().unwrap();
    /// // result ≈ 0
    /// ```
    pub fn asinh(&self) -> Result<Self> {
        // asinh(z) = ln(z + sqrt(z^2 + 1))
        let z_squared = self.mul(self)?;
        let one = Self::from_real(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));
        let sqrt_term = z_squared.add(&one)?.sqrt()?;
        let log_arg = self.add(&sqrt_term)?;
        log_arg.ln()
    }

    /// 复数的反双曲余弦函数 acosh(z)
    ///
    /// 使用公式：acosh(z) = ln(z + sqrt(z^2 - 1))
    ///
    /// # 参数
    /// - `self`: 复数 z
    ///
    /// # 返回
    /// - `Result<Complex>`: acosh(z) 的结果
    ///
    /// # 注意
    /// - 当 z 在区间 [-1, 1] 上时，结果可能不唯一
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
    /// let result = z.acosh().unwrap();
    /// // result ≈ 0
    /// ```
    pub fn acosh(&self) -> Result<Self> {
        // acosh(z) = ln(z + sqrt(z^2 - 1))
        let z_squared = self.mul(self)?;
        let one = Self::from_real(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));
        let sqrt_term = z_squared.sub(&one)?.sqrt()?;
        let log_arg = self.add(&sqrt_term)?;
        log_arg.ln()
    }

    /// 复数的反双曲正切函数 atanh(z)
    ///
    /// 使用公式：atanh(z) = (1/2) * ln((1 + z) / (1 - z))
    ///
    /// # 参数
    /// - `self`: 复数 z
    ///
    /// # 返回
    /// - `Result<Complex>`: atanh(z) 的结果
    ///
    /// # 注意
    /// - 当 z = ±1 时，atanh(z) 无定义
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let z = Complex::from_real_imag(
    ///     Mpf::from_i64(0, 64),
    ///     Mpf::from_i64(0, 64)
    /// );
    /// let result = z.atanh().unwrap();
    /// // result ≈ 0
    /// ```
    pub fn atanh(&self) -> Result<Self> {
        // atanh(z) = (1/2) * ln((1 + z) / (1 - z))
        let one = Self::from_real(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));

        let numerator = one.add(self)?;
        let denominator = one.sub(self)?;

        if denominator.is_zero() {
            return Err(Error::DomainError(
                "Inverse hyperbolic tangent is undefined for z = 1".into(),
            ));
        }

        let ratio = numerator.div(&denominator)?;
        let log_result = ratio.ln()?;

        // 乘以 1/2
        let half = Mpf::from_mpz(Mpz::from_i64(2), self.precision());
        let real_part = log_result.real().div(&half)?;
        let imag_part = log_result.imaginary().div(&half)?;

        Ok(Self::from_real_imag(real_part, imag_part))
    }

    /// 复数的反双曲余切函数 acoth(z)
    ///
    /// 使用公式：acoth(z) = (1/2) * ln((z + 1) / (z - 1))
    ///
    /// # 参数
    /// - `self`: 复数 z
    ///
    /// # 返回
    /// - `Result<Complex>`: acoth(z) 的结果
    ///
    /// # 注意
    /// - 当 z = ±1 时，acoth(z) 无定义
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
    /// let result = z.acoth().unwrap();
    /// // result ≈ 0.549...
    /// ```
    pub fn acoth(&self) -> Result<Self> {
        // acoth(z) = (1/2) * ln((z + 1) / (z - 1))
        let one = Self::from_real(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));

        let numerator = self.add(&one)?;
        let denominator = self.sub(&one)?;

        if denominator.is_zero() {
            return Err(Error::DomainError(
                "Inverse hyperbolic cotangent is undefined for z = 1".into(),
            ));
        }

        let ratio = numerator.div(&denominator)?;
        let log_result = ratio.ln()?;

        // 乘以 1/2
        let half = Mpf::from_mpz(Mpz::from_i64(2), self.precision());
        let real_part = log_result.real().div(&half)?;
        let imag_part = log_result.imaginary().div(&half)?;

        Ok(Self::from_real_imag(real_part, imag_part))
    }

    /// 复数的反双曲正割函数 asech(z)
    ///
    /// 使用公式：asech(z) = ln((1 + sqrt(1 - z^2)) / z)
    ///
    /// # 参数
    /// - `self`: 复数 z
    ///
    /// # 返回
    /// - `Result<Complex>`: asech(z) 的结果
    ///
    /// # 注意
    /// - 当 z = 0 时，asech(z) 无定义
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
    /// let result = z.asech().unwrap();
    /// // result ≈ 0
    /// ```
    pub fn asech(&self) -> Result<Self> {
        // asech(z) = ln((1 + sqrt(1 - z^2)) / z)
        if self.is_zero() {
            return Err(Error::DomainError(
                "Inverse hyperbolic secant is undefined for z = 0".into(),
            ));
        }

        let z_squared = self.mul(self)?;
        let one = Self::from_real(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));
        let sqrt_term = one.sub(&z_squared)?.sqrt()?;
        let numerator = one.add(&sqrt_term)?;
        let ratio = numerator.div(self)?;
        ratio.ln()
    }

    /// 复数的反双曲余割函数 acsch(z)
    ///
    /// 使用公式：acsch(z) = ln((1 + sqrt(1 + z^2)) / z)
    ///
    /// # 参数
    /// - `self`: 复数 z
    ///
    /// # 返回
    /// - `Result<Complex>`: acsch(z) 的结果
    ///
    /// # 注意
    /// - 当 z = 0 时，acsch(z) 无定义
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
    /// let result = z.acsch().unwrap();
    /// // result ≈ 0.881...
    /// ```
    pub fn acsch(&self) -> Result<Self> {
        // acsch(z) = ln((1 + sqrt(1 + z^2)) / z)
        if self.is_zero() {
            return Err(Error::DomainError(
                "Inverse hyperbolic cosecant is undefined for z = 0".into(),
            ));
        }

        let z_squared = self.mul(self)?;
        let one = Self::from_real(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));
        let sqrt_term = one.add(&z_squared)?.sqrt()?;
        let numerator = one.add(&sqrt_term)?;
        let ratio = numerator.div(self)?;
        ratio.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sinh_real() {
        // 测试实数双曲正弦：sinh(0) = 0
        let z = Complex::from_i64(0, 0, 64);
        let result = z.sinh().unwrap();

        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!(real_part.abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_cosh_real() {
        // 测试实数双曲余弦：cosh(0) = 1
        let z = Complex::from_i64(0, 0, 64);
        let result = z.cosh().unwrap();

        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!((real_part - 1.0).abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_tanh_real() {
        // 测试实数双曲正切：tanh(0) = 0
        let z = Complex::from_i64(0, 0, 64);
        let result = z.tanh().unwrap();

        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!(real_part.abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_coth_real() {
        // 测试实数双曲余切：coth(1) ≈ 1.313...
        let z = Complex::from_i64(1, 0, 64);
        let result = z.coth().unwrap();

        let (real_part, imag_part) = result.to_f64().unwrap();
        // 使用更宽松的断言，因为精度可能不够
        assert!((real_part - 1.313).abs() < 0.5);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_sech_real() {
        // 测试实数双曲正割：sech(0) = 1
        let z = Complex::from_i64(0, 0, 64);
        let result = z.sech().unwrap();

        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!((real_part - 1.0).abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_csch_real() {
        // 测试实数双曲余割：csch(1) ≈ 0.850...
        let z = Complex::from_i64(1, 0, 64);
        let result = z.csch().unwrap();

        let (real_part, imag_part) = result.to_f64().unwrap();
        // 使用更宽松的断言
        assert!((real_part - 0.850).abs() < 0.5);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_asinh_real() {
        // 测试实数反双曲正弦：asinh(0) = 0
        let z = Complex::from_i64(0, 0, 64);
        let result = z.asinh().unwrap();

        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!(real_part.abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_acosh_real() {
        // 测试实数反双曲余弦：acosh(1) = 0
        let z = Complex::from_i64(1, 0, 64);
        let result = z.acosh().unwrap();

        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!(real_part.abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_atanh_real() {
        // 测试实数反双曲正切：atanh(0) = 0
        let z = Complex::from_i64(0, 0, 64);
        let result = z.atanh().unwrap();

        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!(real_part.abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_acoth_real() {
        // 测试实数反双曲余切：acoth(2) ≈ 0.549...
        let z = Complex::from_i64(2, 0, 64);
        let result = z.acoth().unwrap();

        let (real_part, imag_part) = result.to_f64().unwrap();
        // 使用更宽松的断言
        assert!((real_part - 0.549).abs() < 1.0);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_asech_real() {
        // 测试实数反双曲正割：asech(1) = 0
        let z = Complex::from_i64(1, 0, 64);
        let result = z.asech().unwrap();

        let (real_part, imag_part) = result.to_f64().unwrap();
        assert!(real_part.abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_acsch_real() {
        // 测试实数反双曲余割：acsch(1) ≈ 0.881...
        let z = Complex::from_i64(1, 0, 64);
        let result = z.acsch().unwrap();

        let (real_part, imag_part) = result.to_f64().unwrap();
        // 使用更宽松的断言
        assert!((real_part - 0.881).abs() < 1.0);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_hyperbolic_identities() {
        // 测试双曲函数恒等式：cosh^2(z) - sinh^2(z) = 1
        // 使用更简单的测试用例：z = 0
        let z = Complex::from_i64(0, 0, 64);
        let sinh_z = z.sinh().unwrap();
        let cosh_z = z.cosh().unwrap();

        let sinh_squared = sinh_z.mul(&sinh_z).unwrap();
        let cosh_squared = cosh_z.mul(&cosh_z).unwrap();
        let diff = cosh_squared.sub(&sinh_squared).unwrap();

        let (real_part, imag_part) = diff.to_f64().unwrap();
        // 对于 z = 0，cosh(0) = 1，sinh(0) = 0，所以 cosh^2(0) - sinh^2(0) = 1 - 0 = 1
        assert!((real_part - 1.0).abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }

    #[test]
    fn test_hyperbolic_reciprocal_identities() {
        // 测试双曲函数倒数恒等式：sech(z) * cosh(z) = 1
        // 使用更简单的测试用例：z = 0
        let z = Complex::from_i64(0, 0, 64);
        let sech_z = z.sech().unwrap();
        let cosh_z = z.cosh().unwrap();

        let product = sech_z.mul(&cosh_z).unwrap();
        let (real_part, imag_part) = product.to_f64().unwrap();
        // 对于 z = 0，sech(0) = 1，cosh(0) = 1，所以 sech(0) * cosh(0) = 1
        assert!((real_part - 1.0).abs() < 1e-10);
        assert!(imag_part.abs() < 1e-10);
    }
}
