//! 复数多项式运算模块
//!
//! 包含复系数多项式的各种运算，如加法、乘法、求值、求根等

use crate::complex::core::Complex;
use crate::error::{Error, Result};
use crate::mpf::Mpf;

/// 复数多项式结构
///
/// 表示形式：P(z) = a₀ + a₁z + a₂z² + ... + aₙzⁿ
/// 其中 aᵢ 是复数系数
#[derive(Debug, Clone)]
pub struct ComplexPolynomial {
    /// 多项式系数，从低次项到高次项
    coefficients: Vec<Complex>,
}

impl ComplexPolynomial {
    /// 创建零多项式
    pub fn zero() -> Self {
        Self {
            coefficients: vec![Complex::new()],
        }
    }

    /// 创建常数多项式
    pub fn constant(c: Complex) -> Self {
        Self {
            coefficients: vec![c],
        }
    }

    /// 从系数向量创建多项式
    ///
    /// # 参数
    /// - `coeffs`: 系数向量，索引i对应z^i项的系数
    ///
    /// # 返回
    /// - `Result<Self>`: 多项式或错误
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::polynomial::ComplexPolynomial;
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let coeffs = vec![
    ///     Complex::from_real(Mpf::from_i64(1, 64)),  // 常数项
    ///     Complex::from_real(Mpf::from_i64(2, 64)),  // z项系数
    ///     Complex::from_real(Mpf::from_i64(1, 64)),  // z²项系数
    /// ];
    /// let poly = ComplexPolynomial::from_coefficients(coeffs).unwrap();
    /// // 表示多项式 1 + 2z + z²
    /// ```
    pub fn from_coefficients(coeffs: Vec<Complex>) -> Result<Self> {
        if coeffs.is_empty() {
            return Err(Error::InvalidInput("Coefficients cannot be empty".into()));
        }

        // 去除尾部的零系数（除非是零多项式）
        let mut trimmed = coeffs;
        while trimmed.len() > 1 && trimmed.last().unwrap().is_zero() {
            trimmed.pop();
        }

        Ok(Self {
            coefficients: trimmed,
        })
    }

    /// 获取多项式的次数（最高次项的指数）
    pub fn degree(&self) -> usize {
        if self.coefficients.len() <= 1 {
            0
        } else {
            self.coefficients.len() - 1
        }
    }

    /// 获取指定次项的系数
    pub fn coefficient(&self, degree: usize) -> Option<&Complex> {
        self.coefficients.get(degree)
    }

    /// 获取所有系数
    pub fn coefficients(&self) -> &[Complex] {
        &self.coefficients
    }

    /// 计算多项式在指定点的值
    ///
    /// 使用霍纳法则（Horner's method）进行高效计算
    ///
    /// # 参数
    /// - `z`: 要计算的点
    ///
    /// # 返回
    /// - `Result<Complex>`: 多项式在该点的值
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::polynomial::ComplexPolynomial;
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    ///
    /// let poly = ComplexPolynomial::from_coefficients(vec![
    ///     Complex::from_real(Mpf::from_i64(1, 64)),  // 常数项
    ///     Complex::from_real(Mpf::from_i64(2, 64)),  // z项系数
    ///     Complex::from_real(Mpf::from_i64(1, 64)),  // z²项系数
    /// ]).unwrap();
    ///
    /// let z = Complex::from_real(Mpf::from_i64(3, 64));
    /// let value = poly.evaluate(&z).unwrap();
    /// // value = 1 + 2*3 + 3² = 1 + 6 + 9 = 16
    /// ```
    pub fn evaluate(&self, z: &Complex) -> Result<Complex> {
        if self.coefficients.is_empty() {
            return Ok(Complex::new());
        }

        let mut result = self.coefficients.last().unwrap().clone();

        // 使用霍纳法则：((aₙz + aₙ₋₁)z + aₙ₋₂)z + ... + a₀
        for i in (0..self.coefficients.len() - 1).rev() {
            result = result.mul(z)?.add(&self.coefficients[i])?;
        }

        Ok(result)
    }

    /// 多项式加法
    ///
    /// # 参数
    /// - `other`: 另一个多项式
    ///
    /// # 返回
    /// - `Result<Self>`: 和多项式
    pub fn add(&self, other: &ComplexPolynomial) -> Result<Self> {
        let max_degree = std::cmp::max(self.degree(), other.degree());
        let mut result_coeffs = Vec::with_capacity(max_degree + 1);

        for i in 0..=max_degree {
            let self_coeff = self.coefficient(i).unwrap_or(&Complex::new()).clone();
            let other_coeff = other.coefficient(i).unwrap_or(&Complex::new()).clone();
            let sum = self_coeff.add(&other_coeff)?;
            result_coeffs.push(sum);
        }

        ComplexPolynomial::from_coefficients(result_coeffs)
    }

    /// 多项式减法
    ///
    /// # 参数
    /// - `other`: 另一个多项式
    ///
    /// # 返回
    /// - `Result<Self>`: 差多项式
    pub fn sub(&self, other: &ComplexPolynomial) -> Result<Self> {
        let max_degree = std::cmp::max(self.degree(), other.degree());
        let mut result_coeffs = Vec::with_capacity(max_degree + 1);

        for i in 0..=max_degree {
            let self_coeff = self.coefficient(i).unwrap_or(&Complex::new()).clone();
            let other_coeff = other.coefficient(i).unwrap_or(&Complex::new()).clone();
            let diff = self_coeff.sub(&other_coeff)?;
            result_coeffs.push(diff);
        }

        ComplexPolynomial::from_coefficients(result_coeffs)
    }

    /// 多项式乘法
    ///
    /// 使用卷积算法计算多项式乘积
    ///
    /// # 参数
    /// - `other`: 另一个多项式
    ///
    /// # 返回
    /// - `Result<Self>`: 积多项式
    pub fn mul(&self, other: &ComplexPolynomial) -> Result<Self> {
        let result_degree = self.degree() + other.degree();
        let mut result_coeffs = vec![Complex::new(); result_degree + 1];

        // 卷积计算：cᵢ = Σⱼ₌₀ⁱ aⱼ × bᵢ₋ⱼ
        for i in 0..=self.degree() {
            for j in 0..=other.degree() {
                let product = self
                    .coefficient(i)
                    .unwrap()
                    .mul(other.coefficient(j).unwrap())?;
                let target_idx = i + j;
                result_coeffs[target_idx] = result_coeffs[target_idx].add(&product)?;
            }
        }

        ComplexPolynomial::from_coefficients(result_coeffs)
    }

    /// 多项式除法（带余数）
    ///
    /// 使用长除法算法
    ///
    /// # 参数
    /// - `divisor`: 除数多项式
    ///
    /// # 返回
    /// - `Result<(Self, Self)>`: (商多项式, 余数多项式)
    ///
    /// # 注意
    /// - 除数不能为零多项式
    /// - 如果被除数的次数小于除数的次数，商为零多项式，余数为被除数
    pub fn div_with_remainder(&self, divisor: &ComplexPolynomial) -> Result<(Self, Self)> {
        if divisor.is_zero() {
            return Err(Error::DivisionByZero);
        }

        if self.degree() < divisor.degree() {
            return Ok((ComplexPolynomial::zero(), self.clone()));
        }

        let mut dividend = self.clone();
        let mut quotient_coeffs = vec![Complex::new(); self.degree() - divisor.degree() + 1];

        while dividend.degree() >= divisor.degree() {
            let current_degree = dividend.degree();
            let divisor_degree = divisor.degree();
            let degree_diff = current_degree - divisor_degree;

            // 计算当前项的商
            let dividend_lead = dividend.coefficient(current_degree).unwrap();
            let divisor_lead = divisor.coefficient(divisor_degree).unwrap();
            let current_quotient = dividend_lead.div(divisor_lead)?;

            // 更新商多项式
            quotient_coeffs[degree_diff] = current_quotient.clone();

            // 从被除数中减去当前项
            let mut term_coeffs = vec![Complex::new(); current_degree + 1];
            term_coeffs[degree_diff] = current_quotient;
            let term_poly = ComplexPolynomial::from_coefficients(term_coeffs)?;
            let product = term_poly.mul(divisor)?;
            dividend = dividend.sub(&product)?;
        }

        let quotient = ComplexPolynomial::from_coefficients(quotient_coeffs)?;
        Ok((quotient, dividend))
    }

    /// 多项式除法
    ///
    /// # 参数
    /// - `divisor`: 除数多项式
    ///
    /// # 返回
    /// - `Result<Self>`: 商多项式
    pub fn div(&self, divisor: &ComplexPolynomial) -> Result<Self> {
        let (quotient, _) = self.div_with_remainder(divisor)?;
        Ok(quotient)
    }

    /// 多项式求余
    ///
    /// # 参数
    /// - `divisor`: 除数多项式
    ///
    /// # 返回
    /// - `Result<Self>`: 余数多项式
    pub fn rem(&self, divisor: &ComplexPolynomial) -> Result<Self> {
        let (_, remainder) = self.div_with_remainder(divisor)?;
        Ok(remainder)
    }

    /// 计算多项式的导数
    ///
    /// # 返回
    /// - `Result<Self>`: 导数多项式
    pub fn derivative(&self) -> Result<Self> {
        if self.degree() == 0 {
            return Ok(ComplexPolynomial::zero());
        }

        let mut deriv_coeffs = Vec::with_capacity(self.degree());

        for i in 1..=self.degree() {
            let coeff = self.coefficient(i).unwrap();
            let power = Mpf::from_i64(i as i64, coeff.precision());
            let deriv_coeff = coeff.mul(&Complex::from_real(power))?;
            deriv_coeffs.push(deriv_coeff);
        }

        ComplexPolynomial::from_coefficients(deriv_coeffs)
    }

    /// 计算多项式的不定积分
    ///
    /// # 参数
    /// - `constant`: 积分常数
    ///
    /// # 返回
    /// - `Result<Self>`: 积分多项式
    pub fn integral(&self, constant: &Complex) -> Result<Self> {
        let mut integral_coeffs = Vec::with_capacity(self.degree() + 2);
        integral_coeffs.push(constant.clone());

        for i in 0..=self.degree() {
            let coeff = self.coefficient(i).unwrap();
            let power = Mpf::from_i64((i + 1) as i64, coeff.precision());
            let integral_coeff = coeff.div(&Complex::from_real(power))?;
            integral_coeffs.push(integral_coeff);
        }

        ComplexPolynomial::from_coefficients(integral_coeffs)
    }

    /// 检查是否为零多项式
    pub fn is_zero(&self) -> bool {
        self.coefficients.len() == 1 && self.coefficients[0].is_zero()
    }
}

impl std::fmt::Display for ComplexPolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = if self.is_zero() {
            "0".to_string()
        } else {
            let mut terms = Vec::new();

            for (i, coeff) in self.coefficients.iter().enumerate() {
                if coeff.is_zero() {
                    continue;
                }

                let term = if i == 0 {
                    coeff.to_string_formatted("cartesian")
                } else if i == 1 {
                    if coeff.to_f64().unwrap_or((0.0, 0.0)).0 == 1.0
                        && coeff.to_f64().unwrap_or((0.0, 0.0)).1 == 0.0
                    {
                        "z".to_string()
                    } else if coeff.to_f64().unwrap_or((0.0, 0.0)).0 == -1.0
                        && coeff.to_f64().unwrap_or((0.0, 0.0)).1 == 0.0
                    {
                        "-z".to_string()
                    } else {
                        format!("{}z", coeff.to_string_formatted("cartesian"))
                    }
                } else {
                    if coeff.to_f64().unwrap_or((0.0, 0.0)).0 == 1.0
                        && coeff.to_f64().unwrap_or((0.0, 0.0)).1 == 0.0
                    {
                        format!("z^{}", i)
                    } else if coeff.to_f64().unwrap_or((0.0, 0.0)).0 == -1.0
                        && coeff.to_f64().unwrap_or((0.0, 0.0)).1 == 0.0
                    {
                        format!("-z^{}", i)
                    } else {
                        format!("{}z^{}", coeff.to_string_formatted("cartesian"), i)
                    }
                };

                terms.push(term);
            }

            if terms.is_empty() {
                "0".to_string()
            } else {
                terms.join(" + ")
            }
        };
        write!(f, "{}", s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_creation() {
        let poly = ComplexPolynomial::constant(Complex::from_i64(5, 0, 64));
        assert_eq!(poly.degree(), 0);
        assert_eq!(poly.coefficient(0).unwrap().to_f64().unwrap().0, 5.0);
    }

    #[test]
    fn test_polynomial_from_coefficients() {
        let coeffs = vec![
            Complex::from_i64(1, 0, 64),
            Complex::from_i64(2, 0, 64),
            Complex::from_i64(1, 0, 64),
        ];
        let poly = ComplexPolynomial::from_coefficients(coeffs).unwrap();
        assert_eq!(poly.degree(), 2);
        assert_eq!(poly.coefficient(0).unwrap().to_f64().unwrap().0, 1.0);
        assert_eq!(poly.coefficient(1).unwrap().to_f64().unwrap().0, 2.0);
        assert_eq!(poly.coefficient(2).unwrap().to_f64().unwrap().0, 1.0);
    }

    #[test]
    fn test_polynomial_evaluation() {
        let poly = ComplexPolynomial::from_coefficients(vec![
            Complex::from_i64(1, 0, 64),
            Complex::from_i64(2, 0, 64),
            Complex::from_i64(1, 0, 64),
        ])
        .unwrap();

        let z = Complex::from_i64(3, 0, 64);
        let value = poly.evaluate(&z).unwrap();
        // 1 + 2*3 + 3² = 1 + 6 + 9 = 16
        assert_eq!(value.to_f64().unwrap().0, 16.0);
    }

    #[test]
    fn test_polynomial_addition() {
        let poly1 = ComplexPolynomial::from_coefficients(vec![
            Complex::from_i64(1, 0, 64),
            Complex::from_i64(2, 0, 64),
        ])
        .unwrap();

        let poly2 = ComplexPolynomial::from_coefficients(vec![
            Complex::from_i64(3, 0, 64),
            Complex::from_i64(4, 0, 64),
        ])
        .unwrap();

        let sum = poly1.add(&poly2).unwrap();
        assert_eq!(sum.coefficient(0).unwrap().to_f64().unwrap().0, 4.0);
        assert_eq!(sum.coefficient(1).unwrap().to_f64().unwrap().0, 6.0);
    }

    #[test]
    fn test_polynomial_multiplication() {
        let poly1 = ComplexPolynomial::from_coefficients(vec![
            Complex::from_i64(1, 0, 64),
            Complex::from_i64(1, 0, 64),
        ])
        .unwrap();

        let poly2 = ComplexPolynomial::from_coefficients(vec![
            Complex::from_i64(1, 0, 64),
            Complex::from_i64(1, 0, 64),
        ])
        .unwrap();

        let product = poly1.mul(&poly2).unwrap();
        // (1 + z)(1 + z) = 1 + 2z + z²
        assert_eq!(product.degree(), 2);
        assert_eq!(product.coefficient(0).unwrap().to_f64().unwrap().0, 1.0);
        assert_eq!(product.coefficient(1).unwrap().to_f64().unwrap().0, 2.0);
        assert_eq!(product.coefficient(2).unwrap().to_f64().unwrap().0, 1.0);
    }

    #[test]
    fn test_polynomial_derivative() {
        let poly = ComplexPolynomial::from_coefficients(vec![
            Complex::from_i64(1, 0, 64),
            Complex::from_i64(2, 0, 64),
            Complex::from_i64(3, 0, 64),
        ])
        .unwrap();

        let derivative = poly.derivative().unwrap();
        // 导数：2 + 6z
        assert_eq!(derivative.degree(), 1);
        assert_eq!(derivative.coefficient(0).unwrap().to_f64().unwrap().0, 2.0);
        assert_eq!(derivative.coefficient(1).unwrap().to_f64().unwrap().0, 6.0);
    }

    #[test]
    fn test_polynomial_integral() {
        let poly = ComplexPolynomial::from_coefficients(vec![
            Complex::from_i64(2, 0, 64),
            Complex::from_i64(6, 0, 64),
        ])
        .unwrap();

        let constant = Complex::from_i64(1, 0, 64);
        let integral = poly.integral(&constant).unwrap();
        // 积分：1 + 2z + 3z²
        assert_eq!(integral.degree(), 2);
        assert_eq!(integral.coefficient(0).unwrap().to_f64().unwrap().0, 1.0);
        assert_eq!(integral.coefficient(1).unwrap().to_f64().unwrap().0, 2.0);
        assert_eq!(integral.coefficient(2).unwrap().to_f64().unwrap().0, 3.0);
    }

    #[test]
    fn test_polynomial_division() {
        let dividend = ComplexPolynomial::from_coefficients(vec![
            Complex::from_i64(1, 0, 64),
            Complex::from_i64(0, 0, 64),
            Complex::from_i64(-1, 0, 64),
        ])
        .unwrap(); // 1 - z²

        let divisor = ComplexPolynomial::from_coefficients(vec![
            Complex::from_i64(1, 0, 64),
            Complex::from_i64(1, 0, 64),
        ])
        .unwrap(); // 1 + z

        let (quotient, remainder) = dividend.div_with_remainder(&divisor).unwrap();
        // (1 - z²) / (1 + z) = 1 - z, 余数 0
        assert_eq!(quotient.degree(), 1);
        assert_eq!(quotient.coefficient(0).unwrap().to_f64().unwrap().0, 1.0);
        assert_eq!(quotient.coefficient(1).unwrap().to_f64().unwrap().0, -1.0);
        assert!(remainder.is_zero());
    }

    #[test]
    fn test_polynomial_to_string() {
        let poly = ComplexPolynomial::from_coefficients(vec![
            Complex::from_i64(1, 0, 64),
            Complex::from_i64(2, 0, 64),
            Complex::from_i64(1, 0, 64),
        ])
        .unwrap();

        let s = poly.to_string();
        assert!(s.contains("1"));
        assert!(s.contains("2z"));
        assert!(s.contains("z^2"));
    }
}
