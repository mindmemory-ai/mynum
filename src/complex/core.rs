//! 复数核心结构体和基本操作

use super::ComplexConfig;
use crate::error::{Error, Result};
use crate::mpf::Mpf;
use crate::mpz::Mpz;
use core::fmt;

/// 高精度复数结构
///
/// 复数由实部和虚部组成，都是高精度浮点数 Mpf
#[derive(Clone, Debug)]
pub struct Complex {
    /// 实部
    real: Mpf,
    /// 虚部
    imaginary: Mpf,
    /// 精度
    precision: usize,
    /// 配置
    config: ComplexConfig,
}

impl Complex {
    /// 创建新的复数 (0 + 0i)
    pub fn new() -> Self {
        Self {
            real: Mpf::new(),
            imaginary: Mpf::new(),
            precision: 64,
            config: ComplexConfig,
        }
    }

    /// 从实部和虚部创建复数
    pub fn from_real_imag(real: Mpf, imaginary: Mpf) -> Self {
        let precision = real.precision().max(imaginary.precision());
        Self {
            real,
            imaginary,
            precision,
            config: ComplexConfig,
        }
    }

    /// 从实数创建复数 (实数 + 0i)
    pub fn from_real(real: Mpf) -> Self {
        let precision = real.precision();
        Self {
            real,
            imaginary: Mpf::new(),
            precision,
            config: ComplexConfig,
        }
    }

    /// 从虚数创建复数 (0 + 虚数i)
    pub fn from_imag(imaginary: Mpf) -> Self {
        let precision = imaginary.precision();
        Self {
            real: Mpf::new(),
            imaginary,
            precision,
            config: ComplexConfig,
        }
    }

    /// 从两个整数创建复数
    pub fn from_i64(real: i64, imag: i64, precision: usize) -> Self {
        let mut real_mpf = Mpf::from_mpz(Mpz::from_i64(real.abs()), precision);
        if real < 0 {
            real_mpf = real_mpf.neg();
        }

        let mut imag_mpf = Mpf::from_mpz(Mpz::from_i64(imag.abs()), precision);
        if imag < 0 {
            imag_mpf = imag_mpf.neg();
        }

        Self {
            real: real_mpf,
            imaginary: imag_mpf,
            precision,
            config: ComplexConfig,
        }
    }

    /// 从字符串解析复数 "a+bi", "a-bi", "a", "bi"
    pub fn from_str(s: &str, precision: usize) -> Result<Self> {
        let s = s.trim().replace(" ", "");

        // 处理纯虚数 "bi" 或 "i"
        if s == "i" {
            return Ok(Self::from_imag(Mpf::from_mpz(Mpz::from_i64(1), precision)));
        }
        if s == "-i" {
            return Ok(Self::from_imag(Mpf::from_mpz(Mpz::from_i64(-1), precision)));
        }
        if s.ends_with('i') && !s.contains('+') && !s.contains('-') {
            let imag_str = &s[..s.len() - 1];
            if imag_str.is_empty() {
                return Ok(Self::from_imag(Mpf::from_mpz(Mpz::from_i64(1), precision)));
            }
            let imag = Mpf::from_str(imag_str, 10)?;
            return Ok(Self::from_imag(imag));
        }

        // 处理纯实数
        if !s.contains('i') {
            let real = Mpf::from_str(&s, 10)?;
            return Ok(Self::from_real(real));
        }

        // 处理复数 "a+bi" 或 "a-bi"
        let (real_str, imag_str) = if let Some(pos) = s.rfind('+') {
            (&s[..pos], &s[pos + 1..])
        } else if let Some(pos) = s.rfind('-') {
            if pos == 0 {
                // 负号在开头，可能是 "-a+bi" 或 "-bi"
                if let Some(pos2) = s[1..].find('+') {
                    (&s[..pos2 + 1], &s[pos2 + 2..])
                } else if let Some(pos2) = s[1..].find('-') {
                    (&s[..pos2 + 1], &s[pos2 + 1..])
                } else {
                    // 纯虚数 "-bi"
                    ("0", s.as_str())
                }
            } else {
                (&s[..pos], &s[pos..])
            }
        } else {
            return Err(Error::ParseError("Invalid complex number format".into()));
        };

        let real = if real_str.is_empty() || real_str == "+" {
            Mpf::new()
        } else {
            Mpf::from_str(real_str, 10)?
        };

        let imag = if imag_str == "i" || imag_str == "+i" {
            Mpf::from_mpz(Mpz::from_i64(1), precision)
        } else if imag_str == "-i" {
            Mpf::from_mpz(Mpz::from_i64(-1), precision)
        } else if imag_str.ends_with('i') {
            let imag_num_str = imag_str.strip_suffix('i').unwrap();
            if imag_num_str.is_empty() || imag_num_str == "+" {
                Mpf::from_mpz(Mpz::from_i64(1), precision)
            } else if imag_num_str == "-" {
                Mpf::from_mpz(Mpz::from_i64(-1), precision)
            } else {
                Mpf::from_str(imag_num_str, 10)?
            }
        } else {
            return Err(Error::ParseError("Invalid complex number format".into()));
        };

        Ok(Self::from_real_imag(real, imag))
    }

    /// 从配置创建复数
    pub fn with_config(mut self, config: ComplexConfig) -> Self {
        self.config = config;
        self
    }

    /// 获取实部
    pub fn real(&self) -> &Mpf {
        &self.real
    }

    /// 获取虚部
    pub fn imaginary(&self) -> &Mpf {
        &self.imaginary
    }

    /// 获取精度
    pub fn precision(&self) -> usize {
        self.precision
    }

    /// 获取配置
    pub fn config(&self) -> &ComplexConfig {
        &self.config
    }

    /// 设置精度
    pub fn with_precision(mut self, precision: usize) -> Self {
        self.precision = precision;
        self
    }

    /// 判断是否为零
    pub fn is_zero(&self) -> bool {
        self.real.is_zero() && self.imaginary.is_zero()
    }

    /// 判断是否为实数
    pub fn is_real(&self) -> bool {
        self.imaginary.is_zero()
    }

    /// 判断是否为纯虚数
    pub fn is_imaginary(&self) -> bool {
        self.real().is_zero()
    }

    /// 判断是否为纯实数（包括零）
    pub fn is_pure_real(&self) -> bool {
        self.imaginary.is_zero()
    }

    /// 判断是否为纯虚数（包括零）
    pub fn is_pure_imaginary(&self) -> bool {
        self.real.is_zero()
    }

    /// 判断是否为单位复数（模长为1）
    pub fn is_unit(&self) -> Result<bool> {
        let mag = self.magnitude()?;
        Ok(mag == Mpf::from_i64(1, self.precision))
    }

    /// 判断是否为实数1
    pub fn is_one(&self) -> bool {
        self.real() == &Mpf::from_i64(1, self.precision) && self.imaginary().is_zero()
    }

    /// 判断是否为实数-1
    pub fn is_neg_one(&self) -> bool {
        self.real() == &Mpf::from_i64(-1, self.precision) && self.imaginary().is_zero()
    }

    /// 获取复数的类型描述
    pub fn type_description(&self) -> &'static str {
        if self.is_zero() {
            "zero"
        } else if self.is_real() {
            "real"
        } else if self.is_imaginary() {
            "pure imaginary"
        } else {
            "complex"
        }
    }

    /// 转换为极坐标形式 (r, θ)
    pub fn to_polar(&self) -> Result<(Mpf, Mpf)> {
        let r = self.magnitude()?;
        let theta = self.argument()?;
        Ok((r, theta))
    }

    /// 从极坐标创建复数
    pub fn from_polar(r: Mpf, theta: Mpf, _precision: usize) -> Result<Self> {
        let real = r.mul(&theta.cos()?);
        let imaginary = r.mul(&theta.sin()?);
        Ok(Self::from_real_imag(real, imaginary))
    }

    /// 获取复数的实部引用（可变）
    pub fn real_mut(&mut self) -> &mut Mpf {
        &mut self.real
    }

    /// 获取复数的虚部引用（可变）
    pub fn imaginary_mut(&mut self) -> &mut Mpf {
        &mut self.imaginary
    }

    /// 设置实部
    pub fn set_real(&mut self, real: Mpf) {
        let precision = real.precision();
        self.real = real;
        self.precision = self.precision.max(precision);
    }

    /// 设置虚部
    pub fn set_imaginary(&mut self, imaginary: Mpf) {
        let precision = imaginary.precision();
        self.imaginary = imaginary;
        self.precision = self.precision.max(precision);
    }

    /// 设置实部和虚部
    pub fn set_real_imag(&mut self, real: Mpf, imaginary: Mpf) {
        let real_precision = real.precision();
        let imag_precision = imaginary.precision();
        self.real = real;
        self.imaginary = imaginary;
        self.precision = self.precision.max(real_precision).max(imag_precision);
    }

    /// 交换实部和虚部
    pub fn swap_real_imag(&mut self) {
        std::mem::swap(&mut self.real, &mut self.imaginary);
    }

    /// 获取复数的位长度（用于精度估计）
    pub fn bit_length(&self) -> usize {
        // 使用实部和虚部的精度作为位长度估计
        self.precision
    }

    /// 检查精度是否足够
    pub fn has_sufficient_precision(&self, required_precision: usize) -> bool {
        self.precision >= required_precision
    }

    /// 提升精度到指定值
    pub fn with_enhanced_precision(mut self, new_precision: usize) -> Self {
        if new_precision > self.precision {
            self.precision = new_precision;
            // 这里可以添加精度提升的逻辑
        }
        self
    }

    /// 转换为f64类型（如果可能）
    ///
    /// 返回一个元组 (real_part, imaginary_part)，如果转换失败则返回None
    pub fn to_f64(&self) -> Option<(f64, f64)> {
        let real_f64 = self.real.to_f64()?;
        let imag_f64 = self.imaginary.to_f64()?;
        Some((real_f64, imag_f64))
    }
}

impl fmt::Display for Complex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = if self.is_zero() {
            "0".to_string()
        } else if self.is_real() {
            self.real().to_string(10)
        } else if self.is_imaginary() {
            if self.imaginary() == &Mpf::from_i64(1, self.precision) {
                "i".to_string()
            } else if self.imaginary() == &Mpf::from_i64(-1, self.precision) {
                "-i".to_string()
            } else {
                format!("{}i", self.imaginary().to_string(10))
            }
        } else {
            let real_str = if self.real().is_zero() {
                "".to_string()
            } else {
                self.real().to_string(10)
            };

            let imag_str = if self.imaginary() == &Mpf::from_i64(1, self.precision) {
                "i".to_string()
            } else if self.imaginary() == &Mpf::from_i64(-1, self.precision) {
                "-i".to_string()
            } else {
                let sign = if self.imaginary().is_negative() {
                    ""
                } else {
                    "+"
                };
                format!("{}{}i", sign, self.imaginary().to_string(10))
            };

            if real_str.is_empty() {
                imag_str
            } else {
                format!("{}{}", real_str, imag_str)
            }
        };
        write!(f, "{}", s)
    }
}

impl Default for Complex {
    fn default() -> Self {
        Self::new()
    }
}

impl PartialEq for Complex {
    fn eq(&self, other: &Self) -> bool {
        self.real == other.real && self.imaginary == other.imaginary
    }
}

impl Eq for Complex {}

impl PartialOrd for Complex {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        // 复数没有自然的全序，但可以比较模长
        let mag1 = self.magnitude().ok()?;
        let mag2 = other.magnitude().ok()?;
        mag1.partial_cmp(&mag2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complex_creation() {
        let zero = Complex::new();
        assert!(zero.is_zero());
        assert!(zero.is_real());
        // 零既是实数又是纯虚数
        assert!(zero.is_imaginary());

        let real = Complex::from_real(Mpf::from_i64(42, 64));
        assert!(real.is_real());
        assert!(!real.is_imaginary());
        assert_eq!(real.real().to_i64(), Some(42));

        let imag = Complex::from_imag(Mpf::from_i64(7, 64));
        assert!(!imag.is_real());
        assert!(imag.is_imaginary());
        assert_eq!(imag.imaginary().to_i64(), Some(7));
    }

    #[test]
    fn test_complex_from_str() {
        let z1 = Complex::from_str("3+4i", 64).unwrap();
        assert_eq!(z1.real().to_i64(), Some(3));
        assert_eq!(z1.imaginary().to_i64(), Some(4));

        let z2 = Complex::from_str("5-2i", 64).unwrap();
        assert_eq!(z2.real().to_i64(), Some(5));
        assert_eq!(z2.imaginary().to_i64(), Some(-2));

        let z3 = Complex::from_str("7", 64).unwrap();
        assert!(z3.is_real());
        assert_eq!(z3.real().to_i64(), Some(7));

        let z4 = Complex::from_str("3i", 64).unwrap();
        assert!(z4.is_imaginary());
        assert_eq!(z4.imaginary().to_i64(), Some(3));
    }

    #[test]
    fn test_complex_properties() {
        let z = Complex::from_str("3+4i", 64).unwrap();
        assert!(!z.is_zero());
        assert!(!z.is_real());
        assert!(!z.is_imaginary());
        assert_eq!(z.type_description(), "complex");

        let zero = Complex::new();
        assert_eq!(zero.type_description(), "zero");

        let real = Complex::from_real(Mpf::from_i64(5, 64));
        assert_eq!(real.type_description(), "real");
    }

    #[test]
    fn test_complex_to_string() {
        let z1 = Complex::from_str("3+4i", 64).unwrap();
        assert_eq!(z1.to_string(), "3+4i");

        let z2 = Complex::from_str("5-2i", 64).unwrap();
        assert_eq!(z2.to_string(), "5-2i");

        let z3 = Complex::from_str("7", 64).unwrap();
        assert_eq!(z3.to_string(), "7");

        let z4 = Complex::from_str("3i", 64).unwrap();
        assert_eq!(z4.to_string(), "3i");

        let z5 = Complex::from_str("-2i", 64).unwrap();
        assert_eq!(z5.to_string(), "-2i");
    }
}
