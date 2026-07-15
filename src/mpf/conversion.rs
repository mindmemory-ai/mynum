//! Mpf 转换功能实现

use super::core::Mpf;
use crate::error::Result;
use crate::mpz::Mpz;

impl Mpf {
    /// 从字符串创建（支持科学计数法）
    pub fn from_str_radix(s: &str, base: u32) -> Result<Self> {
        Self::from_str(s, base)
    }

    /// 转换为字符串（指定进制）
    pub fn to_string_radix(&self, base: u32) -> String {
        self.to_string(base)
    }

    /// 转换为 f32
    pub fn to_f32(&self) -> Option<f32> {
        self.to_f64().map(|f| f as f32)
    }

    // to_f64 方法已经在 core.rs 中定义

    /// 转换为 Mpz（截断小数部分）
    pub fn to_mpz(&self) -> Mpz {
        if self.exponent() >= 0 {
            // 指数非负，直接返回尾数
            self.mantissa().clone()
        } else {
            // 指数为负，需要右移
            let shift = (-self.exponent()) as usize;
            self.mantissa().shr(shift)
        }
    }

    // from_mpz 方法已经在 core.rs 中定义

    /// 从整数创建
    pub fn from_i64(n: i64, precision: usize) -> Self {
        let (negative, abs_n) = if n < 0 {
            (true, n.wrapping_neg() as u64)
        } else {
            (false, n as u64)
        };

        let mut result = Mpf::from_mpz(Mpz::from_u64(abs_n), precision);
        result.set_negative(negative);
        result
    }

    /// 从无符号整数创建
    pub fn from_u64(n: u64, precision: usize) -> Self {
        Mpf::from_mpz(Mpz::from_u64(n), precision)
    }

    /// 转换为 i64（可能丢失精度）
    pub fn to_i64(&self) -> Option<i64> {
        if self.is_zero() {
            return Some(0);
        }

        if self.is_infinity() || self.is_nan() {
            return None;
        }

        let mpz = self.to_mpz();
        let result = mpz.to_i64()?;
        Some(if self.is_negative() { -result } else { result })
    }

    /// 转换为 u64（可能丢失精度）
    pub fn to_u64(&self) -> Option<u64> {
        if self.is_zero() {
            return Some(0);
        }

        if self.is_negative() || self.is_infinity() || self.is_nan() {
            return None;
        }

        self.to_mpz().to_u64()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_to_i64() {
        let mpf = Mpf::from_i64(123, 64);
        assert_eq!(mpf.to_i64(), Some(123));

        let mpf_neg = Mpf::from_i64(-456, 64);
        assert_eq!(mpf_neg.to_i64(), Some(-456));
    }

    #[test]
    fn test_from_to_u64() {
        let mpf = Mpf::from_u64(789, 64);
        assert_eq!(mpf.to_u64(), Some(789));
    }

    #[test]
    fn test_from_to_mpz() {
        let mpz = Mpz::from_i64(123);
        let mpf = Mpf::from_mpz(mpz.clone(), 64);
        let back_to_mpz = mpf.to_mpz();
        assert_eq!(back_to_mpz, mpz);
    }

    #[test]
    fn test_to_f64() {
        let mpf = Mpf::from_i64(123, 64);
        assert_eq!(mpf.to_f64(), Some(123.0));

        let mpf_neg = Mpf::from_i64(-456, 64);
        assert_eq!(mpf_neg.to_f64(), Some(-456.0));
    }
}
