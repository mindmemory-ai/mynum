//! Mpz核心实现
//!
//! 定义了大整数的内部表示和基础操作。

#[cfg(not(feature = "std"))]
use alloc::{string::String, vec::Vec};

use crate::error::{Error, Result};

/// 大整数类型
///
/// 使用小端序存储，limbs\[0\]是最低位
/// 符号单独存储，支持任意精度计算
#[derive(Debug, Clone)]
pub struct Mpz {
    /// 数据limbs，64位为单位，小端序
    limbs: Vec<u64>,
    /// 符号：false为正数或零，true为负数
    negative: bool,
}

impl Mpz {
    /// 创建零值
    pub fn new() -> Self {
        Self {
            limbs: vec![0],
            negative: false,
        }
    }

    /// 创建指定值的大整数
    pub fn from_limbs(limbs: Vec<u64>, negative: bool) -> Self {
        let mut result = Self { limbs, negative };
        result.normalize();
        result
    }

    /// 从字符串创建（支持十进制、十六进制、二进制）
    pub fn from_str(s: &str, base: u32) -> Result<Self> {
        if !(2..=36).contains(&base) {
            return Err(Error::InvalidInput("Base must be between 2 and 36".into()));
        }

        let s = s.trim();
        if s.is_empty() {
            return Err(Error::InvalidInput("Empty string".into()));
        }

        let (negative, s) = if let Some(rest) = s.strip_prefix('-') {
            (true, rest)
        } else if let Some(rest) = s.strip_prefix('+') {
            (false, rest)
        } else {
            (false, s)
        };

        if s.is_empty() {
            return Err(Error::InvalidInput("No digits after sign".into()));
        }

        // 处理不同进制的前缀
        let (actual_base, s) = if base == 16 && (s.starts_with("0x") || s.starts_with("0X")) {
            (16, &s[2..])
        } else if base == 2 && (s.starts_with("0b") || s.starts_with("0B")) {
            (2, &s[2..])
        } else if base == 8 && s.starts_with("0o") {
            (8, &s[2..])
        } else {
            (base, s)
        };

        if s.is_empty() {
            return Err(Error::InvalidInput("No digits after prefix".into()));
        }

        let mut result = Self::new();
        let base_mpz = Self::from_u64(actual_base as u64);

        for ch in s.chars() {
            let digit = ch.to_digit(actual_base).ok_or_else(|| {
                Error::InvalidInput(format!("Invalid digit '{}' for base {}", ch, actual_base))
            })?;

            result = result.mul(&base_mpz).add(&Self::from_u64(digit as u64));
        }

        result.negative = negative && !result.is_zero();
        Ok(result)
    }

    /// 从64位有符号整数创建
    pub fn from_i64(n: i64) -> Self {
        if n == 0 {
            return Self::new();
        }

        let (negative, abs_n) = if n < 0 {
            (true, n.wrapping_neg() as u64)
        } else {
            (false, n as u64)
        };

        Self {
            limbs: vec![abs_n],
            negative,
        }
    }

    /// 从64位无符号整数创建
    pub fn from_u64(n: u64) -> Self {
        if n == 0 {
            Self::new()
        } else {
            Self {
                limbs: vec![n],
                negative: false,
            }
        }
    }

    /// 从128位有符号整数创建
    pub fn from_i128(n: i128) -> Self {
        if n == 0 {
            return Self::new();
        }

        let (negative, abs_n) = if n < 0 {
            (true, n.wrapping_neg() as u128)
        } else {
            (false, n as u128)
        };

        let low = abs_n as u64;
        let high = (abs_n >> 64) as u64;

        let limbs = if high == 0 {
            vec![low]
        } else {
            vec![low, high]
        };

        Self { limbs, negative }
    }

    /// 从128位无符号整数创建
    pub fn from_u128(n: u128) -> Self {
        if n == 0 {
            return Self::new();
        }

        let low = n as u64;
        let high = (n >> 64) as u64;

        let limbs = if high == 0 {
            vec![low]
        } else {
            vec![low, high]
        };

        Self {
            limbs,
            negative: false,
        }
    }

    /// 转换为字符串
    pub fn to_string(&self, base: u32) -> String {
        if !(2..=36).contains(&base) {
            return "Invalid base".to_string();
        }

        if self.is_zero() {
            return "0".to_string();
        }

        let mut digits = Vec::new();
        let mut temp = self.abs();
        let base_mpz = Self::from_u64(base as u64);

        while !temp.is_zero() {
            let (quotient, remainder) = temp
                .div_rem(&base_mpz)
                .unwrap_or_else(|_| (Self::new(), Self::new()));
            let digit = remainder.to_u64().unwrap_or(0) as u8;

            let ch = if digit < 10 {
                (b'0' + digit) as char
            } else {
                (b'a' + digit - 10) as char
            };

            digits.push(ch);
            temp = quotient;
        }

        digits.reverse();
        let mut result: String = digits.into_iter().collect();

        if self.negative {
            result.insert(0, '-');
        }

        result
    }

    /// 转换为64位有符号整数
    pub fn to_i64(&self) -> Option<i64> {
        if self.limbs.len() > 1 {
            return None;
        }

        let abs_val = self.limbs[0];

        if self.negative {
            if abs_val > (i64::MAX as u64) + 1 {
                None
            } else if abs_val == (i64::MAX as u64) + 1 {
                Some(i64::MIN)
            } else {
                Some(-(abs_val as i64))
            }
        } else {
            if abs_val > i64::MAX as u64 {
                None
            } else {
                Some(abs_val as i64)
            }
        }
    }

    /// 转换为64位无符号整数
    pub fn to_u64(&self) -> Option<u64> {
        if self.negative || self.limbs.len() > 1 {
            return None;
        }

        Some(self.limbs[0])
    }

    /// 转换为128位有符号整数
    pub fn to_i128(&self) -> Option<i128> {
        if self.limbs.len() > 2 {
            return None;
        }

        let abs_val = if self.limbs.len() == 1 {
            self.limbs[0] as u128
        } else {
            (self.limbs[1] as u128) << 64 | (self.limbs[0] as u128)
        };

        if self.negative {
            if abs_val > (i128::MAX as u128) + 1 {
                None
            } else if abs_val == (i128::MAX as u128) + 1 {
                Some(i128::MIN)
            } else {
                Some(-(abs_val as i128))
            }
        } else {
            if abs_val > i128::MAX as u128 {
                None
            } else {
                Some(abs_val as i128)
            }
        }
    }

    /// 转换为128位无符号整数
    pub fn to_u128(&self) -> Option<u128> {
        if self.negative || self.limbs.len() > 2 {
            return None;
        }

        let result = if self.limbs.len() == 1 {
            self.limbs[0] as u128
        } else {
            (self.limbs[1] as u128) << 64 | (self.limbs[0] as u128)
        };

        Some(result)
    }

    /// 是否为零
    pub fn is_zero(&self) -> bool {
        self.limbs.len() == 1 && self.limbs[0] == 0
    }

    /// 是否为正数
    pub fn is_positive(&self) -> bool {
        !self.negative && !self.is_zero()
    }

    /// 是否为负数
    pub fn is_negative(&self) -> bool {
        self.negative
    }

    /// 获取绝对值
    pub fn abs(&self) -> Self {
        let mut result = self.clone();
        result.negative = false;
        result
    }

    /// 获取符号 (-1, 0, 1)
    pub fn sign(&self) -> i32 {
        if self.is_zero() {
            0
        } else if self.negative {
            -1
        } else {
            1
        }
    }

    /// 获取位长度
    pub fn bit_length(&self) -> usize {
        if self.is_zero() {
            return 0;
        }

        let highest_limb = self.limbs[self.limbs.len() - 1];
        let limb_bits = 64 - highest_limb.leading_zeros() as usize;
        (self.limbs.len() - 1) * 64 + limb_bits
    }

    /// 标准化：移除前导零
    pub(crate) fn normalize(&mut self) {
        while self.limbs.len() > 1 && self.limbs[self.limbs.len() - 1] == 0 {
            self.limbs.pop();
        }

        // 如果结果为零，确保符号为正
        if self.is_zero() {
            self.negative = false;
        }
    }

    /// 获取limb数量
    pub fn limb_count(&self) -> usize {
        self.limbs.len()
    }

    /// 获取指定位置的limb
    pub fn get_limb(&self, index: usize) -> u64 {
        self.limbs.get(index).copied().unwrap_or(0)
    }

    /// 获取limbs的引用（内部使用）
    pub(crate) fn limbs(&self) -> &[u64] {
        &self.limbs
    }

    /// 获取可变limbs的引用（内部使用）
    pub(crate) fn limbs_mut(&mut self) -> &mut Vec<u64> {
        &mut self.limbs
    }

    /// 获取符号（内部使用）
    #[allow(dead_code)]
    pub(crate) fn get_negative(&self) -> bool {
        self.negative
    }

    /// 设置符号
    pub(crate) fn set_negative(&mut self, negative: bool) {
        self.negative = negative;
        if self.is_zero() {
            self.negative = false;
        }
    }
}

impl Default for Mpz {
    fn default() -> Self {
        Self::new()
    }
}

impl PartialEq for Mpz {
    fn eq(&self, other: &Self) -> bool {
        // 零的特殊处理
        if self.is_zero() && other.is_zero() {
            return true;
        }

        self.negative == other.negative && self.limbs == other.limbs
    }
}

impl Eq for Mpz {}

impl std::fmt::Display for Mpz {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string(10))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_creation() {
        let zero = Mpz::new();
        assert!(zero.is_zero());
        assert!(!zero.is_positive());
        assert!(!zero.is_negative());

        let pos = Mpz::from_i64(42);
        assert!(!pos.is_zero());
        assert!(pos.is_positive());
        assert!(!pos.is_negative());

        let neg = Mpz::from_i64(-42);
        assert!(!neg.is_zero());
        assert!(!neg.is_positive());
        assert!(neg.is_negative());
    }

    #[test]
    fn test_conversion() {
        let n = Mpz::from_i64(12345);
        assert_eq!(n.to_i64(), Some(12345));
        assert_eq!(n.to_u64(), Some(12345));

        let neg = Mpz::from_i64(-12345);
        assert_eq!(neg.to_i64(), Some(-12345));
        assert_eq!(neg.to_u64(), None);
    }

    #[test]
    fn test_from_string() {
        let n = Mpz::from_str("12345", 10).unwrap();
        assert_eq!(n.to_i64(), Some(12345));

        let hex = Mpz::from_str("0xFF", 16).unwrap();
        assert_eq!(hex.to_i64(), Some(255));

        let bin = Mpz::from_str("0b1010", 2).unwrap();
        assert_eq!(bin.to_i64(), Some(10));
    }

    #[test]
    fn test_to_string() {
        let n = Mpz::from_i64(255);
        assert_eq!(n.to_string(10), "255");
        assert_eq!(n.to_string(16), "ff");
        assert_eq!(n.to_string(2), "11111111");
    }

    #[test]
    fn test_large_numbers() {
        // 测试大数值的字符串解析
        let large_str = "3141592653589793";
        let large_mpz = Mpz::from_str(large_str, 10).unwrap();

        println!("测试大数值解析:");
        println!("  输入字符串: {}", large_str);
        println!("  解析结果: {}", large_mpz);
        println!("  位数: {}", large_mpz.bit_length());
        println!("  能否转换为 u64: {:?}", large_mpz.to_u64());
        println!("  能否转换为 u128: {:?}", large_mpz.to_u128());
        println!("  能否转换为 i64: {:?}", large_mpz.to_i64());
        println!("  能否转换为 i128: {:?}", large_mpz.to_i128());

        // 验证解析结果
        assert_eq!(large_mpz.to_string(10), large_str);
        assert_eq!(large_mpz.bit_length(), 52); // 3141592653589793 的二进制位数

        // 测试边界情况
        let max_u64_str = "18446744073709551615"; // u64::MAX
        let max_u64_mpz = Mpz::from_str(max_u64_str, 10).unwrap();
        assert_eq!(max_u64_mpz.to_u64(), Some(u64::MAX));

        let beyond_u64_str = "18446744073709551616"; // u64::MAX + 1
        let beyond_u64_mpz = Mpz::from_str(beyond_u64_str, 10).unwrap();
        assert_eq!(beyond_u64_mpz.to_u64(), None); // 应该无法转换为 u64
        assert_eq!(beyond_u64_mpz.to_u128(), Some(18446744073709551616));
    }

    #[test]
    fn test_pi_related_numbers() {
        // 测试与 π 相关的数值
        let pi_mantissa_str = "3141592653589793";
        let pi_mantissa = Mpz::from_str(pi_mantissa_str, 10).unwrap();

        println!("测试 π 相关数值:");
        println!("  π 尾数字符串: {}", pi_mantissa_str);
        println!("  解析结果: {}", pi_mantissa);
        println!("  位数: {}", pi_mantissa.bit_length());
        println!("  能否转换为 u64: {:?}", pi_mantissa.to_u64());
        println!("  能否转换为 u128: {:?}", pi_mantissa.to_u128());

        // 验证 π 尾数的正确性
        assert_eq!(pi_mantissa.to_string(10), pi_mantissa_str);

        // 测试 π/4 相关的数值
        let pi_over_4_mantissa_str = "785398163397448";
        let pi_over_4_mantissa = Mpz::from_str(pi_over_4_mantissa_str, 10).unwrap();

        println!("  π/4 尾数字符串: {}", pi_over_4_mantissa_str);
        println!("  解析结果: {}", pi_over_4_mantissa);
        println!("  位数: {}", pi_over_4_mantissa.bit_length());
        println!("  能否转换为 u64: {:?}", pi_over_4_mantissa.to_u64());

        assert_eq!(pi_over_4_mantissa.to_string(10), pi_over_4_mantissa_str);
    }

    #[test]
    fn test_arithmetic_operations() {
        // 测试大数值的算术运算
        let a = Mpz::from_str("3141592653589793", 10).unwrap();
        let b = Mpz::from_str("1000000000000000", 10).unwrap();

        let sum = a.add(&b);
        let expected_sum = "4141592653589793";
        assert_eq!(sum.to_string(10), expected_sum);

        let diff = a.sub(&b);
        let expected_diff = "2141592653589793";
        assert_eq!(diff.to_string(10), expected_diff);

        let product = a.mul(&b);
        let expected_product = "3141592653589793000000000000000";
        assert_eq!(product.to_string(10), expected_product);
    }

    #[test]
    fn test_edge_cases() {
        // 测试边界情况
        let zero = Mpz::from_str("0", 10).unwrap();
        assert!(zero.is_zero());
        assert_eq!(zero.to_u64(), Some(0));

        let one = Mpz::from_str("1", 10).unwrap();
        assert_eq!(one.to_u64(), Some(1));
        assert_eq!(one.bit_length(), 1);

        let minus_one = Mpz::from_str("-1", 10).unwrap();
        assert!(minus_one.is_negative());
        assert_eq!(minus_one.to_i64(), Some(-1));
        assert_eq!(minus_one.to_u64(), None);

        // 测试空字符串和无效输入
        assert!(Mpz::from_str("", 10).is_err());
        assert!(Mpz::from_str("abc", 10).is_err());
        assert!(Mpz::from_str("123abc", 10).is_err());
    }
}
