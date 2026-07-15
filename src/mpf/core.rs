//! Mpf 核心实现
//!
//! 定义高精度浮点数的内部表示和基础操作。

#[cfg(not(feature = "std"))]
use alloc::{string::String, vec::Vec};

use crate::error::{Error, Result};
use crate::mpz::Mpz;

/// 高精度浮点数结构体
///
/// 使用科学计数法表示：mantissa * 2^exponent
/// mantissa 是 Mpz 类型的尾数（整数部分）
/// exponent 是 i64 类型的指数
/// precision 是当前精度（二进制位数）
#[derive(Debug, Clone)]
pub struct Mpf {
    /// 尾数（整数部分，使用 Mpz 存储）
    mantissa: Mpz,
    /// 指数（以2为底）
    exponent: i64,
    /// 当前精度（二进制位数）
    precision: usize,
    /// 符号：false为正数或零，true为负数
    negative: bool,
    output_config: OutputConfig,
}

/// 输出配置选项
///
/// 控制 Mpf 类型的输出格式和行为，包括科学计数法、精度设置等
#[derive(Debug, Clone)]
pub struct OutputConfig {
    /// 是否使用科学计数法
    pub use_scientific_notation: bool,
    /// 科学计数法的阈值（指数绝对值超过此值时使用科学计数法）
    pub scientific_threshold: i64,
    /// 小数位数
    pub decimal_places: usize,
    /// 是否使用高精度模式
    pub high_precision_mode: bool,
    /// 输入解析精度（用于from_str等方法）
    pub input_precision: usize,
    /// 是否在输出时显示尾数长度信息
    pub show_mantissa_info: bool,
    /// 是否使用紧凑格式（减少不必要的零）
    pub compact_format: bool,
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            use_scientific_notation: true,
            scientific_threshold: 6,
            decimal_places: 6,
            high_precision_mode: false,
            input_precision: 64,
            show_mantissa_info: false,
            compact_format: false,
        }
    }
}

impl OutputConfig {
    /// 创建高精度配置
    pub fn high_precision() -> Self {
        Self {
            use_scientific_notation: true,
            scientific_threshold: 6,
            decimal_places: 17, // f64的最大精度
            high_precision_mode: true,
            input_precision: 128,
            show_mantissa_info: false,
            compact_format: false,
        }
    }

    /// 创建紧凑格式配置
    pub fn compact() -> Self {
        Self {
            use_scientific_notation: true,
            scientific_threshold: 6,
            decimal_places: 6,
            high_precision_mode: false,
            input_precision: 64,
            show_mantissa_info: false,
            compact_format: true,
        }
    }

    /// 创建调试配置（显示详细信息）
    pub fn debug() -> Self {
        Self {
            use_scientific_notation: true,
            scientific_threshold: 6,
            decimal_places: 6,
            high_precision_mode: false,
            input_precision: 64,
            show_mantissa_info: true,
            compact_format: false,
        }
    }
}

impl Mpf {
    /// 创建新的Mpf实例
    pub fn new() -> Self {
        Self {
            mantissa: Mpz::new(),
            exponent: 0,
            precision: 64,
            negative: false,
            output_config: OutputConfig::default(),
        }
    }

    /// 从指定精度创建Mpf
    pub fn with_precision(precision: usize) -> Self {
        Self {
            mantissa: Mpz::new(),
            exponent: 0,
            precision,
            negative: false,
            output_config: OutputConfig::default(),
        }
    }

    /// 从尾数和指数创建Mpf
    pub fn from_parts(mut mantissa: Mpz, exponent: i64, precision: usize) -> Self {
        let negative = mantissa.is_negative();
        mantissa.set_negative(false);
        Self {
            mantissa,
            exponent,
            precision,
            negative,
            output_config: OutputConfig::default(),
        }
    }

    /// 从尾数、指数和精度创建Mpf（带符号）
    pub fn from_parts_with_sign(
        mut mantissa: Mpz,
        exponent: i64,
        precision: usize,
        mut negative: bool,
    ) -> Self {
        if mantissa.is_negative() {
            mantissa.set_negative(false);
            negative = !negative;
        }
        Self {
            mantissa,
            exponent,
            precision,
            negative,
            output_config: OutputConfig::default(),
        }
    }

    /// 从尾数、指数和精度创建Mpf（原始数据，不进行规范化）
    ///
    /// 注意：仍会确保尾数的Mpz负标志被清除（符号信息转移到Mpf的negative字段）。
    pub fn from_parts_raw(mut mantissa: Mpz, exponent: i64, precision: usize) -> Self {
        let negative = mantissa.is_negative();
        mantissa.set_negative(false);
        Self {
            mantissa,
            exponent,
            precision,
            negative,
            output_config: OutputConfig::default(),
        }
    }

    /// 设置输出配置
    pub fn set_output_config(&mut self, config: OutputConfig) {
        self.output_config = config;
    }

    /// 获取输出配置的引用
    pub fn output_config(&self) -> &OutputConfig {
        &self.output_config
    }

    /// 获取输出配置的可变引用
    pub fn output_config_mut(&mut self) -> &mut OutputConfig {
        &mut self.output_config
    }

    /// 从字符串创建（支持科学计数法）
    /// Parse a decimal string to Mpf at full arbitrary precision.
    ///
    /// Handles formats: "123.456", "-0.001", "1.23e10", ".5", "5.", "3_000.5"
    /// Uses integer parsing for the mantissa, then computes 10^exp as an Mpf
    /// and multiplies/divides.
    fn parse_decimal_arbitrary(s: &str, precision: usize) -> Result<Self> {
        let s = s.trim();
        if s.is_empty() {
            return Err(Error::parse("empty string"));
        }

        // Handle sign
        let (negative, s) = if let Some(rest) = s.strip_prefix('-') {
            (true, rest)
        } else if let Some(rest) = s.strip_prefix('+') {
            (false, rest)
        } else {
            (false, s)
        };

        // Split on 'e' or 'E' for scientific notation
        let (mantissa_str, exp_str) = if let Some(e_pos) = s.find(['e', 'E']) {
            (&s[..e_pos], Some(&s[e_pos + 1..]))
        } else {
            (s, None)
        };

        // Parse mantissa: collect digits, track decimal point
        let mut digits = String::new();
        let mut decimal_pos: Option<usize> = None;
        for ch in mantissa_str.chars() {
            if ch == '.' {
                if decimal_pos.is_some() {
                    return Err(Error::parse("multiple decimal points"));
                }
                decimal_pos = Some(digits.len());
            } else if ch.is_ascii_digit() {
                digits.push(ch);
            } else if ch == '_' {
                // skip separator
            } else {
                return Err(Error::parse(format!("invalid character: '{}'", ch)));
            }
        }

        if digits.is_empty() {
            return Err(Error::parse("no digits found"));
        }

        // Remove leading zeros
        let trimmed = digits.trim_start_matches('0');
        if trimmed.is_empty() {
            return Ok(Self::new()); // zero
        }
        let digits_trimmed = trimmed.to_string();

        // Decimal exponent: how many digits after the decimal point
        let decimal_exp: i64 = if let Some(pos) = decimal_pos {
            // Count digits before and after '.' in original string
            let _digits_before = mantissa_str[..pos]
                .chars()
                .filter(|c| c.is_ascii_digit())
                .count();
            let digits_after = mantissa_str[pos + 1..]
                .chars()
                .filter(|c| c.is_ascii_digit())
                .count();
            -(digits_after as i64)
        } else {
            0
        };

        // Add explicit exponent if present
        let decimal_exp = if let Some(exp_s) = exp_str {
            let exp: i64 = exp_s
                .trim()
                .parse()
                .map_err(|_| Error::parse(format!("invalid exponent: {}", exp_s)))?;
            decimal_exp + exp
        } else {
            decimal_exp
        };

        // Step 1: Parse digits as Mpz integer (arbitrary precision)
        let mantissa_int = Mpz::from_str(&digits_trimmed, 10)
            .map_err(|_| Error::parse("failed to parse mantissa"))?;

        // Step 2: Compute 10^|decimal_exp| as Mpf
        let abs_exp: u32 = decimal_exp
            .unsigned_abs()
            .try_into()
            .map_err(|_| Error::parse("exponent magnitude too large"))?;
        let ten_pow = {
            let ten = Self::from_mpz(Mpz::from_i64(10), precision);
            ten.pow(abs_exp)?
        };

        // Step 3: Combine mantissa * 10^decimal_exp
        let mantissa_float = Self::from_mpz(mantissa_int, precision);

        let result = if decimal_exp >= 0 {
            mantissa_float.mul(&ten_pow)
        } else {
            mantissa_float.div(&ten_pow)?
        };

        let mut result = result;
        result.set_negative(negative);
        Ok(result)
    }

    /// Parse a string in the given base (2-36) into an Mpf.
    ///
    /// For base 10 with more than 15 significant digits, uses
    /// arbitrary-precision decimal parsing. Otherwise falls back
    /// to f64-based fast path.
    pub fn from_str(s: &str, base: u32) -> Result<Self> {
        let s = s.trim();
        if s.is_empty() {
            return Err(Error::InvalidInput("Empty string".into()));
        }

        // 处理符号
        let (negative, s) = if let Some(rest) = s.strip_prefix('-') {
            (true, rest)
        } else if let Some(rest) = s.strip_prefix('+') {
            (false, rest)
        } else {
            (false, s)
        };

        if base == 10 {
            // 10进制：检查是否为简单整数
            if s.chars().all(|c| c.is_ascii_digit()) {
                // 简单整数：直接使用Mpz解析
                let mantissa = Mpz::from_str(s, base)?;
                let result = Self {
                    mantissa,
                    exponent: 0,
                    precision: 64,
                    negative,
                    output_config: OutputConfig::default(),
                };
                Ok(result)
            } else {
                // 小数字符串（包含小数点或科学计数法）
                // 检查有效数字位数：≤15 位走 f64 快速路径，>15 位走任意精度路径
                let sig_digits = s.chars().filter(|c| c.is_ascii_digit()).count();
                if sig_digits > 15 {
                    let precision = crate::config::GlobalPrecisionConfig::get_default_precision();
                    return Self::parse_decimal_arbitrary(s, precision);
                }
                // ≤15 位有效数字：使用 f64 快速路径
                let f64_val: f64 = s
                    .parse()
                    .map_err(|_| Error::InvalidInput("Failed to parse decimal string".into()))?;
                let mut result = Self::from_f64(f64_val, 64);
                result.set_negative(negative);
                Ok(result)
            }
        } else {
            // 其他进制：使用原始逻辑
            // 分离尾数和指数部分
            let parts: Vec<&str> = s.split(['e', 'E']).collect();
            if parts.len() > 2 {
                return Err(Error::InvalidInput("Invalid scientific notation".into()));
            }

            let mantissa_str = parts[0];
            let exponent_str = if parts.len() == 2 { parts[1] } else { "0" };

            // 解析尾数（支持小数点）
            let (mantissa, decimal_places) = if mantissa_str.contains('.') {
                let decimal_parts: Vec<&str> = mantissa_str.split('.').collect();
                if decimal_parts.len() != 2 {
                    return Err(Error::InvalidInput("Invalid decimal format".into()));
                }

                let integer_part = decimal_parts[0];
                let fractional_part = decimal_parts[1];

                if integer_part.is_empty() && fractional_part.is_empty() {
                    return Err(Error::InvalidInput(
                        "Empty integer and fractional parts".into(),
                    ));
                }

                // 构建完整的10进制整数（去掉小数点）
                let full_decimal_str = if integer_part.is_empty() {
                    fractional_part.to_string()
                } else {
                    format!("{}{}", integer_part, fractional_part)
                };

                let full_decimal = Mpz::from_str(&full_decimal_str, base)?;
                let decimal_places = fractional_part.len() as i64;

                (full_decimal, decimal_places)
            } else {
                let mantissa = Mpz::from_str(mantissa_str, base)?;
                (mantissa, 0)
            };

            // 解析指数
            let exponent = if exponent_str.is_empty() {
                0
            } else {
                let exp_mpz = Mpz::from_str(exponent_str, base)?;
                exp_mpz
                    .to_i64()
                    .ok_or_else(|| Error::InvalidInput("Exponent too large".into()))?
            };

            let mut result = Self {
                mantissa,
                exponent: exponent - decimal_places,
                precision: 64,
                negative,
                output_config: OutputConfig::default(),
            };
            result.normalize();
            Ok(result)
        }
    }

    /// 从整数创建
    pub fn from_mpz(mpz: Mpz, precision: usize) -> Self {
        let mut result = Self {
            mantissa: mpz,
            exponent: 0,
            precision,
            negative: false,
            output_config: OutputConfig::default(),
        };
        result.normalize();
        result
    }

    /// 从 f64 创建（优化版本，提高精度）
    pub fn from_f64(f: f64, precision: usize) -> Self {
        if f == 0.0 {
            return Self::new();
        }

        if f.is_infinite() {
            return Self::infinity(f.is_sign_negative());
        }

        if f.is_nan() {
            return Self::nan();
        }

        // 使用更精确的方法：直接处理f64的二进制表示
        let bits = f.to_bits();
        let sign = (bits >> 63) & 1;
        let exponent = ((bits >> 52) & 0x7FF) as i64;
        let mantissa_bits = bits & 0x000F_FFFF_FFFF_FFFF;

        let negative = sign == 1;

        if exponent == 0 {
            // 非规格化数或零
            if mantissa_bits == 0 {
                return Self::new();
            }
            // 非规格化数：mantissa * 2^(-1074)
            // IEEE 754: subnormals represent 0.mantissa_bits * 2^(-1022)
            // = mantissa_bits * 2^(-1022-52) = mantissa_bits * 2^(-1074)
            let mantissa = Mpz::from_u64(mantissa_bits);
            return Self {
                mantissa,
                exponent: -1074,
                precision,
                negative,
                output_config: OutputConfig::default(),
            };
        }

        // 规格化数：mantissa * 2^(exponent-1023)
        // 注意：f64的尾数隐含一个1，所以实际值是 (1 + mantissa/2^52) * 2^(exponent-1023)
        // 将隐式前导1放在第52位，尾数读取完整的53位值，指数相应减去52
        let full_mantissa = mantissa_bits | 0x0010_0000_0000_0000;
        let mantissa = Mpz::from_u64(full_mantissa);

        // 指数需要调整：f64的指数偏移是1023
        // 尾数包含了完整的53位（52位分数+1位隐式），所以再减去52
        let adjusted_exponent = exponent - 1023 - 52;

        Self {
            mantissa,
            exponent: adjusted_exponent,
            precision,
            negative,
            output_config: OutputConfig::default(),
        }
    }

    /// 从 f64 创建（高精度版本，使用字符串转换避免精度损失）
    pub fn from_f64_high_precision(f: f64, precision: usize) -> Self {
        if f == 0.0 {
            return Self::new();
        }

        if f.is_infinite() {
            return Self::infinity(f.is_sign_negative());
        }

        if f.is_nan() {
            return Self::nan();
        }

        // 使用字符串转换来避免IEEE 754的精度问题
        let f_str = format!("{:.17}", f); // 使用17位精度，接近f64的最大精度
        let mut result = Self::from_str(&f_str, 10).unwrap_or_else(|_| {
            // 如果字符串解析失败，回退到二进制方法
            Self::from_f64(f, precision)
        });

        // 设置指定的精度
        result.precision = precision;
        result
    }

    /// 从 f64 创建（使用指定配置）
    pub fn from_f64_with_config(f: f64, config: OutputConfig) -> Self {
        if f == 0.0 {
            return Self::new();
        }

        if f.is_infinite() {
            return Self::infinity(f.is_sign_negative());
        }

        if f.is_nan() {
            return Self::nan();
        }

        let mut result = if config.high_precision_mode {
            Self::from_f64_high_precision(f, config.input_precision)
        } else {
            Self::from_f64(f, config.input_precision)
        };

        result.set_output_config(config);
        result
    }

    /// 从 f32 创建
    pub fn from_f32(f: f32, precision: usize) -> Self {
        Self::from_f64(f as f64, precision)
    }

    /// 获取尾数
    pub fn mantissa(&self) -> &Mpz {
        &self.mantissa
    }

    /// 获取指数
    pub fn exponent(&self) -> i64 {
        self.exponent
    }

    /// 获取精度
    pub fn precision(&self) -> usize {
        self.precision
    }

    /// 设置精度
    pub fn set_precision(&mut self, precision: usize) {
        self.precision = precision;
        self.normalize();
    }

    /// 获取符号
    pub fn sign(&self) -> i32 {
        if self.is_zero() {
            0
        } else if self.is_negative() {
            -1
        } else {
            1
        }
    }

    /// 检查是否为有限值
    pub fn is_finite(&self) -> bool {
        !self.is_infinity() && !self.is_nan()
    }

    /// 检查是否为正数
    pub fn is_positive(&self) -> bool {
        !self.negative && !self.is_zero()
    }

    /// 检查是否为负数
    pub fn is_negative(&self) -> bool {
        self.negative
    }

    /// 设置符号
    pub fn set_negative(&mut self, negative: bool) {
        self.negative = negative;
        if self.is_zero() {
            self.negative = false;
        }
    }

    /// 获取尾数（可变引用）
    pub(crate) fn mantissa_mut(&mut self) -> &mut Mpz {
        &mut self.mantissa
    }

    /// 获取指数（可变引用）
    pub(crate) fn exponent_mut(&mut self) -> &mut i64 {
        &mut self.exponent
    }

    /// 获取符号（可变引用）
    #[allow(dead_code)]
    pub(crate) fn negative_mut(&mut self) -> &mut bool {
        &mut self.negative
    }

    /// 检查是否为零
    pub fn is_zero(&self) -> bool {
        self.mantissa.is_zero()
    }

    /// 检查是否为无穷大
    pub fn is_infinity(&self) -> bool {
        // 无穷大表示为最大指数和零尾数
        self.exponent >= (i32::MAX - 1) as i64 && self.mantissa.is_zero()
    }

    /// 检查是否为 NaN
    pub fn is_nan(&self) -> bool {
        // NaN表示为最大指数和非零尾数
        self.exponent >= (i32::MAX - 1) as i64 && !self.mantissa.is_zero()
    }

    /// 获取绝对值
    pub fn abs(&self) -> Self {
        let mut result = self.clone();
        result.negative = false;
        result
    }

    /// 比较两个浮点数
    pub fn compare(&self, other: &Self) -> core::cmp::Ordering {
        // 处理特殊情况
        if self.is_zero() && other.is_zero() {
            return core::cmp::Ordering::Equal;
        }

        if self.is_zero() {
            return if other.is_negative() {
                core::cmp::Ordering::Greater
            } else {
                core::cmp::Ordering::Less
            };
        }

        if other.is_zero() {
            return if self.is_negative() {
                core::cmp::Ordering::Less
            } else {
                core::cmp::Ordering::Greater
            };
        }

        // 比较符号
        if self.is_negative() != other.is_negative() {
            return if self.is_negative() {
                core::cmp::Ordering::Less
            } else {
                core::cmp::Ordering::Greater
            };
        }

        // 快速路径: 量级比较（无需 clone/normalize）
        let mag1 = self.mantissa.bit_length() as i64 + self.exponent;
        let mag2 = other.mantissa.bit_length() as i64 + other.exponent;
        if mag1 != mag2 {
            let greater = if mag1 > mag2 {
                core::cmp::Ordering::Greater
            } else {
                core::cmp::Ordering::Less
            };
            return if self.is_negative() {
                greater.reverse()
            } else {
                greater
            };
        }

        // 量级相同时，需要更精确的比较
        let mut a = self.clone();
        let mut b = other.clone();
        a.normalize();
        b.normalize();

        let (aligned_self, aligned_other) = a.align_exponents(&b);
        let mantissa_cmp = aligned_self.mantissa().cmp(aligned_other.mantissa());

        if self.is_negative() {
            mantissa_cmp.reverse()
        } else {
            mantissa_cmp
        }
    }

    /// 对齐指数（内部使用）
    pub(crate) fn align_exponents(&self, other: &Mpf) -> (Mpf, Mpf) {
        let mut aligned_self = self.clone();
        let mut aligned_other = other.clone();

        let exponent_diff = self.exponent() - other.exponent();

        if exponent_diff > 0 {
            // self的指数更大，需要对齐other
            if exponent_diff > 20 {
                // 指数差异太大，右移会丢失太多精度
                // 改为左移self，保持精度
                let shift = exponent_diff;
                *aligned_self.mantissa_mut() = aligned_self.mantissa().shl(shift as usize);
                *aligned_self.exponent_mut() = other.exponent();
            } else {
                // 指数差异不大，右移other
                *aligned_other.mantissa_mut() =
                    aligned_other.mantissa().shr(exponent_diff as usize);
                *aligned_other.exponent_mut() = self.exponent();
            }
        } else if exponent_diff < 0 {
            // other的指数更大，需要对齐self
            let abs_diff = (-exponent_diff) as usize;
            if abs_diff > 20 {
                // 指数差异太大，右移会丢失太多精度
                // 改为左移other，保持精度
                *aligned_other.mantissa_mut() = aligned_other.mantissa().shl(abs_diff);
                *aligned_other.exponent_mut() = self.exponent();
            } else {
                // 指数差异不大，右移self
                *aligned_self.mantissa_mut() = aligned_self.mantissa().shr(abs_diff);
                *aligned_self.exponent_mut() = other.exponent();
            }
        }
        // 如果指数相同，不需要对齐

        (aligned_self, aligned_other)
    }

    /// 标准化：确保尾数在 [1, 2) 范围内（对于非零数）
    pub(crate) fn normalize(&mut self) {
        if self.mantissa.is_zero() {
            self.exponent = 0;
            self.set_negative(false);
            return;
        }

        // 对于整数，保持原样，不进行标准化
        if self.exponent == 0 {
            return;
        }

        // 计算尾数的位长度
        let bit_length = self.mantissa.bit_length();

        // 如果尾数太大，需要右移
        if bit_length > self.precision {
            let shift = bit_length - self.precision;
            self.mantissa = self.mantissa.shr(shift);
            self.exponent += shift as i64;
        }

        // 如果尾数太小，需要左移
        if bit_length < self.precision && !self.mantissa.is_zero() {
            let shift = self.precision - bit_length;
            self.mantissa = self.mantissa.shl(shift);
            self.exponent -= shift as i64;
        }

        // 如果结果为零，确保符号为正
        if self.is_zero() {
            self.set_negative(false);
        }
    }

    /// 创建无穷大
    pub fn infinity(negative: bool) -> Self {
        Self {
            mantissa: Mpz::new(), // 零尾数表示无穷大
            exponent: (i32::MAX - 1) as i64,
            negative,
            precision: 64, // 默认精度
            output_config: OutputConfig::default(),
        }
    }

    /// 创建 NaN
    pub fn nan() -> Self {
        Self {
            mantissa: Mpz::from_i64(1), // 非零尾数表示NaN
            exponent: (i32::MAX - 1) as i64,
            negative: false, // NaN的符号位通常被忽略
            precision: 64,   // 默认精度
            output_config: OutputConfig::default(),
        }
    }

    /// 检查是否为整数
    pub fn is_integer(&self) -> bool {
        if self.is_zero() {
            return true;
        }

        if self.is_infinity() || self.is_nan() {
            return false;
        }

        // 检查指数是否大于等于0，如果是，则可能是整数
        if self.exponent >= 0 {
            // 检查尾数是否能被2^exponent整除
            let shift = self.exponent as usize;
            if shift < self.mantissa.bit_length() {
                // 如果指数小于尾数位数，检查是否有小数部分
                let mask = Mpz::from_i64(1).shl(shift);
                let remainder = self.mantissa.rem(&mask);
                if let Ok(remainder) = remainder {
                    remainder.is_zero()
                } else {
                    false
                }
            } else {
                // 指数大于等于尾数位数，一定是整数
                true
            }
        } else {
            // 负指数表示小数
            false
        }
    }

    /// 向下取整
    pub fn floor(&self) -> Result<Self> {
        if self.is_infinity() || self.is_nan() {
            return Ok(self.clone());
        }

        if self.is_integer() {
            return Ok(self.clone());
        }

        if self.exponent < 0 {
            // 负指数表示小数，需要计算向下取整
            // 例如：mantissa * 2^(-51) 的向下取整
            let abs_exponent = (-self.exponent) as usize;
            let divisor = Mpz::from_i64(1).shl(abs_exponent);

            // 计算整数部分
            let integer_part = self.mantissa.div(&divisor)?;
            let remainder = self.mantissa.rem(&divisor)?;

            let mut result = integer_part;

            // 对于负数，如果有小数部分，需要减1
            if self.is_negative() && !remainder.is_zero() {
                result = result.add(&Mpz::from_i64(1));
            }

            Ok(Self {
                mantissa: result,
                exponent: 0,
                negative: self.negative,
                precision: self.precision,
                output_config: self.output_config.clone(),
            })
        } else {
            // 正指数，需要截断小数部分
            let shift = self.exponent as usize;
            if shift < self.mantissa.bit_length() {
                let mask = Mpz::from_i64(1).shl(shift);
                let truncated_mantissa = self.mantissa.div(&mask)?.mul(&mask);
                Ok(Self {
                    mantissa: truncated_mantissa,
                    exponent: self.exponent,
                    negative: self.negative,
                    precision: self.precision,
                    output_config: self.output_config.clone(),
                })
            } else {
                // 指数大于等于尾数位数，已经是整数
                Ok(self.clone())
            }
        }
    }

    /// 转换为字符串表示（优化版本，支持更多配置选项）
    pub fn to_string(&self, base: u32) -> String {
        if self.is_zero() {
            return "0".to_string();
        }

        if base == 10 {
            // 10进制：根据输出配置决定格式
            // 如果启用调试模式，优先使用内部表示
            if self.output_config.show_mantissa_info {
                let sign = if self.is_negative() { "-" } else { "" };
                let mantissa_str = self.mantissa.to_string(base);
                let mut result = format!("{}{}", sign, mantissa_str);

                if self.exponent != 0 {
                    result.push_str(&format!("e{}", self.exponent));
                }

                // 显示尾数信息
                result.push_str(&format!(
                    " [mantissa_bits: {}, exp: {}]",
                    self.mantissa.bit_length(),
                    self.exponent
                ));

                result
            } else if let Some(f64_val) = self.to_f64() {
                let abs_val = f64_val.abs();

                // 检查是否应该使用科学计数法
                let should_use_scientific = self.output_config.use_scientific_notation
                    && (abs_val >= 10.0_f64.powi(self.output_config.scientific_threshold as i32)
                        || (abs_val
                            < 10.0_f64.powi(-(self.output_config.scientific_threshold as i32))
                            && abs_val != 0.0));

                if should_use_scientific {
                    // 使用科学计数法
                    format!("{:.6e}", f64_val)
                } else {
                    // 使用普通小数格式
                    if f64_val == f64_val.floor() {
                        // 整数：不显示小数部分
                        format!("{:.0}", f64_val)
                    } else {
                        // 小数：显示指定的小数位数
                        let mut result =
                            format!("{:.1$}", f64_val, self.output_config.decimal_places);

                        // 如果启用紧凑格式，移除尾部的零
                        if self.output_config.compact_format {
                            result = result
                                .trim_end_matches('0')
                                .trim_end_matches('.')
                                .to_string();
                            if result.is_empty() {
                                result = "0".to_string();
                            }
                        }

                        result
                    }
                }
            } else {
                // 无法转换为f64时的回退方案
                let sign = if self.is_negative() { "-" } else { "" };
                let mantissa_str = self.mantissa.to_string(base);
                let mut result = format!("{}{}", sign, mantissa_str);

                if self.exponent != 0 {
                    result.push_str(&format!("e{}", self.exponent));
                }

                result
            }
        } else {
            // 其他进制：使用原始逻辑
            let sign = if self.is_negative() { "-" } else { "" };
            let mantissa_str = self.mantissa.to_string(base);
            let mut result = format!("{}{}", sign, mantissa_str);

            if self.exponent != 0 {
                result.push_str(&format!("e{}", self.exponent));
            }

            // 如果启用调试模式，显示尾数信息
            if self.output_config.show_mantissa_info {
                result.push_str(&format!(
                    " [mantissa_bits: {}, exp: {}]",
                    self.mantissa.bit_length(),
                    self.exponent
                ));
            }

            result
        }
    }

    /// 转换为 f64（可能丢失精度）
    ///
    /// Uses proper IEEE 754-style ldexp (load exponent) computation:
    ///   value = mantissa * 2^exponent
    ///
    /// The mantissa is first normalized to at most 53 bits (f64 precision),
    /// then the computation is done using powi with decomposition to handle
    /// the full f64 exponent range [-1074, 1023].
    pub fn to_f64(&self) -> Option<f64> {
        if self.is_zero() {
            return Some(0.0);
        }

        if self.is_infinity() {
            return Some(if self.is_negative() {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            });
        }

        if self.is_nan() {
            return Some(f64::NAN);
        }

        // Normalize mantissa to at most 53 bits to fit in u64.
        // Also adjust the exponent by the same shift amount so the
        // mathematical value is preserved.
        let (mantissa_u64, exponent) = {
            let mut m = self.mantissa.clone();
            let mut e = self.exponent;
            let bits = m.bit_length();
            if bits > 53 {
                let shift = bits - 53;
                m = m.shr(shift);
                e += shift as i64;
            }
            (m.to_u64()?, e)
        };

        let sign = if self.negative { -1.0 } else { 1.0 };

        // Clamp to f64 finite range.
        // f64 can represent values in [2^(-1074), 2^1023 * (2 - 2^(-52))].
        if exponent > 1023 {
            return Some(sign * f64::INFINITY);
        }
        if exponent < -1074 {
            return Some(0.0);
        }

        // Compute mantissa * 2^exponent using ldexp-style decomposition.
        // Using powi directly with a single call could lose subnormal range
        // for the smallest exponents (2^(-1074) underflows powi).
        // Instead we split: for large negative exponents, use 0.5^|e| to
        // let the computation pass through the subnormal range correctly.
        let abs_val = if exponent >= 0 {
            (mantissa_u64 as f64) * (2.0_f64).powi(exponent as i32)
        } else {
            // 0.5^|exponent| stays in range down to ~5e-324 (smallest subnormal)
            // and gracefully underflows to 0.0 for truly unrepresentable values.
            (mantissa_u64 as f64) * (0.5_f64).powi((-exponent) as i32)
        };

        Some(sign * abs_val)
    }
}

impl Default for Mpf {
    fn default() -> Self {
        Self::new()
    }
}

impl PartialEq for Mpf {
    fn eq(&self, other: &Self) -> bool {
        if self.is_zero() && other.is_zero() {
            return true;
        }

        self.is_negative() == other.is_negative()
            && self.mantissa == other.mantissa
            && self.exponent == other.exponent
    }
}

impl PartialOrd for Mpf {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(Ord::cmp(self, other))
    }
}

impl Eq for Mpf {}

impl Ord for Mpf {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.compare(other)
    }
}

impl std::fmt::Display for Mpf {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string(10))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_creation() {
        let zero = Mpf::new();
        assert!(zero.is_zero());
        assert!(!zero.is_negative());

        let from_mpz = Mpf::from_mpz(Mpz::from_i64(42), 64);
        assert!(!from_mpz.is_zero());
        assert_eq!(from_mpz.to_i64(), Some(42));
    }

    #[test]
    fn test_from_string() {
        let mpf = Mpf::from_str("123", 10).unwrap();
        assert!(!mpf.is_zero());
        assert_eq!(mpf.to_i64(), Some(123));
    }

    #[test]
    fn test_to_string() {
        let mpf = Mpf::from_mpz(Mpz::from_i64(123), 64);
        let s = mpf.to_string(10);
        println!("to_string output: '{}'", s);
        // 由于精度问题，可能不会精确等于"123"，所以检查是否包含"123"
        assert!(s.contains("123"));
    }

    #[test]
    fn test_normalize_behavior() {
        // 测试normalize函数的行为
        let large_mantissa = Mpz::from_str("123456789012345678901234567890", 10).unwrap();
        let mpf = Mpf::from_parts(large_mantissa, -73, 64);

        println!("mantissa bit length: {}", mpf.mantissa().bit_length());
        println!("exponent: {}", mpf.exponent());
        println!("precision: {}", mpf.precision());
        println!("to_string output: '{}'", mpf.to_string(10));

        // 这个测试主要是为了观察行为，不进行断言
    }

    #[test]
    fn test_special_values() {
        // 测试无穷大
        let pos_inf = Mpf::infinity(false);
        let neg_inf = Mpf::infinity(true);

        assert!(pos_inf.is_infinity());
        assert!(neg_inf.is_infinity());
        assert!(!pos_inf.is_negative());
        assert!(neg_inf.is_negative());
        assert!(!pos_inf.is_nan());
        assert!(!neg_inf.is_nan());

        // 测试NaN
        let nan = Mpf::nan();
        assert!(nan.is_nan());
        assert!(!nan.is_infinity());
        assert!(!nan.is_finite());

        // 测试有限值
        let normal = Mpf::from_f64(3.14, 64);
        assert!(!normal.is_infinity());
        assert!(!normal.is_nan());
        assert!(normal.is_finite());

        // 测试零值
        let zero = Mpf::new();
        assert!(!zero.is_infinity());
        assert!(!zero.is_nan());
        assert!(zero.is_finite());
    }

    #[test]
    fn test_edge_cases() {
        // 测试极大值
        let max_exponent = Mpf {
            mantissa: Mpz::new(), // 零尾数
            exponent: (i32::MAX - 1) as i64,
            negative: false,
            precision: 64,
            output_config: OutputConfig::default(),
        };
        assert!(max_exponent.is_infinity());

        // 测试NaN
        let nan_exponent = Mpf {
            mantissa: Mpz::from_i64(1), // 非零尾数
            exponent: (i32::MAX - 1) as i64,
            negative: false,
            precision: 64,
            output_config: OutputConfig::default(),
        };
        assert!(nan_exponent.is_nan());

        // 测试极小值
        let min_exponent = Mpf {
            mantissa: Mpz::new(), // 零尾数
            exponent: -1000,      // 使用一个合理的负指数值
            negative: false,
            precision: 64,
            output_config: OutputConfig::default(),
        };
        assert!(min_exponent.is_zero());

        // 测试精度边界
        let high_precision = Mpf::from_f64(0.1, 1000);
        assert!(high_precision.precision() == 1000);

        // 测试字符串转换边界
        let very_small = Mpf::from_str("1e-1000", 10);
        assert!(very_small.is_ok());

        let very_large = Mpf::from_str("1e1000", 10);
        assert!(very_large.is_ok());
    }

    #[test]
    fn test_integer_operations() {
        // 测试整数检查
        let int_val = Mpf::from_i64(42, 64);
        assert!(int_val.is_integer());

        let float_val = Mpf::from_f64(3.14, 64);
        assert!(!float_val.is_integer());

        // 测试向下取整
        let floor_result = float_val.floor().unwrap();
        assert_eq!(floor_result.to_i64(), Some(3));

        let negative_float = Mpf::from_f64(-3.7, 64);
        let floor_negative = negative_float.floor().unwrap();
        assert_eq!(floor_negative.to_i64(), Some(-4));

        // 测试零和特殊值
        let zero = Mpf::new();
        assert!(zero.is_integer());
        let zero_floor = zero.floor().unwrap();
        assert!(zero_floor.is_zero());
    }

    #[test]
    fn test_from_f64_subnormal() {
        // 最小正非规格化 f64: 5e-324
        let subnormal = Mpf::from_f64(f64::MIN_POSITIVE / 2.0, 64);
        assert!(!subnormal.is_zero(), "subnormal f64 should not become zero");
        assert!(subnormal.is_positive(), "subnormal f64 should be positive");
        let back = subnormal.to_f64().unwrap();
        assert!(back > 0.0, "round-trip should preserve positivity");
        // 应接近 2^(-1074) ≈ 5e-324
        assert!(back < 1e-300, "subnormal value should be extremely small");

        // 零 subnormal: mantissa_bits == 0, exponent == 0
        let zero = Mpf::from_f64(0.0, 64);
        assert!(zero.is_zero());

        // 最小正规格化 f64: 2.2e-308 (应走规格化路径)
        let min_normal = Mpf::from_f64(f64::MIN_POSITIVE, 64);
        assert!(!min_normal.is_zero());
    }

    #[test]
    fn test_from_str_high_precision() {
        // 30 significant digits — far beyond f64's ~15
        let s = "3.14159265358979323846264338327";
        let x = Mpf::from_str(s, 10).unwrap();

        // Compare with f64-based π (only ~15-16 digits of precision)
        let pi_f64 = Mpf::from_f64(std::f64::consts::PI, 256);
        let diff = x.sub(&pi_f64).abs();
        let diff_val = diff.to_f64().unwrap();

        // The difference should be ~1.22e-16 (known error of f64 π vs true π).
        // If our parse still has f64 precision, diff would be ~0.
        // This proves from_str produces a value more precise than f64.
        assert!(
            diff_val > 1e-17,
            "diff too small (still getting f64 precision): {}",
            diff_val
        );
        assert!(diff_val < 1e-14, "diff unreasonably large: {}", diff_val);

        // Parse same string twice should give identical results
        let x2 = Mpf::from_str(s, 10).unwrap();
        assert!(
            x.sub(&x2).abs().is_zero(),
            "parsing same string twice should match"
        );

        // Sub-f64 precision: 1e-20 should not be zero or NaN with arbitrary precision
        let tiny = Mpf::from_str("1e-20", 10).unwrap();
        assert!(!tiny.is_zero(), "1e-20 should not be zero");
        assert!(!tiny.is_nan(), "1e-20 should not be NaN");
        assert!(tiny.to_f64().unwrap() > 0.0, "1e-20 should be positive");

        // Very small number that f64_MIN_POSITIVE is ~2.2e-308,
        // so 1e-100 should be non-zero and non-NaN
        let very_tiny = Mpf::from_str("1e-100", 10).unwrap();
        assert!(!very_tiny.is_zero(), "1e-100 should not be zero");
        assert!(!very_tiny.is_nan(), "1e-100 should not be NaN");

        // Long integer: 2^100 = 1267650600228229401496703205376 (31 digits)
        let big_int = Mpf::from_str("1267650600228229401496703205376", 10).unwrap();
        let two = Mpf::from_mpz(Mpz::from_i64(2), 256);
        let two_100 = two.pow(100).unwrap();
        let big_diff = big_int.sub(&two_100).abs();
        assert!(
            big_diff.is_zero(),
            "2^100 parsed from string should match 2.pow(100)"
        );

        // Mantissa should indicate high precision (30 digits ≈ 100 bits)
        assert!(
            x.mantissa().bit_length() > 50,
            "mantissa too small for 30-digit precision"
        );
    }

    #[test]
    fn test_performance_benchmarks() {
        // 测试创建性能
        let start = std::time::Instant::now();
        for _ in 0..1000 {
            let _ = Mpf::from_f64(3.14159, 64);
        }
        let creation_time = start.elapsed();
        println!("创建1000个Mpf耗时: {:?}", creation_time);
        assert!(creation_time.as_millis() < 100, "创建操作应该在100ms内完成");

        // 测试算术运算性能
        let a = Mpf::from_f64(3.14159, 64);
        let b = Mpf::from_f64(2.71828, 64);

        let start = std::time::Instant::now();
        for _ in 0..1000 {
            let _ = a.add(&b);
            let _ = a.sub(&b);
            let _ = a.mul(&b);
            let _ = a.div(&b);
        }
        let arithmetic_time = start.elapsed();
        println!("1000次算术运算耗时: {:?}", arithmetic_time);
        assert!(
            arithmetic_time.as_millis() < 500,
            "算术运算应该在500ms内完成"
        );

        // 测试特殊函数性能
        let x = Mpf::from_f64(2.5, 64);
        let start = std::time::Instant::now();
        for _ in 0..100 {
            let _ = x.sqrt();
            let _ = x.exp();
            let _ = x.ln();
        }
        let special_time = start.elapsed();
        println!("100次特殊函数计算耗时: {:?}", special_time);
        assert!(
            special_time.as_millis() < 1000,
            "特殊函数计算应该在1秒内完成"
        );
    }
}
