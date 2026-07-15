//! 基本算术运算实现
//!
//! Newton's method for integer square root uses iteration control constants
//! shared with [`crate::mpf::arithmetic`] via [`crate::algorithm`].
//!
//! See also: [`crate::mpf::arithmetic`] for the floating-point variant.

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

use super::core::Mpz;
use crate::error::{Error, Result};
use crate::mpz::MultiplicationBackend;

impl Mpz {
    /// 加法运算
    pub fn add(&self, other: &Mpz) -> Mpz {
        // 处理符号
        if self.is_negative() != other.is_negative() {
            // 符号不同，转换为减法
            if self.is_negative() {
                // -a + b = b - a
                other.sub_abs(self)
            } else {
                // a + (-b) = a - b
                self.sub_abs(other)
            }
        } else {
            // 符号相同，进行绝对值加法
            let mut result = self.add_abs(other);
            result.set_negative(self.is_negative());
            result
        }
    }

    /// 减法运算
    pub fn sub(&self, other: &Mpz) -> Mpz {
        // 处理符号
        if self.is_negative() != other.is_negative() {
            // 符号不同，转换为加法
            let mut result = self.add_abs(other);
            result.set_negative(self.is_negative());
            result
        } else {
            // 符号相同，进行绝对值减法
            if self.is_negative() {
                // -a - (-b) = b - a
                other.sub_abs(self)
            } else {
                // a - b
                self.sub_abs(other)
            }
        }
    }

    /// 乘法运算
    pub fn mul(&self, other: &Mpz) -> Mpz {
        if self.is_zero() || other.is_zero() {
            return Mpz::new();
        }

        let result = self.mul_abs(other);
        let negative = self.is_negative() != other.is_negative();

        let mut final_result = result;
        final_result.set_negative(negative);
        final_result
    }

    /// 除法运算
    pub fn div(&self, other: &Mpz) -> Result<Mpz> {
        if other.is_zero() {
            return Err(Error::DivisionByZero);
        }

        if self.is_zero() {
            return Ok(Mpz::new());
        }

        let (quotient, _) = self.div_rem_abs(other)?;
        let negative = self.is_negative() != other.is_negative();

        let mut result = quotient;
        result.set_negative(negative);
        Ok(result)
    }

    /// 取余运算
    pub fn rem(&self, other: &Mpz) -> Result<Mpz> {
        if other.is_zero() {
            return Err(Error::DivisionByZero);
        }

        if self.is_zero() {
            return Ok(Mpz::new());
        }

        let (_, remainder) = self.div_rem_abs(other)?;
        let mut result = remainder;
        result.set_negative(self.is_negative());
        Ok(result)
    }

    /// rem函数的别名，用于取模运算
    pub fn mod_(&self, other: &Mpz) -> Result<Mpz> {
        self.rem(other)
    }

    /// 除法和取余运算
    pub fn div_rem(&self, other: &Mpz) -> Result<(Mpz, Mpz)> {
        if other.is_zero() {
            return Err(Error::DivisionByZero);
        }

        if self.is_zero() {
            return Ok((Mpz::new(), Mpz::new()));
        }

        let (quotient, remainder) = self.div_rem_abs(other)?;
        let negative = self.is_negative() != other.is_negative();

        let mut final_quotient = quotient;
        final_quotient.set_negative(negative);

        let mut final_remainder = remainder;
        final_remainder.set_negative(self.is_negative());

        Ok((final_quotient, final_remainder))
    }

    /// 取负
    pub fn neg(&self) -> Mpz {
        let mut result = self.clone();
        result.set_negative(!result.is_negative());
        result
    }

    /// Raise self to a u32 exponent (fast path for small exponents).
    pub fn pow_u32(&self, exp: u32) -> Mpz {
        if exp == 0 {
            return Mpz::from_i64(1);
        }

        if self.is_zero() {
            return Mpz::new();
        }

        let mut result = Mpz::from_i64(1);
        let mut base = self.clone();
        let mut exponent = exp;

        while exponent > 0 {
            if exponent % 2 == 1 {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            exponent /= 2;
        }

        result
    }

    /// Raise self to an arbitrary-precision integer exponent (binary exponentiation).
    pub fn pow(&self, exp: &Mpz) -> Mpz {
        if exp.is_zero() {
            return Mpz::from_i64(1);
        }
        if self.is_zero() {
            return Mpz::new();
        }
        if exp.is_negative() {
            // 负指数在整数域中未定义
            return Mpz::new();
        }

        let mut result = Mpz::from_i64(1);
        let mut base = self.clone();
        let mut e = exp.clone();
        let two = Mpz::from_i64(2);
        let one = Mpz::from_i64(1);

        while !e.is_zero() {
            if e.rem(&two).unwrap() == one {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            e = e.div(&two).unwrap();
        }
        result
    }

    /// Raise self to an arbitrary-precision exponent.
    ///
    /// Replaced by [`Mpz::pow`]. This function is deprecated.
    #[deprecated(note = "use pow(&exp) instead")]
    pub fn pow_big(&self, exp: &Mpz) -> Mpz {
        self.pow(exp)
    }

    /// 绝对值加法
    fn add_abs(&self, other: &Mpz) -> Mpz {
        let a_limbs = self.limbs();
        let b_limbs = other.limbs();

        let max_len = a_limbs.len().max(b_limbs.len());
        let mut result = Vec::with_capacity(max_len + 1);
        let mut carry = 0u64;

        for i in 0..max_len {
            let a_limb = a_limbs.get(i).copied().unwrap_or(0);
            let b_limb = b_limbs.get(i).copied().unwrap_or(0);

            let (sum1, overflow1) = a_limb.overflowing_add(b_limb);
            let (sum2, overflow2) = sum1.overflowing_add(carry);

            result.push(sum2);
            carry = if overflow1 || overflow2 { 1 } else { 0 };
        }

        if carry != 0 {
            result.push(carry);
        }

        Mpz::from_limbs(result, false)
    }

    /// 绝对值减法
    fn sub_abs(&self, other: &Mpz) -> Mpz {
        let self_cmp = self.cmp_abs(other);

        let (a_limbs, b_limbs, will_be_negative) = match self_cmp {
            core::cmp::Ordering::Greater => (self.limbs(), other.limbs(), false),
            core::cmp::Ordering::Less => (other.limbs(), self.limbs(), true),
            core::cmp::Ordering::Equal => return Mpz::new(),
        };

        let mut result = Vec::with_capacity(a_limbs.len());
        let mut borrow = 0u64;

        for (i, &a_limb) in a_limbs.iter().enumerate() {
            let b_limb = b_limbs.get(i).copied().unwrap_or(0);

            let (diff1, underflow1) = a_limb.overflowing_sub(b_limb);
            let (diff2, underflow2) = diff1.overflowing_sub(borrow);

            result.push(diff2);
            borrow = if underflow1 || underflow2 { 1 } else { 0 };
        }

        let mut final_result = Mpz::from_limbs(result, will_be_negative);
        final_result.normalize();
        final_result
    }

    /// 绝对值乘法（使用自适应算法）
    fn mul_abs(&self, other: &Mpz) -> Mpz {
        // 使用自适应乘法后端，自动选择最优算法
        self.mul_with_backend(other, MultiplicationBackend::Adaptive)
    }

    /// 绝对值除法和取余
    fn div_rem_abs(&self, other: &Mpz) -> Result<(Mpz, Mpz)> {
        let self_limbs = self.limbs();
        let other_limbs = other.limbs();

        // 单limb除法优化
        if other_limbs.len() == 1 && self_limbs.len() == 1 {
            let quotient = self_limbs[0] / other_limbs[0];
            let remainder = self_limbs[0] % other_limbs[0];
            return Ok((Mpz::from_u64(quotient), Mpz::from_u64(remainder)));
        }

        // 使用长除法算法
        self.long_division(other)
    }

    /// 长除法算法（优化版本）
    fn long_division(&self, divisor: &Mpz) -> Result<(Mpz, Mpz)> {
        let dividend = self;
        if divisor.is_zero() {
            return Err(Error::DivisionByZero);
        }

        if dividend.is_zero() {
            return Ok((Mpz::new(), Mpz::new()));
        }

        let dividend_cmp = dividend.cmp_abs(divisor);
        if dividend_cmp == core::cmp::Ordering::Less {
            return Ok((Mpz::new(), dividend.clone()));
        }

        if dividend_cmp == core::cmp::Ordering::Equal {
            return Ok((Mpz::from_i64(1), Mpz::new()));
        }

        // 单limb除法优化
        if divisor.limb_count() == 1 {
            return self.single_limb_division(divisor);
        }

        // 多limb长除法
        self.multi_limb_division(divisor)
    }

    /// 单limb除法（快速路径）
    fn single_limb_division(&self, divisor: &Mpz) -> Result<(Mpz, Mpz)> {
        let divisor_limb = divisor.limbs()[0];
        let mut quotient_limbs = Vec::new();
        let mut remainder = 0u128;

        // 从高位到低位进行除法
        for &dividend_limb in self.limbs().iter().rev() {
            let current = (remainder << 64) | dividend_limb as u128;
            let q = current / divisor_limb as u128;
            remainder = current % divisor_limb as u128;

            if !quotient_limbs.is_empty() || q != 0 {
                quotient_limbs.push(q as u64);
            }
        }

        // 反转结果（因为我们是从高位到低位处理的）
        quotient_limbs.reverse();

        let quotient = if quotient_limbs.is_empty() {
            Mpz::new()
        } else {
            Mpz::from_limbs(quotient_limbs, false)
        };

        let remainder_mpz = if remainder == 0 {
            Mpz::new()
        } else {
            Mpz::from_u64(remainder as u64)
        };

        Ok((quotient, remainder_mpz))
    }

    /// 多limb长除法（Knuth算法）
    fn multi_limb_division(&self, divisor: &Mpz) -> Result<(Mpz, Mpz)> {
        let dividend_limbs = self.limbs();
        let divisor_limbs = divisor.limbs();

        let n = divisor_limbs.len();
        let m = dividend_limbs.len().saturating_sub(n);

        if m == 0 && dividend_limbs.len() < n {
            return Ok((Mpz::new(), self.clone()));
        }

        // 标准化除数（确保最高位limb >= 2^63）
        let (normalized_divisor, shift) = self.normalize_divisor(divisor);
        let normalized_dividend = if shift > 0 {
            self.shl(shift)
        } else {
            self.clone()
        };

        let mut quotient_limbs = vec![0u64; m + 1];
        let mut remainder_limbs = normalized_dividend.limbs().to_vec();

        // 确保remainder_limbs有足够的长度以避免越界
        if remainder_limbs.len() < m + n + 1 {
            remainder_limbs.resize(m + n + 1, 0);
        }

        // 执行长除法
        for j in (0..=m).rev() {
            // 估算商
            let q_hat =
                self.estimate_quotient_digit(&remainder_limbs, normalized_divisor.limbs(), j);

            // 乘以除数并减去
            let mut borrow = 0i128;
            for i in 0..n {
                let product = q_hat as u128 * normalized_divisor.limbs()[i] as u128;
                let current = remainder_limbs[j + i] as u128;
                let diff = current as i128 - (product & 0xFFFFFFFFFFFFFFFF) as i128 - borrow;

                remainder_limbs[j + i] = diff as u64;
                borrow = (product >> 64) as i128 - (diff >> 64);
            }

            // 处理借位
            if borrow > 0 {
                remainder_limbs[j + n] = remainder_limbs[j + n].wrapping_sub(borrow as u64);
            }

            // 如果结果为负，需要修正
            if remainder_limbs[j + n] > normalized_divisor.limbs()[n - 1] {
                self.correct_quotient_digit(&mut remainder_limbs, normalized_divisor.limbs(), j);
            }

            quotient_limbs[j] = q_hat;
        }

        // 构建结果
        let quotient = Mpz::from_limbs(quotient_limbs, false);
        let remainder = if shift > 0 {
            Mpz::from_limbs(remainder_limbs, false).shr(shift)
        } else {
            Mpz::from_limbs(remainder_limbs, false)
        };

        Ok((quotient, remainder))
    }

    /// 标准化除数（确保最高位limb >= 2^63）
    fn normalize_divisor(&self, divisor: &Mpz) -> (Mpz, usize) {
        let divisor_limbs = divisor.limbs();
        let highest_limb = divisor_limbs[divisor_limbs.len() - 1];

        if highest_limb >= 0x8000000000000000 {
            return (divisor.clone(), 0);
        }

        // 计算需要左移的位数
        let shift = highest_limb.leading_zeros() as usize;
        let normalized = divisor.shl(shift);

        (normalized, shift)
    }

    /// 估算商位
    fn estimate_quotient_digit(
        &self,
        remainder_limbs: &[u64],
        divisor_limbs: &[u64],
        j: usize,
    ) -> u64 {
        let n = divisor_limbs.len();

        if j + n >= remainder_limbs.len() {
            return 0;
        }

        let r_high = remainder_limbs[j + n] as u128;
        let r_mid = remainder_limbs[j + n - 1] as u128;
        let d_high = divisor_limbs[n - 1] as u128;

        // 使用Knuth的估算公式
        let q_hat = if r_high >= d_high {
            0xFFFFFFFFFFFFFFFF
        } else {
            let temp = ((r_high << 64) | r_mid) / d_high;
            temp.min(0xFFFFFFFFFFFFFFFF)
        };

        q_hat as u64
    }

    /// 修正商位（如果需要）
    fn correct_quotient_digit(&self, remainder_limbs: &mut [u64], divisor_limbs: &[u64], j: usize) {
        let n = divisor_limbs.len();

        // 添加除数
        let mut carry = 0u128;
        for i in 0..n {
            let sum = remainder_limbs[j + i] as u128 + divisor_limbs[i] as u128 + carry;
            remainder_limbs[j + i] = sum as u64;
            carry = sum >> 64;
        }

        if j + n < remainder_limbs.len() {
            remainder_limbs[j + n] = remainder_limbs[j + n].wrapping_add(carry as u64);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_addition() {
        let a = Mpz::from_i64(123);
        let b = Mpz::from_i64(456);
        let result = a.add(&b);
        assert_eq!(result.to_i64(), Some(579));

        // 测试负数
        let c = Mpz::from_i64(-123);
        let d = Mpz::from_i64(456);
        let result2 = c.add(&d);
        assert_eq!(result2.to_i64(), Some(333));
    }

    #[test]
    fn test_subtraction() {
        let a = Mpz::from_i64(456);
        let b = Mpz::from_i64(123);
        let result = a.sub(&b);
        assert_eq!(result.to_i64(), Some(333));

        // 测试结果为负数
        let result2 = b.sub(&a);
        assert_eq!(result2.to_i64(), Some(-333));
    }

    #[test]
    fn test_multiplication() {
        let a = Mpz::from_i64(123);
        let b = Mpz::from_i64(456);
        let result = a.mul(&b);
        assert_eq!(result.to_i64(), Some(56088));

        // 测试负数
        let c = Mpz::from_i64(-123);
        let result2 = c.mul(&b);
        assert_eq!(result2.to_i64(), Some(-56088));
    }

    #[test]
    fn test_division() {
        let a = Mpz::from_i64(456);
        let b = Mpz::from_i64(123);
        let result = a.div(&b).unwrap();
        assert_eq!(result.to_i64(), Some(3));

        // 测试除零
        let zero = Mpz::new();
        assert!(a.div(&zero).is_err());
    }

    #[test]
    fn test_remainder() {
        let a = Mpz::from_i64(456);
        let b = Mpz::from_i64(123);
        let result = a.rem(&b).unwrap();
        assert_eq!(result.to_i64(), Some(87));

        // 测试mod_别名
        let result2 = a.mod_(&b).unwrap();
        assert_eq!(result, result2);
    }

    #[test]
    fn test_div_rem() {
        let a = Mpz::from_i64(456);
        let b = Mpz::from_i64(123);
        let (quotient, remainder) = a.div_rem(&b).unwrap();
        assert_eq!(quotient.to_i64(), Some(3));
        assert_eq!(remainder.to_i64(), Some(87));
    }

    #[test]
    fn test_large_division() {
        // 测试大数除法
        let a = Mpz::from_str("123456789012345678901234567890", 10).unwrap();
        let b = Mpz::from_str("987654321", 10).unwrap();
        let (quotient, remainder) = a.div_rem(&b).unwrap();

        // 验证: a = b * quotient + remainder
        let product = b.mul(&quotient);
        let sum = product.add(&remainder);
        assert_eq!(sum, a);

        // 验证余数小于除数
        assert!(remainder.cmp_abs(&b) == core::cmp::Ordering::Less);
    }

    #[test]
    fn test_division_edge_cases() {
        // 测试边界情况
        let zero = Mpz::new();
        let one = Mpz::from_i64(1);
        let large = Mpz::from_str("12345678901234567890", 10).unwrap();

        // 0 / 1 = 0
        let (q, r) = zero.div_rem(&one).unwrap();
        assert_eq!(q, zero);
        assert_eq!(r, zero);

        // 1 / 1 = 1
        let (q, r) = one.div_rem(&one).unwrap();
        assert_eq!(q, one);
        assert_eq!(r, zero);

        // large / 1 = large
        let (q, r) = large.div_rem(&one).unwrap();
        assert_eq!(q, large);
        assert_eq!(r, zero);

        // large / large = 1
        let (q, r) = large.div_rem(&large).unwrap();
        assert_eq!(q, one);
        assert_eq!(r, zero);
    }

    #[test]
    fn test_division_performance() {
        // 测试除法性能
        let a = Mpz::from_str("2", 10).unwrap().pow_u32(1000);
        let b = Mpz::from_str("3", 10).unwrap().pow_u32(500);

        let start = std::time::Instant::now();
        let (quotient, remainder) = a.div_rem(&b).unwrap();
        let duration = start.elapsed();

        // 验证结果正确性
        let product = b.mul(&quotient);
        let sum = product.add(&remainder);
        assert_eq!(sum, a);

        println!("大数除法耗时: {:?}", duration);
        assert!(duration.as_millis() < 1000); // 应该在1秒内完成
    }
}
