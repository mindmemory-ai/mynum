//! Mpf 算术运算实现
//!
//! Shared Newton-iteration constants (max iterations, convergence threshold)
//! are defined in [`crate::algorithm`] — used by both this module and
//! [`crate::mpz::arithmetic`].

use super::core::Mpf;
use crate::config::RoundingMode;
use crate::error::{Error, Result};
use crate::mpz::Mpz;

impl Mpf {
    /// 加法运算
    pub fn add(&self, other: &Mpf) -> Mpf {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }

        // 对齐指数
        let (aligned_self, aligned_other) = self.align_exponents(other);

        // 执行加法
        let result_mantissa = if aligned_self.is_negative() == aligned_other.is_negative() {
            // 同号相加
            aligned_self.mantissa().add(aligned_other.mantissa())
        } else {
            // 异号相减
            if aligned_self.mantissa() >= aligned_other.mantissa() {
                aligned_self.mantissa().sub(aligned_other.mantissa())
            } else {
                aligned_other.mantissa().sub(aligned_self.mantissa())
            }
        };

        let result_negative = if aligned_self.is_negative() == aligned_other.is_negative() {
            aligned_self.is_negative()
        } else {
            // 异号相减时，符号取决于绝对值较大的数
            if aligned_self.mantissa() >= aligned_other.mantissa() {
                aligned_self.is_negative()
            } else {
                aligned_other.is_negative()
            }
        };

        let result = Mpf::from_parts_with_sign(
            result_mantissa,
            aligned_self.exponent(),
            self.precision(),
            result_negative,
        );
        // 不调用normalize，保持精度
        result
    }

    /// 减法运算
    pub fn sub(&self, other: &Mpf) -> Mpf {
        let neg_other = other.neg();
        self.add(&neg_other)
    }

    /// 乘法运算
    pub fn mul(&self, other: &Mpf) -> Mpf {
        if self.is_zero() || other.is_zero() {
            return Mpf::new();
        }

        let result_mantissa = self.mantissa().mul(other.mantissa());
        let result_exponent = self.exponent() + other.exponent();
        let result_negative = self.is_negative() != other.is_negative();

        let mut result = Mpf::from_parts_with_sign(
            result_mantissa,
            result_exponent,
            self.precision(),
            result_negative,
        );
        result.normalize();
        result
    }

    /// 除法运算
    pub fn div(&self, other: &Mpf) -> Result<Mpf> {
        if other.is_zero() {
            return Err(Error::DivisionByZero);
        }

        if self.is_zero() {
            return Ok(Mpf::new());
        }

        // 对齐指数，确保分子和分母有相同的指数
        let (aligned_self, aligned_other) = self.align_exponents(other);

        let mut dividend = aligned_self.mantissa().clone();
        let divisor = aligned_other.mantissa();

        // Shift dividend left by precision bits to get sufficient fractional
        // precision in the integer division result. Without this shift, integer
        // division truncates the fractional part (e.g., 1/3 = 0 instead of 0.333...).
        let precision = self.precision();
        dividend = dividend.shl(precision);

        // 执行除法
        let (quotient, _remainder) = dividend.div_rem(divisor)?;

        // Adjust exponent: we shifted dividend left by `precision` bits,
        // so we must subtract `precision` from the result exponent.
        let result_exponent = aligned_self.exponent() - aligned_other.exponent() - precision as i64;
        let result_negative = aligned_self.is_negative() != aligned_other.is_negative();

        let mut result =
            Mpf::from_parts_with_sign(quotient, result_exponent, self.precision(), result_negative);
        result.normalize();
        Ok(result)
    }

    /// 取负
    pub fn neg(&self) -> Mpf {
        let mut result = self.clone();
        result.set_negative(!result.is_negative());
        result
    }

    /// 幂运算
    pub fn pow(&self, exp: u32) -> Result<Mpf> {
        if exp == 0 {
            return Ok(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));
        }

        if self.is_zero() {
            return Ok(Mpf::new());
        }

        let mut result = Mpf::from_mpz(Mpz::from_i64(1), self.precision());
        let mut base = self.clone();
        let mut exponent = exp;

        while exponent > 0 {
            if exponent % 2 == 1 {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            exponent /= 2;
        }

        Ok(result)
    }

    /// Raise self to an arbitrary-precision integer exponent (binary exponentiation).
    ///
    /// Uses repeated squaring. The exponent must be non-negative; negative
    /// exponents are rejected with `DomainError`.
    pub fn pow_big(&self, exp: &Mpz) -> Result<Mpf> {
        if exp.is_zero() {
            return Ok(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));
        }
        if self.is_zero() {
            return Ok(Mpf::new());
        }
        if exp.is_negative() {
            return Err(Error::DomainError(
                "negative exponent not supported for Mpf::pow_big".into(),
            ));
        }

        let mut result = Mpf::from_mpz(Mpz::from_i64(1), self.precision());
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
        Ok(result)
    }

    /// 舍入到指定精度
    pub fn round(&self, precision: usize, mode: RoundingMode) -> Mpf {
        if self.precision() <= precision {
            return self.clone();
        }

        let mut result = self.clone();
        let shift = self.precision() - precision;

        // 根据舍入模式进行舍入
        let remainder = self
            .mantissa()
            .rem(&Mpz::from_i64(1).shl(shift))
            .expect("Failed to compute remainder for rounding");

        let should_round_up = match mode {
            RoundingMode::TowardNearest => {
                let half_way = Mpz::from_i64(1).shl(shift - 1);
                remainder > half_way || (remainder == half_way && self.mantissa().test_bit(shift))
            }
            RoundingMode::TowardZero => false,
            RoundingMode::TowardPositive => !self.is_negative(),
            RoundingMode::TowardNegative => self.is_negative(),
        };

        if should_round_up {
            *result.mantissa_mut() = result.mantissa().add(&Mpz::from_i64(1).shl(shift));
        }

        *result.mantissa_mut() = result.mantissa().shr(shift);
        *result.exponent_mut() += shift as i64;
        result.set_precision(precision);
        result.normalize();
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_addition() {
        let a = Mpf::from_mpz(Mpz::from_i64(123), 64);
        let b = Mpf::from_mpz(Mpz::from_i64(456), 64);
        let result = a.add(&b);
        assert_eq!(result.to_i64(), Some(579));
    }

    #[test]
    fn test_subtraction() {
        let a = Mpf::from_mpz(Mpz::from_i64(456), 64);
        let b = Mpf::from_mpz(Mpz::from_i64(123), 64);
        let result = a.sub(&b);
        assert_eq!(result.to_i64(), Some(333));
    }

    #[test]
    fn test_multiplication() {
        let a = Mpf::from_mpz(Mpz::from_i64(12), 64);
        let b = Mpf::from_mpz(Mpz::from_i64(34), 64);
        let result = a.mul(&b);
        assert_eq!(result.to_i64(), Some(408));
    }

    #[test]
    fn test_division() {
        let a = Mpf::from_mpz(Mpz::from_i64(100), 64);
        let b = Mpf::from_mpz(Mpz::from_i64(25), 64);
        let result = a.div(&b).unwrap();
        assert_eq!(result.to_i64(), Some(4));
    }
}
