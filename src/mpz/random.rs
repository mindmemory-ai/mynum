//! 随机数生成实现

use super::core::Mpz;
use crate::error::{Error, Result};

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

impl Mpz {
    /// 生成指定位数的随机大整数
    pub fn random_bits(bits: usize) -> Result<Mpz> {
        if bits == 0 {
            return Ok(Mpz::new());
        }

        let mut rng = rand::thread_rng();
        Self::random_with_rng(&mut rng, bits)
    }

    /// 生成指定范围内的随机大整数 [min, max)
    pub fn random_range(min: &Mpz, max: &Mpz) -> Result<Mpz> {
        if min >= max {
            return Err(Error::InvalidInput("min must be less than max".into()));
        }

        let mut rng = rand::thread_rng();
        Self::random_range_with_rng(&mut rng, min, max)
    }

    /// 生成指定范围内的随机大整数（包含上下界） [min, max]
    pub fn random_range_inclusive(min: &Mpz, max: &Mpz) -> Result<Mpz> {
        if min > max {
            return Err(Error::InvalidInput(
                "min must be less than or equal to max".into(),
            ));
        }

        // 对于包含范围，我们生成 [min, max+1) 然后检查是否等于 max+1
        let max_plus_one = max.add(&Mpz::from_i64(1));
        let result = Self::random_range(min, &max_plus_one)?;

        if result == max_plus_one {
            Ok(max.clone())
        } else {
            Ok(result)
        }
    }

    /// 生成小于指定值的随机大整数 [0, max)
    pub fn random_below(max: &Mpz) -> Result<Mpz> {
        if max.is_zero() || max.is_negative() {
            return Err(Error::InvalidInput("max must be positive".into()));
        }

        Self::random_range(&Mpz::new(), max)
    }

    /// 生成随机奇数（指定位数）
    pub fn random_odd_bits(bits: usize) -> Result<Mpz> {
        if bits == 0 {
            return Ok(Mpz::new());
        }

        let mut rng = rand::thread_rng();
        let mut result = Self::random_with_rng(&mut rng, bits)?;

        // 确保是奇数（设置最低位为1）
        if !result.is_zero() {
            result.set_bit(0);
        }

        Ok(result)
    }

    /// 使用指定的随机数生成器生成随机数
    pub fn random_with_rng<R: rand::Rng>(rng: &mut R, bits: usize) -> Result<Mpz> {
        if bits == 0 {
            return Ok(Mpz::new());
        }

        let limbs_needed = bits.div_ceil(64);
        let mut limbs = Vec::with_capacity(limbs_needed);

        // 生成完整的limbs
        for _ in 0..limbs_needed - 1 {
            limbs.push(rng.gen::<u64>());
        }

        // 生成最后一个limb，确保不超过指定位数
        let remaining_bits = bits % 64;
        let last_limb = if remaining_bits == 0 {
            rng.gen::<u64>()
        } else {
            let mask = (1u64 << remaining_bits) - 1;
            rng.gen::<u64>() & mask
        };

        limbs.push(last_limb);

        // 确保最高位不为0（除非是0）
        if bits > 0 && last_limb == 0 && limbs.len() > 1 {
            let highest_bit = 1u64 << (remaining_bits - 1);
            let last_index = limbs.len() - 1;
            limbs[last_index] = highest_bit;
        }

        Ok(Mpz::from_limbs(limbs, false))
    }

    /// 使用指定的随机数生成器生成范围内的随机数
    pub fn random_range_with_rng<R: rand::Rng>(rng: &mut R, min: &Mpz, max: &Mpz) -> Result<Mpz> {
        if min >= max {
            return Err(Error::InvalidInput("min must be less than max".into()));
        }

        let range = max.sub(min);
        let range_bits = range.bit_length();

        loop {
            let candidate = Self::random_with_rng(rng, range_bits)?;

            if candidate < range {
                return Ok(min.add(&candidate));
            }
        }
    }

    /// 生成密码学安全的随机数（指定位数）
    pub fn cryptographically_secure_random(bits: usize) -> Result<Mpz> {
        if bits == 0 {
            return Ok(Mpz::new());
        }

        // 使用系统随机数生成器
        let mut rng = rand::thread_rng();
        Self::random_with_rng(&mut rng, bits)
    }

    /// 生成密码学安全的随机素数
    pub fn cryptographically_secure_prime(bits: usize) -> Result<Mpz> {
        if bits < 2 {
            return Err(Error::InvalidInput("Prime must be at least 2 bits".into()));
        }

        let mut rng = rand::thread_rng();

        loop {
            // 生成随机奇数
            let mut candidate = Self::random_with_rng(&mut rng, bits)?;
            if candidate.is_zero() {
                candidate = Mpz::from_i64(1);
            }

            // 确保是奇数
            if !candidate.test_bit(0) {
                candidate.set_bit(0);
            }

            // 确保最高位为1
            if !candidate.test_bit(bits - 1) {
                candidate.set_bit(bits - 1);
            }

            // 简单素数测试
            if candidate.is_probably_prime(25) {
                return Ok(candidate);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random_bits() {
        let result = Mpz::random_bits(64).unwrap();
        assert!(result.bit_length() <= 64);
        assert!(!result.is_negative());
    }

    #[test]
    fn test_random_range() {
        let min = Mpz::from_i64(100);
        let max = Mpz::from_i64(200);

        for _ in 0..100 {
            let result = Mpz::random_range(&min, &max).unwrap();
            assert!(result >= min);
            assert!(result < max);
        }
    }

    #[test]
    fn test_random_below() {
        let max = Mpz::from_i64(100);

        for _ in 0..100 {
            let result = Mpz::random_below(&max).unwrap();
            assert!(result >= Mpz::new());
            assert!(result < max);
        }
    }

    #[test]
    fn test_random_odd_bits() {
        let result = Mpz::random_odd_bits(64).unwrap();
        assert!(result.test_bit(0)); // 最低位应该是1
        assert!(result.bit_length() <= 64);
    }

    #[test]
    fn test_cryptographically_secure_random() {
        let result = Mpz::cryptographically_secure_random(128).unwrap();
        assert!(result.bit_length() <= 128);
        assert!(!result.is_negative());
    }
}
