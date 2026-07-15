//! 位运算和位操作实现

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

use super::core::Mpz;

impl Mpz {
    /// 按位与运算
    pub fn bitwise_and(&self, other: &Mpz) -> Mpz {
        if self.is_zero() || other.is_zero() {
            return Mpz::new();
        }

        let self_limbs = self.limbs();
        let other_limbs = other.limbs();
        let min_len = self_limbs.len().min(other_limbs.len());

        let mut result = Vec::with_capacity(min_len);

        for i in 0..min_len {
            result.push(self_limbs[i] & other_limbs[i]);
        }

        // 对于负数的处理需要更复杂的逻辑，这里先实现正数版本
        Mpz::from_limbs(result, false)
    }

    /// 按位或运算
    pub fn bitwise_or(&self, other: &Mpz) -> Mpz {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }

        let self_limbs = self.limbs();
        let other_limbs = other.limbs();
        let max_len = self_limbs.len().max(other_limbs.len());

        let mut result = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let self_limb = self_limbs.get(i).copied().unwrap_or(0);
            let other_limb = other_limbs.get(i).copied().unwrap_or(0);
            result.push(self_limb | other_limb);
        }

        Mpz::from_limbs(result, false)
    }

    /// 按位异或运算
    pub fn bitwise_xor(&self, other: &Mpz) -> Mpz {
        let self_limbs = self.limbs();
        let other_limbs = other.limbs();
        let max_len = self_limbs.len().max(other_limbs.len());

        let mut result = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let self_limb = self_limbs.get(i).copied().unwrap_or(0);
            let other_limb = other_limbs.get(i).copied().unwrap_or(0);
            result.push(self_limb ^ other_limb);
        }

        Mpz::from_limbs(result, false)
    }

    /// 按位取反运算
    pub fn bitwise_not(&self) -> Mpz {
        if self.is_zero() {
            // ~0 = -1，但这需要确定位长度
            // 这里返回一个表示 -1 的值
            return Mpz::from_i64(-1);
        }

        let self_limbs = self.limbs();
        let mut result = Vec::with_capacity(self_limbs.len());

        for &limb in self_limbs {
            result.push(!limb);
        }

        // 按位取反的符号处理比较复杂，这里先简化实现
        Mpz::from_limbs(result, !self.is_negative())
    }

    /// 左移运算
    pub fn shl(&self, bits: usize) -> Mpz {
        if bits == 0 || self.is_zero() {
            return self.clone();
        }

        let limb_shift = bits / 64;
        let bit_shift = bits % 64;

        let old_limbs = self.limbs();
        let mut new_limbs = vec![0u64; old_limbs.len() + limb_shift + 1];

        if bit_shift == 0 {
            // 只有limb级别的移位
            new_limbs[limb_shift..limb_shift + old_limbs.len()].copy_from_slice(old_limbs);
        } else {
            // 需要位级别的移位
            let mut carry = 0u64;
            for i in 0..old_limbs.len() {
                let limb = old_limbs[i];
                new_limbs[i + limb_shift] = (limb << bit_shift) | carry;
                carry = limb >> (64 - bit_shift);
            }
            if carry != 0 {
                new_limbs[old_limbs.len() + limb_shift] = carry;
            }
        }

        Mpz::from_limbs(new_limbs, self.is_negative())
    }

    /// 右移运算
    pub fn shr(&self, bits: usize) -> Mpz {
        if bits == 0 || self.is_zero() {
            return self.clone();
        }

        let limb_shift = bits / 64;
        let bit_shift = bits % 64;

        let old_limbs = self.limbs();

        // 如果移位超过了数的位长，结果为0（正数）或-1（负数）
        if limb_shift >= old_limbs.len() {
            return if self.is_negative() {
                Mpz::from_i64(-1)
            } else {
                Mpz::new()
            };
        }

        let new_len = old_limbs.len() - limb_shift;
        let mut new_limbs = Vec::with_capacity(new_len);

        if bit_shift == 0 {
            // 只有limb级别的移位
            new_limbs.extend_from_slice(&old_limbs[limb_shift..]);
        } else {
            // 需要位级别的移位
            for i in limb_shift..old_limbs.len() {
                let mut limb = old_limbs[i] >> bit_shift;
                if i + 1 < old_limbs.len() {
                    limb |= old_limbs[i + 1] << (64 - bit_shift);
                }
                new_limbs.push(limb);
            }
        }

        if new_limbs.is_empty() {
            new_limbs.push(0);
        }

        Mpz::from_limbs(new_limbs, self.is_negative())
    }

    /// 测试指定位是否为1
    pub fn test_bit(&self, bit: usize) -> bool {
        let limb_index = bit / 64;
        let bit_index = bit % 64;

        if limb_index >= self.limbs().len() {
            // 超出范围的位，对于正数是0，对于负数需要考虑二进制补码
            self.is_negative()
        } else {
            (self.limbs()[limb_index] & (1u64 << bit_index)) != 0
        }
    }

    /// 设置指定位为1
    pub fn set_bit(&mut self, bit: usize) {
        let limb_index = bit / 64;
        let bit_index = bit % 64;

        // 确保有足够的limbs
        while self.limbs().len() <= limb_index {
            self.limbs_mut().push(0);
        }

        self.limbs_mut()[limb_index] |= 1u64 << bit_index;

        // 重新标准化
        self.normalize();
    }

    /// 清除指定位（设为0）
    pub fn clear_bit(&mut self, bit: usize) {
        let limb_index = bit / 64;
        let bit_index = bit % 64;

        if limb_index < self.limbs().len() {
            self.limbs_mut()[limb_index] &= !(1u64 << bit_index);
            self.normalize();
        }
    }

    /// 翻转指定位
    pub fn flip_bit(&mut self, bit: usize) {
        let limb_index = bit / 64;
        let bit_index = bit % 64;

        // 确保有足够的limbs
        while self.limbs().len() <= limb_index {
            self.limbs_mut().push(0);
        }

        self.limbs_mut()[limb_index] ^= 1u64 << bit_index;

        // 重新标准化
        self.normalize();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitwise_and() {
        let a = Mpz::from_u64(0b1010);
        let b = Mpz::from_u64(0b1100);
        let result = a.bitwise_and(&b);
        assert_eq!(result.to_u64(), Some(0b1000));
    }

    #[test]
    fn test_bitwise_or() {
        let a = Mpz::from_u64(0b1010);
        let b = Mpz::from_u64(0b1100);
        let result = a.bitwise_or(&b);
        assert_eq!(result.to_u64(), Some(0b1110));
    }

    #[test]
    fn test_bitwise_xor() {
        let a = Mpz::from_u64(0b1010);
        let b = Mpz::from_u64(0b1100);
        let result = a.bitwise_xor(&b);
        assert_eq!(result.to_u64(), Some(0b0110));
    }

    #[test]
    fn test_shifts() {
        let a = Mpz::from_u64(0b1010);

        let left = a.shl(2);
        assert_eq!(left.to_u64(), Some(0b101000));

        let right = left.shr(2);
        assert_eq!(right.to_u64(), Some(0b1010));
    }

    #[test]
    fn test_bit_operations() {
        let mut a = Mpz::from_u64(0b1010);

        assert!(a.test_bit(1));
        assert!(!a.test_bit(0));
        assert!(a.test_bit(3));
        assert!(!a.test_bit(2));

        a.set_bit(0);
        assert_eq!(a.to_u64(), Some(0b1011));

        a.clear_bit(3);
        assert_eq!(a.to_u64(), Some(0b0011));

        a.flip_bit(2);
        assert_eq!(a.to_u64(), Some(0b0111));
    }

    #[test]
    fn test_bit_length() {
        assert_eq!(Mpz::new().bit_length(), 0);
        assert_eq!(Mpz::from_u64(1).bit_length(), 1);
        assert_eq!(Mpz::from_u64(2).bit_length(), 2);
        assert_eq!(Mpz::from_u64(3).bit_length(), 2);
        assert_eq!(Mpz::from_u64(255).bit_length(), 8);
        assert_eq!(Mpz::from_u64(256).bit_length(), 9);
    }
}
