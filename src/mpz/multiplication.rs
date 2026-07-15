//! 乘法后端实现
//!
//! 包含不同的乘法算法实现

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

use crate::mpz::parallel::{get_global_parallel, ParallelExecutor};
use crate::mpz::{Mpz, MpzMultiplicationConfig, MultiplicationBackend};

impl Mpz {
    /// 使用指定后端进行乘法运算
    pub fn mul_with_backend(&self, other: &Mpz, backend: MultiplicationBackend) -> Mpz {
        match backend {
            MultiplicationBackend::Schoolbook => self.mul_schoolbook(other),
            MultiplicationBackend::Karatsuba => self.mul_karatsuba(other),
            MultiplicationBackend::ToomCook3 => self.mul_toom_cook(other),
            MultiplicationBackend::ToomCook4 => self.mul_toom_cook(other),
            MultiplicationBackend::FFT => self.mul_fft(other),
            MultiplicationBackend::NTT => self.mul_fft(other), // 暂时使用FFT
            MultiplicationBackend::Adaptive => self.mul_adaptive(other),
            MultiplicationBackend::ParallelKaratsuba => self.mul_parallel_karatsuba(other),
            MultiplicationBackend::ParallelFFT => self.mul_parallel_fft(other),
            MultiplicationBackend::Parallel => self.mul_parallel(other),
            MultiplicationBackend::Custom => self.mul_schoolbook(other), // 默认回退
        }
    }

    /// 自适应乘法（根据操作数大小自动选择最优算法）
    pub fn mul_adaptive(&self, other: &Mpz) -> Mpz {
        // 检查是否应该使用并行算法
        if MpzMultiplicationConfig::is_parallel_enabled()
            && self.bit_length() >= MpzMultiplicationConfig::get_parallel_threshold()
            && other.bit_length() >= MpzMultiplicationConfig::get_parallel_threshold()
        {
            return self.mul_parallel(other);
        }

        // 使用建议的后端
        let suggested_backend =
            MpzMultiplicationConfig::suggest_backend(self.bit_length(), other.bit_length());
        self.mul_with_backend(other, suggested_backend)
    }

    /// 并行乘法（智能选择并行Karatsuba或并行FFT）
    /// 使用 ParallelExecutor 统一后端进行子任务并行化。
    pub fn mul_parallel(&self, other: &Mpz) -> Mpz {
        let max_bits = self.bit_length().max(other.bit_length());

        if max_bits <= 4096 {
            // 使用 ParallelExecutor 并行计算 Karatsuba 三个子乘积
            self.mul_parallel_karatsuba_v2(other)
        } else {
            // 大数字使用并行FFT
            self.mul_parallel_fft(other)
        }
    }

    /// 使用 ParallelExecutor::map_reduce 并行计算 Karatsuba 子乘积
    fn mul_parallel_karatsuba_v2(&self, other: &Mpz) -> Mpz {
        let max_bits = self.bit_length().max(other.bit_length());

        // 数字太小，不需要并行
        if max_bits < 512 {
            return self.mul(other);
        }

        let half = max_bits.div_ceil(2);
        if half < 256 {
            return self.mul(other);
        }

        // 分割数字为高位和低位
        let mask = Mpz::from_i64(1).shl(half).sub(&Mpz::from_i64(1));
        let a_low = self.bitwise_and(&mask);
        let a_high = self.shr(half);
        let b_low = other.bitwise_and(&mask);
        let b_high = other.shr(half);

        // 三个子任务：(a_low, b_low), (a_high, b_high), (a_low+a_high, b_low+b_high)
        let a_sum = a_low.add(&a_high);
        let b_sum = b_low.add(&b_high);

        let tasks = [
            (&a_low, &b_low),   // p0 = a_low * b_low
            (&a_high, &b_high), // p2 = a_high * b_high
            (&a_sum, &b_sum),   // p1_full = (a_low+a_high) * (b_low+b_high)
        ];

        let executor = ParallelExecutor::new();
        let results: Vec<Mpz> = executor.map_reduce(&tasks, |(a, b)| a.mul(b));

        let p0 = &results[0];
        let p2 = &results[1];
        let p1_full = &results[2];
        let p1 = p1_full.sub(p0).sub(p2);

        // 使用位运算优化移位操作
        let shift_half = Mpz::from_i64(1).shl(half);
        let term1 = p1.mul(&shift_half);
        let term2 = p2.mul(&shift_half).mul(&shift_half);
        p0.add(&term1).add(&term2)
    }

    /// 并行Karatsuba乘法
    pub fn mul_parallel_karatsuba(&self, other: &Mpz) -> Mpz {
        get_global_parallel()
            .karatsuba_parallel_optimized(self, other)
            .unwrap_or_else(|_| self.mul_karatsuba(other))
    }

    /// 并行FFT乘法
    pub fn mul_parallel_fft(&self, other: &Mpz) -> Mpz {
        get_global_parallel()
            .parallel_fft_multiply(self, other)
            .unwrap_or_else(|_| self.mul_fft(other))
    }

    /// 使用基础乘法算法（学校乘法）
    pub fn mul_schoolbook(&self, other: &Mpz) -> Mpz {
        let a_limbs = self.limbs();
        let b_limbs = other.limbs();

        // 单limb优化
        if a_limbs.len() == 1 && b_limbs.len() == 1 {
            let result = (a_limbs[0] as u128) * (b_limbs[0] as u128);
            let low = result as u64;
            let high = (result >> 64) as u64;

            if high == 0 {
                return Mpz::from_limbs(vec![low], false);
            } else {
                return Mpz::from_limbs(vec![low, high], false);
            }
        }

        let mut result = vec![0u64; a_limbs.len() + b_limbs.len()];

        for i in 0..a_limbs.len() {
            let mut carry = 0u128;

            for j in 0..b_limbs.len() {
                let product =
                    (a_limbs[i] as u128) * (b_limbs[j] as u128) + (result[i + j] as u128) + carry;

                result[i + j] = product as u64;
                carry = product >> 64;
            }

            if carry != 0 {
                result[i + b_limbs.len()] = carry as u64;
            }
        }

        Mpz::from_limbs(result, false)
    }

    /// 使用Karatsuba算法
    pub fn mul_karatsuba(&self, other: &Mpz) -> Mpz {
        let a_limbs = self.limbs();
        let b_limbs = other.limbs();

        // 对于小数使用基础算法
        if a_limbs.len() <= 32 || b_limbs.len() <= 32 {
            return self.mul_schoolbook(other);
        }

        // 确保两个数长度相等，用零填充
        let max_len = a_limbs.len().max(b_limbs.len());
        let mut a_padded = a_limbs.to_vec();
        let mut b_padded = b_limbs.to_vec();

        a_padded.resize(max_len, 0);
        b_padded.resize(max_len, 0);

        let result = self.karatsuba_multiply(&a_padded, &b_padded);
        Mpz::from_limbs(result, false)
    }

    /// Karatsuba乘法的递归实现
    fn karatsuba_multiply(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        let n = a.len();

        // 基础情况：小数组使用基础乘法
        if n <= 32 {
            return self.schoolbook_multiply_limbs(a, b);
        }

        // 分割数组
        let m = n / 2;
        let (a_low, a_high) = a.split_at(m);
        let (b_low, b_high) = b.split_at(m);

        // 递归计算三个乘积
        let z0 = self.karatsuba_multiply(a_low, b_low); // a_low * b_low
        let z2 = self.karatsuba_multiply(a_high, b_high); // a_high * b_high

        // 计算 (a_low + a_high) * (b_low + b_high)
        let a_sum = self.add_limb_arrays(a_low, a_high);
        let b_sum = self.add_limb_arrays(b_low, b_high);
        let z1_full = self.karatsuba_multiply(&a_sum, &b_sum);

        // z1 = z1_full - z0 - z2
        let z1 = self.subtract_and_shift(&z1_full, &z0, &z2);

        // 合并结果: z0 + z1 * base^m + z2 * base^(2m)
        self.combine_karatsuba_results(&z0, &z1, &z2, m)
    }

    /// 基础乘法（用于limb数组）
    fn schoolbook_multiply_limbs(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        let mut result = vec![0u64; a.len() + b.len()];

        for i in 0..a.len() {
            let mut carry = 0u128;

            for j in 0..b.len() {
                let product = (a[i] as u128) * (b[j] as u128) + (result[i + j] as u128) + carry;

                result[i + j] = product as u64;
                carry = product >> 64;
            }

            if carry != 0 && i + b.len() < result.len() {
                result[i + b.len()] = carry as u64;
            }
        }

        result
    }

    /// 添加两个limb数组
    fn add_limb_arrays(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        let max_len = a.len().max(b.len());
        let mut result = Vec::with_capacity(max_len + 1);
        let mut carry = 0u64;

        for i in 0..max_len {
            let a_val = a.get(i).copied().unwrap_or(0);
            let b_val = b.get(i).copied().unwrap_or(0);

            let (sum1, overflow1) = a_val.overflowing_add(b_val);
            let (sum2, overflow2) = sum1.overflowing_add(carry);

            result.push(sum2);
            carry = if overflow1 || overflow2 { 1 } else { 0 };
        }

        if carry != 0 {
            result.push(carry);
        }

        result
    }

    /// 计算 z1_full - z0 - z2
    fn subtract_and_shift(&self, z1_full: &[u64], z0: &[u64], z2: &[u64]) -> Vec<u64> {
        // 简化实现：先计算z1_full，然后减去z0和z2
        let mut result = z1_full.to_vec();

        // 减去z0
        let mut borrow = 0u64;
        for i in 0..z0.len().min(result.len()) {
            let (diff1, underflow1) = result[i].overflowing_sub(z0[i]);
            let (diff2, underflow2) = diff1.overflowing_sub(borrow);
            result[i] = diff2;
            borrow = if underflow1 || underflow2 { 1 } else { 0 };
        }

        // 减去z2
        borrow = 0;
        for i in 0..z2.len().min(result.len()) {
            let (diff1, underflow1) = result[i].overflowing_sub(z2[i]);
            let (diff2, underflow2) = diff1.overflowing_sub(borrow);
            result[i] = diff2;
            borrow = if underflow1 || underflow2 { 1 } else { 0 };
        }

        result
    }

    /// 合并Karatsuba结果
    fn combine_karatsuba_results(&self, z0: &[u64], z1: &[u64], z2: &[u64], m: usize) -> Vec<u64> {
        let result_len = z0.len().max(z1.len() + m).max(z2.len() + 2 * m);
        let mut result = vec![0u64; result_len];

        // 添加 z0
        for (i, &val) in z0.iter().enumerate() {
            if i < result.len() {
                let (sum, overflow) = result[i].overflowing_add(val);
                result[i] = sum;
                if overflow && i + 1 < result.len() {
                    result[i + 1] = result[i + 1].wrapping_add(1);
                }
            }
        }

        // 添加 z1 * base^m
        for (i, &val) in z1.iter().enumerate() {
            let pos = i + m;
            if pos < result.len() {
                let (sum, overflow) = result[pos].overflowing_add(val);
                result[pos] = sum;
                if overflow && pos + 1 < result.len() {
                    result[pos + 1] = result[pos + 1].wrapping_add(1);
                }
            }
        }

        // 添加 z2 * base^(2m)
        for (i, &val) in z2.iter().enumerate() {
            let pos = i + 2 * m;
            if pos < result.len() {
                let (sum, overflow) = result[pos].overflowing_add(val);
                result[pos] = sum;
                if overflow && pos + 1 < result.len() {
                    result[pos + 1] = result[pos + 1].wrapping_add(1);
                }
            }
        }

        result
    }

    /// 使用Toom-Cook算法（简化版本）
    pub fn mul_toom_cook(&self, other: &Mpz) -> Mpz {
        // 对于中等大小的数，暂时使用Karatsuba
        // TODO: 实现真正的Toom-Cook算法
        if self.limb_count() < 1024 && other.limb_count() < 1024 {
            self.mul_karatsuba(other)
        } else {
            self.mul_schoolbook(other)
        }
    }

    /// 使用FFT算法（占位符）
    pub fn mul_fft(&self, other: &Mpz) -> Mpz {
        // TODO: 实现FFT乘法算法
        // 对于非常大的数，这将提供最佳性能
        self.mul_karatsuba(other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_schoolbook_multiplication() {
        let a = Mpz::from_u64(12345);
        let b = Mpz::from_u64(67890);
        let result = a.mul_schoolbook(&b);
        let expected = a.mul(&b);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_karatsuba_multiplication() {
        let a = Mpz::from_str("123456789012345678901234567890", 10).unwrap();
        let b = Mpz::from_str("987654321098765432109876543210", 10).unwrap();
        let result = a.mul_karatsuba(&b);
        let expected = a.mul(&b);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_mul_with_backend() {
        let a = Mpz::from_u64(12345);
        let b = Mpz::from_u64(67890);

        let schoolbook = a.mul_with_backend(&b, MultiplicationBackend::Schoolbook);
        let karatsuba = a.mul_with_backend(&b, MultiplicationBackend::Karatsuba);
        let adaptive = a.mul_with_backend(&b, MultiplicationBackend::Adaptive);

        // 所有算法应该得到相同的结果
        assert_eq!(schoolbook, karatsuba);
        assert_eq!(karatsuba, adaptive);
        assert_eq!(adaptive.to_u64(), Some(838102050));
    }

    #[test]
    fn test_large_karatsuba() {
        // 测试大数的Karatsuba算法
        let a = Mpz::from_str("2", 10).unwrap().pow_u32(100);
        let b = Mpz::from_str("3", 10).unwrap().pow_u32(100);

        let result = a.mul_karatsuba(&b);
        let expected = Mpz::from_str("6", 10).unwrap().pow_u32(100);

        assert_eq!(result, expected);
    }
}
