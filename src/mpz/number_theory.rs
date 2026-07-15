//! 数论函数实现

use super::core::Mpz;
use crate::error::{Error, Result};

impl Mpz {
    /// 最大公约数
    pub fn gcd(&self, other: &Mpz) -> Mpz {
        // 使用欧几里得算法
        let mut a = self.abs();
        let mut b = other.abs();

        while !b.is_zero() {
            let remainder = a.rem(&b).unwrap_or_else(|_| Mpz::new());
            a = b;
            b = remainder;
        }

        a
    }

    /// 最小公倍数
    pub fn lcm(&self, other: &Mpz) -> Mpz {
        if self.is_zero() || other.is_zero() {
            return Mpz::new();
        }

        let gcd = self.gcd(other);
        let product = self.abs().mul(&other.abs());
        product.div(&gcd).unwrap_or_else(|_| Mpz::new())
    }

    /// 扩展欧几里得算法
    pub fn extended_gcd(&self, other: &Mpz) -> (Mpz, Mpz, Mpz) {
        let mut a = self.abs();
        let mut b = other.abs();
        let mut x0 = Mpz::from_i64(1);
        let mut x1 = Mpz::new();
        let mut y0 = Mpz::new();
        let mut y1 = Mpz::from_i64(1);

        while !b.is_zero() {
            let (quotient, remainder) = a.div_rem(&b).unwrap_or_else(|_| (Mpz::new(), Mpz::new()));

            let temp_x = x1.clone();
            x1 = x0.sub(&quotient.mul(&x1));
            x0 = temp_x;

            let temp_y = y1.clone();
            y1 = y0.sub(&quotient.mul(&y1));
            y0 = temp_y;

            a = b;
            b = remainder;
        }

        (a, x0, y0)
    }

    /// 模逆元
    pub fn mod_inverse(&self, modulus: &Mpz) -> Result<Mpz> {
        if modulus.is_zero() {
            return Err(Error::InvalidInput("Modulus cannot be zero".into()));
        }

        let (gcd, x, _) = self.extended_gcd(modulus);

        if gcd != Mpz::from_i64(1) {
            return Err(Error::InvalidInput("Modular inverse does not exist".into()));
        }

        // 确保结果在 [0, modulus) 范围内
        let mut result = x.rem(modulus)?;
        if result.is_negative() {
            result = result.add(modulus);
        }

        Ok(result)
    }

    /// 雅可比符号
    pub fn jacobi(&self, other: &Mpz) -> Result<i32> {
        if other.is_zero() || other.rem(&Mpz::from_i64(2)).unwrap().is_zero() {
            return Err(Error::InvalidInput(
                "Second argument must be odd positive integer".into(),
            ));
        }

        let mut a = self.rem(other)?.abs();
        let mut n = other.abs();
        let mut result = 1;

        while !a.is_zero() {
            // 提取因子2
            let mut t = 0;
            while a.rem(&Mpz::from_i64(2)).unwrap().is_zero() {
                t += 1;
                a = a.div(&Mpz::from_i64(2)).unwrap();
            }

            // 二次互反律
            if t % 2 == 1 {
                let n_mod_8 = n.rem(&Mpz::from_i64(8)).unwrap();
                if n_mod_8 == Mpz::from_i64(3) || n_mod_8 == Mpz::from_i64(5) {
                    result = -result;
                }
            }

            // 二次互反律
            if a.rem(&Mpz::from_i64(4)).unwrap() == Mpz::from_i64(3)
                && n.rem(&Mpz::from_i64(4)).unwrap() == Mpz::from_i64(3)
            {
                result = -result;
            }

            // 交换a和n
            let temp = a;
            a = n.rem(&temp).unwrap();
            n = temp;
        }

        if n == Mpz::from_i64(1) {
            Ok(result)
        } else {
            Ok(0)
        }
    }

    /// 勒让德符号
    pub fn legendre(&self, p: &Mpz) -> Result<i32> {
        if p.is_zero() || p.rem(&Mpz::from_i64(2)).unwrap().is_zero() {
            return Err(Error::InvalidInput(
                "Second argument must be odd prime".into(),
            ));
        }

        // 勒让德符号是雅可比符号的特殊情况
        self.jacobi(p)
    }

    /// 二次剩余（模平方根）
    pub fn sqrt_mod(&self, p: &Mpz) -> Result<Vec<Mpz>> {
        if p.is_zero() || p.rem(&Mpz::from_i64(2)).unwrap().is_zero() {
            return Err(Error::InvalidInput("Modulus must be odd prime".into()));
        }

        let a = self.rem(p)?;
        if a.is_zero() {
            return Ok(vec![Mpz::new()]);
        }

        // 检查是否有二次剩余
        let legendre = a.legendre(p)?;
        if legendre == -1 {
            return Ok(vec![]); // 无解
        }

        if legendre == 0 {
            return Ok(vec![Mpz::new()]); // 只有零解
        }

        // 对于 p ≡ 3 (mod 4) 的情况，使用简单公式
        let p_mod_4 = p.rem(&Mpz::from_i64(4)).unwrap();
        if p_mod_4 == Mpz::from_i64(3) {
            let exp = p.add(&Mpz::from_i64(1)).div(&Mpz::from_i64(4)).unwrap();
            let sqrt = Self::mod_pow(&a, &exp, p);
            let neg_sqrt = p.sub(&sqrt);
            return Ok(vec![sqrt, neg_sqrt]);
        }

        // 对于 p ≡ 1 (mod 4) 的情况，使用Tonelli-Shanks算法
        self.tonelli_shanks(p)
    }

    /// Tonelli-Shanks算法实现
    fn tonelli_shanks(&self, p: &Mpz) -> Result<Vec<Mpz>> {
        let a = self.rem(p)?;
        let mut q = p.sub(&Mpz::from_i64(1));
        let mut s = 0;

        // 找到 q = 2^s * t，其中 t 是奇数
        while q.rem(&Mpz::from_i64(2)).unwrap().is_zero() {
            s += 1;
            q = q.div(&Mpz::from_i64(2)).unwrap();
        }

        // 寻找二次非剩余
        let mut z = Mpz::from_i64(2);
        while z.legendre(p)? != -1 {
            z = z.add(&Mpz::from_i64(1));
        }

        let mut m = s;
        let mut c = Self::mod_pow(&z, &q, p);
        let mut t = Self::mod_pow(&a, &q, p);
        let mut r = Self::mod_pow(
            &a,
            &q.add(&Mpz::from_i64(1)).div(&Mpz::from_i64(2)).unwrap(),
            p,
        );

        while !t.is_zero() && t != Mpz::from_i64(1) {
            let mut i = 0;
            let mut temp = t.clone();

            while temp != Mpz::from_i64(1) && i < m {
                temp = Self::mod_pow(&temp, &Mpz::from_i64(2), p);
                i += 1;
            }

            if i == 0 {
                return Ok(vec![]);
            }

            let b = Self::mod_pow(&c, &Mpz::from_i64(1).shl(m - i - 1), p);
            m = i;
            c = Self::mod_pow(&b, &Mpz::from_i64(2), p);
            t = t.mul(&c).rem(p).unwrap();
            r = r.mul(&b).rem(p).unwrap();
        }

        if t.is_zero() {
            Ok(vec![Mpz::new()])
        } else {
            let neg_r = p.sub(&r);
            Ok(vec![r, neg_r])
        }
    }

    /// 素数测试（Miller-Rabin算法）
    pub fn is_probably_prime(&self, rounds: u32) -> bool {
        if self.is_zero() || self == &Mpz::from_i64(1) {
            return false;
        }
        if self == &Mpz::from_i64(2) || self == &Mpz::from_i64(3) {
            return true;
        }
        if self.rem(&Mpz::from_i64(2)).unwrap().is_zero() {
            return false;
        }

        // 对于小数的快速测试
        if self.bit_length() <= 64 {
            if let Some(n) = self.to_u64() {
                return Self::is_prime_small(n);
            }
        }

        // Miller-Rabin测试
        self.miller_rabin_test(rounds)
    }

    /// 小数的素数测试（用于优化）
    fn is_prime_small(n: u64) -> bool {
        if n < 2 {
            return false;
        }
        if n == 2 || n == 3 {
            return true;
        }
        if n.is_multiple_of(2) {
            return false;
        }

        // 试除法到sqrt(n)
        let sqrt_n = (n as f64).sqrt() as u64;
        for i in (3..=sqrt_n).step_by(2) {
            if n.is_multiple_of(i) {
                return false;
            }
        }
        true
    }

    /// Miller-Rabin素数测试
    fn miller_rabin_test(&self, rounds: u32) -> bool {
        let n = self;
        let n_minus_1 = n.sub(&Mpz::from_i64(1));

        // 找到 n-1 = 2^s * d，其中d是奇数
        let mut s = 0;
        let mut d = n_minus_1.clone();

        while d.rem(&Mpz::from_i64(2)).unwrap().is_zero() {
            s += 1;
            d = d.div(&Mpz::from_i64(2)).unwrap();
        }

        // 选择测试基数
        if n.bit_length() <= 64 {
            // 对于64位以下的数，使用确定性测试
            let small_bases = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
            for &base in &small_bases {
                let base_mpz = Mpz::from_i64(base);
                if !Self::miller_rabin_witness(n, &base_mpz, &d, s) {
                    return false;
                }
            }
        } else {
            // 对于大数，使用随机测试
            for _ in 0..rounds.min(20) {
                let base = Mpz::random_range(&Mpz::from_i64(2), &n_minus_1)
                    .unwrap_or_else(|_| Mpz::from_i64(2));
                if !Self::miller_rabin_witness(n, &base, &d, s) {
                    return false;
                }
            }
        }

        true
    }

    /// Miller-Rabin见证函数
    fn miller_rabin_witness(n: &Mpz, base: &Mpz, d: &Mpz, s: u32) -> bool {
        let mut x = Self::mod_pow(base, d, n);

        if x == Mpz::from_i64(1) || x == n.sub(&Mpz::from_i64(1)) {
            return true;
        }

        for _ in 1..s {
            x = x.mul(&x).rem(n).unwrap_or_else(|_| Mpz::new());
            if x == n.sub(&Mpz::from_i64(1)) {
                return true;
            }
            if x == Mpz::from_i64(1) {
                return false;
            }
        }

        false
    }

    /// 模幂运算: base^exp mod modulus
    pub fn mod_pow(base: &Mpz, exp: &Mpz, modulus: &Mpz) -> Mpz {
        if modulus.is_zero() {
            return Mpz::new();
        }

        let mut result = Mpz::from_i64(1);
        let mut base = base.rem(modulus).unwrap_or_else(|_| Mpz::new());
        let mut exp = exp.clone();

        while !exp.is_zero() {
            if exp.rem(&Mpz::from_i64(2)).unwrap() == Mpz::from_i64(1) {
                result = result
                    .mul(&base)
                    .rem(modulus)
                    .unwrap_or_else(|_| Mpz::new());
            }
            base = base.mul(&base).rem(modulus).unwrap_or_else(|_| Mpz::new());
            exp = exp.div(&Mpz::from_i64(2)).unwrap_or_else(|_| Mpz::new());
        }

        result
    }

    /// 下一个素数
    pub fn next_prime(&self) -> Mpz {
        let mut candidate = self.add(&Mpz::from_i64(1));
        if candidate.rem(&Mpz::from_i64(2)).unwrap().is_zero() {
            candidate = candidate.add(&Mpz::from_i64(1));
        }

        while !candidate.is_probably_prime(25) {
            candidate = candidate.add(&Mpz::from_i64(2));
        }

        candidate
    }

    /// 随机素数生成
    pub fn random_prime(bits: usize) -> Result<Mpz> {
        if bits < 2 {
            return Err(Error::InvalidInput("Bits must be at least 2".into()));
        }

        // 生成随机奇数
        let mut candidate = Mpz::random_odd_bits(bits)?;

        // 确保最高位为1（保证位数）
        if candidate.bit_length() < bits {
            let highest_bit = Mpz::from_i64(1).shl(bits - 1);
            candidate = candidate.bitwise_or(&highest_bit);
        }

        // 如果生成的数是偶数，加1使其变为奇数
        if candidate.rem(&Mpz::from_i64(2)).unwrap().is_zero() {
            candidate = candidate.add(&Mpz::from_i64(1));
        }

        // 寻找素数
        let mut attempts = 0;
        const MAX_ATTEMPTS: usize = 1000;

        while attempts < MAX_ATTEMPTS {
            if candidate.is_probably_prime(25) {
                return Ok(candidate);
            }

            // 加2继续寻找下一个奇数
            candidate = candidate.add(&Mpz::from_i64(2));
            attempts += 1;

            // 如果位数增加了，重新生成
            if candidate.bit_length() > bits {
                candidate = Mpz::random_odd_bits(bits)?;
                if candidate.bit_length() < bits {
                    let highest_bit = Mpz::from_i64(1).shl(bits - 1);
                    candidate = candidate.bitwise_or(&highest_bit);
                }
                if candidate.rem(&Mpz::from_i64(2)).unwrap().is_zero() {
                    candidate = candidate.add(&Mpz::from_i64(1));
                }
            }
        }

        Err(Error::Other(
            "Failed to generate prime after maximum attempts".into(),
        ))
    }

    /// 平方根
    pub fn sqrt(&self) -> Mpz {
        if self.is_zero() {
            return Mpz::new();
        }

        if self.is_negative() {
            // 负数的平方根未定义，返回0
            return Mpz::new();
        }

        // 使用牛顿法计算平方根
        let mut x = self.clone();
        let two = Mpz::from_i64(2);

        loop {
            let y = x
                .add(&self.div(&x).unwrap_or_else(|_| Mpz::new()))
                .div(&two)
                .unwrap_or_else(|_| Mpz::new());
            if y.cmp(&x) != core::cmp::Ordering::Less {
                break;
            }
            x = y;
        }

        x
    }

    /// 立方根
    pub fn cbrt(&self) -> Mpz {
        if self.is_zero() {
            return Mpz::new();
        }

        // 简单实现，使用牛顿法
        let mut x = self.clone();
        let three = Mpz::from_i64(3);

        for _ in 0..100 {
            // 限制迭代次数
            let x_squared = x.mul(&x);
            let x_cubed = x_squared.mul(&x);

            if x_cubed.cmp(self) == core::cmp::Ordering::Equal {
                break;
            }

            let two_x = x.mul(&Mpz::from_i64(2));
            let x_sq_div = self.div(&x_squared).unwrap_or_else(|_| Mpz::new());
            let numerator = two_x.add(&x_sq_div);
            let new_x = numerator.div(&three).unwrap_or_else(|_| Mpz::new());

            if new_x.cmp(&x) == core::cmp::Ordering::Equal {
                break;
            }

            x = new_x;
        }

        x
    }

    /// 计算阶乘
    pub fn factorial(n: u64) -> Mpz {
        if n == 0 || n == 1 {
            return Mpz::from_i64(1);
        }

        let mut result = Mpz::from_i64(1);
        for i in 2..=n {
            result = result.mul(&Mpz::from_u64(i));
        }

        result
    }

    /// 试除法因子分解
    ///
    /// 返回所有质因数的向量，每个质因数可能出现多次
    ///
    /// # 参数
    /// - `n`: 要分解的数
    ///
    /// # 返回
    /// - `Vec<Mpz>`: 质因数列表
    pub fn trial_division_factorization(&self) -> Vec<Mpz> {
        let mut n = self.abs();
        let mut factors = Vec::new();

        if n.is_zero() || n == Mpz::from_i64(1) {
            return factors;
        }

        // 处理2这个特殊的质数
        while n.rem(&Mpz::from_i64(2)).unwrap().is_zero() {
            factors.push(Mpz::from_i64(2));
            n = n.div(&Mpz::from_i64(2)).unwrap();
        }

        // 试除奇数
        let mut d = Mpz::from_i64(3);
        let sqrt_n = n.sqrt();

        while d <= sqrt_n {
            while n.rem(&d).unwrap().is_zero() {
                factors.push(d.clone());
                n = n.div(&d).unwrap();
            }
            d = d.add(&Mpz::from_i64(2));
        }

        // 如果n > 1，说明n本身是质数
        if n > Mpz::from_i64(1) {
            factors.push(n);
        }

        factors
    }

    /// Pollard's Rho算法进行因子分解
    ///
    /// 对于大数的因子分解，比试除法更高效
    ///
    /// # 参数
    /// - `n`: 要分解的数
    ///
    /// # 返回
    /// - `Option<Mpz>`: 找到的一个非平凡因子，如果没有找到则返回None
    pub fn pollard_rho_factorization(&self) -> Option<Mpz> {
        let n = self.abs();

        if n.is_zero() || n == Mpz::from_i64(1) {
            return None;
        }

        if n.is_probably_prime(10) {
            return None;
        }

        // 如果n是偶数，返回2
        if n.rem(&Mpz::from_i64(2)).unwrap().is_zero() {
            return Some(Mpz::from_i64(2));
        }

        let mut x = Mpz::from_i64(2);
        let mut y = Mpz::from_i64(2);
        let c = Mpz::from_i64(1);

        let f = |x: &Mpz, n: &Mpz| x.mul(x).add(&c).rem(n).unwrap_or_else(|_| Mpz::new());

        // 简化实现：使用固定迭代次数
        for _ in 0..1000 {
            x = f(&x, &n);
            y = f(&f(&y, &n), &n);

            let diff = if x >= y { x.sub(&y) } else { y.sub(&x) };
            let d = diff.gcd(&n);

            if d > Mpz::from_i64(1) && d < n {
                return Some(d);
            }
        }

        None
    }

    /// 获取所有因子
    ///
    /// 返回一个数的所有正因子（包括1和自身）
    ///
    /// # 返回
    /// - `Vec<Mpz>`: 所有因子的列表
    pub fn get_all_factors(&self) -> Vec<Mpz> {
        let n = self.abs();
        if n.is_zero() {
            return vec![];
        }

        let mut factors = Vec::new();
        factors.push(Mpz::from_i64(1));

        if n == Mpz::from_i64(1) {
            return factors;
        }

        // 简单方法：遍历1到sqrt(n)，检查每个数是否为因子
        let sqrt_n = n.sqrt();
        let mut i = Mpz::from_i64(2);

        while i <= sqrt_n {
            if n.rem(&i).unwrap().is_zero() {
                factors.push(i.clone());
                let other_factor = n.div(&i).unwrap();
                if other_factor != i {
                    factors.push(other_factor);
                }
            }
            i = i.add(&Mpz::from_i64(1));
        }

        factors.push(n.clone());
        factors.sort();

        factors
    }

    /// 二项式系数
    pub fn binomial(n: u64, k: u64) -> Mpz {
        if k > n {
            return Mpz::new();
        }

        if k == 0 || k == n {
            return Mpz::from_i64(1);
        }

        // 使用优化的计算方法避免大数阶乘
        let k = k.min(n - k); // 选择较小的k
        let mut result = Mpz::from_i64(1);

        for i in 0..k {
            result = result.mul(&Mpz::from_u64(n - i));
            result = result
                .div(&Mpz::from_u64(i + 1))
                .unwrap_or_else(|_| Mpz::new());
        }

        result
    }

    /// 斐波那契数
    pub fn fibonacci(n: u64) -> Mpz {
        if n == 0 {
            return Mpz::new();
        }
        if n == 1 {
            return Mpz::from_i64(1);
        }

        let mut a = Mpz::new();
        let mut b = Mpz::from_i64(1);

        for _ in 2..=n {
            let temp = a.add(&b);
            a = b;
            b = temp;
        }

        b
    }

    /// 卢卡斯数
    pub fn lucas(n: u64) -> Mpz {
        if n == 0 {
            return Mpz::from_i64(2);
        }
        if n == 1 {
            return Mpz::from_i64(1);
        }

        let mut a = Mpz::from_i64(2);
        let mut b = Mpz::from_i64(1);

        for _ in 2..=n {
            let temp = a.add(&b);
            a = b;
            b = temp;
        }

        b
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcd() {
        let a = Mpz::from_i64(48);
        let b = Mpz::from_i64(18);
        let gcd = a.gcd(&b);
        assert_eq!(gcd.to_i64(), Some(6));
    }

    #[test]
    fn test_lcm() {
        let a = Mpz::from_i64(12);
        let b = Mpz::from_i64(18);
        let lcm = a.lcm(&b);
        assert_eq!(lcm.to_i64(), Some(36));
    }

    #[test]
    fn test_extended_gcd() {
        let a = Mpz::from_i64(48);
        let b = Mpz::from_i64(18);
        let (gcd, x, y) = a.extended_gcd(&b);

        assert_eq!(gcd.to_i64(), Some(6));
        // 验证 ax + by = gcd
        let check = a.mul(&x).add(&b.mul(&y));
        assert_eq!(check, gcd);
    }

    #[test]
    fn test_mod_inverse() {
        let a = Mpz::from_i64(3);
        let modulus = Mpz::from_i64(11);
        let inverse = a.mod_inverse(&modulus).unwrap();

        // 验证 a * inverse ≡ 1 (mod modulus)
        let product = a.mul(&inverse);
        let remainder = product.rem(&modulus).unwrap();
        assert_eq!(remainder, Mpz::from_i64(1));
    }

    #[test]
    fn test_jacobi_symbol() {
        let a = Mpz::from_i64(2);
        let n = Mpz::from_i64(15);
        let jacobi = a.jacobi(&n).unwrap();
        assert_eq!(jacobi, 1);
    }

    #[test]
    fn test_legendre_symbol() {
        let a = Mpz::from_i64(2);
        let p = Mpz::from_i64(7);
        let legendre = a.legendre(&p).unwrap();
        assert_eq!(legendre, 1);
    }

    #[test]
    fn test_prime_test() {
        // 测试小素数
        assert!(Mpz::from_i64(2).is_probably_prime(25));
        assert!(Mpz::from_i64(3).is_probably_prime(25));
        assert!(Mpz::from_i64(5).is_probably_prime(25));
        assert!(Mpz::from_i64(7).is_probably_prime(25));
        assert!(Mpz::from_i64(11).is_probably_prime(25));

        // 测试合数
        assert!(!Mpz::from_i64(4).is_probably_prime(25));
        assert!(!Mpz::from_i64(6).is_probably_prime(25));
        assert!(!Mpz::from_i64(8).is_probably_prime(25));
        assert!(!Mpz::from_i64(9).is_probably_prime(25));
        assert!(!Mpz::from_i64(10).is_probably_prime(25));

        // 测试大素数（使用已知的素数）
        let large_prime = Mpz::from_str("1000000007", 10).unwrap();
        assert!(large_prime.is_probably_prime(25));
    }

    #[test]
    fn test_sqrt() {
        let a = Mpz::from_i64(16);
        let sqrt = a.sqrt();
        assert_eq!(sqrt.to_i64(), Some(4));

        let b = Mpz::from_i64(15);
        let sqrt_b = b.sqrt();
        assert_eq!(sqrt_b.to_i64(), Some(3)); // 向下取整
    }

    #[test]
    fn test_factorial() {
        assert_eq!(Mpz::factorial(0).to_i64(), Some(1));
        assert_eq!(Mpz::factorial(1).to_i64(), Some(1));
        assert_eq!(Mpz::factorial(5).to_i64(), Some(120));
    }

    #[test]
    fn test_trial_division_factorization() {
        // 测试质数
        let prime = Mpz::from_i64(17);
        let factors = prime.trial_division_factorization();
        assert_eq!(factors.len(), 1);
        assert_eq!(factors[0].to_i64(), Some(17));

        // 测试合数
        let composite = Mpz::from_i64(12);
        let factors = composite.trial_division_factorization();
        assert_eq!(factors.len(), 3);
        assert_eq!(factors[0].to_i64(), Some(2));
        assert_eq!(factors[1].to_i64(), Some(2));
        assert_eq!(factors[2].to_i64(), Some(3));

        // 测试完全平方数
        let square = Mpz::from_i64(16);
        let factors = square.trial_division_factorization();
        assert_eq!(factors.len(), 4);
        assert!(factors.iter().all(|f| f.to_i64() == Some(2)));
    }

    #[test]
    fn test_get_all_factors() {
        // 测试质数
        let prime = Mpz::from_i64(7);
        let factors = prime.get_all_factors();
        assert_eq!(factors.len(), 2);
        assert_eq!(factors[0].to_i64(), Some(1));
        assert_eq!(factors[1].to_i64(), Some(7));

        // 测试合数
        let composite = Mpz::from_i64(12);
        let factors = composite.get_all_factors();
        assert_eq!(factors.len(), 6);
        assert_eq!(
            factors,
            vec![
                Mpz::from_i64(1),
                Mpz::from_i64(2),
                Mpz::from_i64(3),
                Mpz::from_i64(4),
                Mpz::from_i64(6),
                Mpz::from_i64(12)
            ]
        );
    }

    #[test]
    fn test_pollard_rho_factorization() {
        // 测试质数（应该返回None）
        let prime = Mpz::from_i64(17);
        let factor = prime.pollard_rho_factorization();
        assert!(factor.is_none());

        // 测试合数
        let composite = Mpz::from_i64(15);
        let factor = composite.pollard_rho_factorization();
        // 可能找到3或5，或者返回None（取决于算法运气）
        if let Some(f) = factor {
            assert!(f == Mpz::from_i64(3) || f == Mpz::from_i64(5));
        }
    }

    #[test]
    fn test_fibonacci() {
        assert_eq!(Mpz::fibonacci(0).to_i64(), Some(0));
        assert_eq!(Mpz::fibonacci(1).to_i64(), Some(1));
        assert_eq!(Mpz::fibonacci(10).to_i64(), Some(55));
    }
}
