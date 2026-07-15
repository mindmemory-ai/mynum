//! Mpf 比较运算实现

use core::cmp::Ordering;

use super::core::Mpf;

impl Mpf {
    /// 绝对值比较
    pub fn cmp_abs(&self, other: &Mpf) -> Ordering {
        let (aligned_self, aligned_other) = self.align_exponents(other);
        aligned_self.mantissa().cmp(aligned_other.mantissa())
    }

    /// 检查是否相等
    pub fn is_eq(&self, other: &Mpf) -> bool {
        self.cmp(other) == Ordering::Equal
    }

    /// 检查是否小于
    pub fn lt(&self, other: &Mpf) -> bool {
        self.cmp(other) == Ordering::Less
    }

    /// 检查是否小于等于
    pub fn le(&self, other: &Mpf) -> bool {
        self.cmp(other) != Ordering::Greater
    }

    /// 检查是否大于
    pub fn gt(&self, other: &Mpf) -> bool {
        self.cmp(other) == Ordering::Greater
    }

    /// 检查是否大于等于
    pub fn ge(&self, other: &Mpf) -> bool {
        self.cmp(other) != Ordering::Less
    }

    /// 获取最大值
    pub fn max(&self, other: &Mpf) -> Mpf {
        if self.cmp(other) == Ordering::Greater {
            self.clone()
        } else {
            other.clone()
        }
    }

    /// 获取最小值
    pub fn min(&self, other: &Mpf) -> Mpf {
        if self.cmp(other) == Ordering::Less {
            self.clone()
        } else {
            other.clone()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mpz::Mpz;

    #[test]
    fn test_comparison() {
        let a = Mpf::from_mpz(Mpz::from_i64(15), 64);
        let b = Mpf::from_mpz(Mpz::from_i64(20), 64);
        let c = Mpf::from_mpz(Mpz::from_i64(15), 64);

        assert!(a < b);
        assert!(b > a);
        assert!(a == c);
        assert!(a <= c);
        assert!(a >= c);
    }

    #[test]
    fn test_equality() {
        let a = Mpf::from_mpz(Mpz::from_i64(15), 64);
        let b = Mpf::from_mpz(Mpz::from_i64(15), 64);
        let c = Mpf::from_mpz(Mpz::from_i64(20), 64);

        assert!(a.is_eq(&b));
        assert!(!a.is_eq(&c));
    }

    #[test]
    fn test_zero_comparison() {
        let zero = Mpf::new();
        let pos = Mpf::from_mpz(Mpz::from_i64(10), 64);
        let mut neg = Mpf::from_mpz(Mpz::from_i64(10), 64);
        neg.set_negative(true);

        assert!(zero == zero);
        assert!(zero < pos);
        assert!(zero > neg);
    }
}
