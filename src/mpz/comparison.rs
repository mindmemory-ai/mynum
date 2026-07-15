//! 比较运算实现

use super::core::Mpz;
use core::cmp::Ordering;

impl Mpz {
    /// 比较函数
    pub fn compare(&self, other: &Mpz) -> Ordering {
        // 处理零的特殊情况
        if self.is_zero() && other.is_zero() {
            return Ordering::Equal;
        }

        // 比较符号
        match (self.is_negative(), other.is_negative()) {
            (false, true) => Ordering::Greater,    // 正数 > 负数
            (true, false) => Ordering::Less,       // 负数 < 正数
            (false, false) => self.cmp_abs(other), // 都是正数，比较绝对值
            (true, true) => other.cmp_abs(self),   // 都是负数，反向比较绝对值
        }
    }

    /// 绝对值比较（内部使用）
    pub(crate) fn cmp_abs(&self, other: &Mpz) -> Ordering {
        let self_limbs = self.limbs();
        let other_limbs = other.limbs();

        // 首先比较长度
        match self_limbs.len().cmp(&other_limbs.len()) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => {
                // 长度相等，从高位开始比较
                for i in (0..self_limbs.len()).rev() {
                    match self_limbs[i].cmp(&other_limbs[i]) {
                        Ordering::Equal => continue,
                        other => return other,
                    }
                }
                Ordering::Equal
            }
        }
    }
}

impl PartialOrd for Mpz {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(Ord::cmp(self, other))
    }
}

impl Ord for Mpz {
    fn cmp(&self, other: &Self) -> Ordering {
        self.compare(other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_comparison() {
        let a = Mpz::from_i64(100);
        let b = Mpz::from_i64(200);
        let c = Mpz::from_i64(-100);

        assert_eq!(a.cmp(&b), Ordering::Less);
        assert_eq!(b.cmp(&a), Ordering::Greater);
        assert_eq!(a.cmp(&a), Ordering::Equal);
        assert_eq!(a.cmp(&c), Ordering::Greater);
        assert_eq!(c.cmp(&a), Ordering::Less);
    }

    #[test]
    fn test_trait_comparison() {
        let a = Mpz::from_i64(100);
        let b = Mpz::from_i64(200);

        assert!(a < b);
        assert!(b > a);
        assert!(a == a);
        assert!(a <= b);
        assert!(b >= a);
    }
}
