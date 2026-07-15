//! Mpz 的 trait 实现
//!
//! 为 Mpz 实现标准的 Rust traits，提升易用性

use super::core::Mpz;
use core::ops::{AddAssign, DivAssign, MulAssign, RemAssign, SubAssign};
use core::ops::{BitAndAssign, BitOrAssign, BitXorAssign, ShlAssign, ShrAssign};

// ============= 算术赋值 traits =============

impl AddAssign<&Mpz> for Mpz {
    fn add_assign(&mut self, rhs: &Mpz) {
        *self = self.add(rhs);
    }
}

impl AddAssign<Mpz> for Mpz {
    fn add_assign(&mut self, rhs: Mpz) {
        *self = self.add(&rhs);
    }
}

impl SubAssign<&Mpz> for Mpz {
    fn sub_assign(&mut self, rhs: &Mpz) {
        *self = self.sub(rhs);
    }
}

impl SubAssign<Mpz> for Mpz {
    fn sub_assign(&mut self, rhs: Mpz) {
        *self = self.sub(&rhs);
    }
}

impl MulAssign<&Mpz> for Mpz {
    fn mul_assign(&mut self, rhs: &Mpz) {
        *self = self.mul(rhs);
    }
}

impl MulAssign<Mpz> for Mpz {
    fn mul_assign(&mut self, rhs: Mpz) {
        *self = self.mul(&rhs);
    }
}

impl DivAssign<&Mpz> for Mpz {
    fn div_assign(&mut self, rhs: &Mpz) {
        *self = self.div(rhs).unwrap_or_else(|_| Mpz::new());
    }
}

impl DivAssign<Mpz> for Mpz {
    fn div_assign(&mut self, rhs: Mpz) {
        *self = self.div(&rhs).unwrap_or_else(|_| Mpz::new());
    }
}

impl RemAssign<&Mpz> for Mpz {
    fn rem_assign(&mut self, rhs: &Mpz) {
        *self = self.rem(rhs).unwrap_or_else(|_| Mpz::new());
    }
}

impl RemAssign<Mpz> for Mpz {
    fn rem_assign(&mut self, rhs: Mpz) {
        *self = self.rem(&rhs).unwrap_or_else(|_| Mpz::new());
    }
}

// ============= 位操作赋值 traits =============

impl ShlAssign<usize> for Mpz {
    fn shl_assign(&mut self, rhs: usize) {
        *self = self.shl(rhs);
    }
}

impl ShrAssign<usize> for Mpz {
    fn shr_assign(&mut self, rhs: usize) {
        *self = self.shr(rhs);
    }
}

impl BitAndAssign<&Mpz> for Mpz {
    fn bitand_assign(&mut self, rhs: &Mpz) {
        *self = self.bitwise_and(rhs);
    }
}

impl BitAndAssign<Mpz> for Mpz {
    fn bitand_assign(&mut self, rhs: Mpz) {
        *self = self.bitwise_and(&rhs);
    }
}

impl BitOrAssign<&Mpz> for Mpz {
    fn bitor_assign(&mut self, rhs: &Mpz) {
        *self = self.bitwise_or(rhs);
    }
}

impl BitOrAssign<Mpz> for Mpz {
    fn bitor_assign(&mut self, rhs: Mpz) {
        *self = self.bitwise_or(&rhs);
    }
}

impl BitXorAssign<&Mpz> for Mpz {
    fn bitxor_assign(&mut self, rhs: &Mpz) {
        *self = self.bitwise_xor(rhs);
    }
}

impl BitXorAssign<Mpz> for Mpz {
    fn bitxor_assign(&mut self, rhs: Mpz) {
        *self = self.bitwise_xor(&rhs);
    }
}

// ============= 与内置类型的赋值操作 =============

impl AddAssign<i64> for Mpz {
    fn add_assign(&mut self, rhs: i64) {
        *self = self.add(&Mpz::from_i64(rhs));
    }
}

impl AddAssign<u64> for Mpz {
    fn add_assign(&mut self, rhs: u64) {
        *self = self.add(&Mpz::from_u64(rhs));
    }
}

impl SubAssign<i64> for Mpz {
    fn sub_assign(&mut self, rhs: i64) {
        *self = self.sub(&Mpz::from_i64(rhs));
    }
}

impl SubAssign<u64> for Mpz {
    fn sub_assign(&mut self, rhs: u64) {
        *self = self.sub(&Mpz::from_u64(rhs));
    }
}

impl MulAssign<i64> for Mpz {
    fn mul_assign(&mut self, rhs: i64) {
        *self = self.mul(&Mpz::from_i64(rhs));
    }
}

impl MulAssign<u64> for Mpz {
    fn mul_assign(&mut self, rhs: u64) {
        *self = self.mul(&Mpz::from_u64(rhs));
    }
}

impl DivAssign<i64> for Mpz {
    fn div_assign(&mut self, rhs: i64) {
        *self = self.div(&Mpz::from_i64(rhs)).unwrap_or_else(|_| Mpz::new());
    }
}

impl DivAssign<u64> for Mpz {
    fn div_assign(&mut self, rhs: u64) {
        *self = self.div(&Mpz::from_u64(rhs)).unwrap_or_else(|_| Mpz::new());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_assign() {
        let mut a = Mpz::from_i64(42);
        let b = Mpz::from_i64(8);

        a += &b;
        assert_eq!(a.to_i64(), Some(50));

        a += 5i64;
        assert_eq!(a.to_i64(), Some(55));

        a += 5u64;
        assert_eq!(a.to_i64(), Some(60));
    }

    #[test]
    fn test_sub_assign() {
        let mut a = Mpz::from_i64(100);
        let b = Mpz::from_i64(30);

        a -= &b;
        assert_eq!(a.to_i64(), Some(70));

        a -= 20i64;
        assert_eq!(a.to_i64(), Some(50));
    }

    #[test]
    fn test_mul_assign() {
        let mut a = Mpz::from_i64(6);
        let b = Mpz::from_i64(7);

        a *= &b;
        assert_eq!(a.to_i64(), Some(42));

        a *= 2i64;
        assert_eq!(a.to_i64(), Some(84));
    }

    #[test]
    fn test_div_assign() {
        let mut a = Mpz::from_i64(84);
        let b = Mpz::from_i64(7);

        a /= &b;
        assert_eq!(a.to_i64(), Some(12));

        a /= 3i64;
        assert_eq!(a.to_i64(), Some(4));
    }

    #[test]
    fn test_bitwise_assign() {
        let mut a = Mpz::from_u64(0b1010);
        let b = Mpz::from_u64(0b1100);

        // AND
        let mut temp = a.clone();
        temp &= &b;
        assert_eq!(temp.to_u64(), Some(0b1000));

        // OR
        let mut temp = a.clone();
        temp |= &b;
        assert_eq!(temp.to_u64(), Some(0b1110));

        // XOR
        let mut temp = a.clone();
        temp ^= &b;
        assert_eq!(temp.to_u64(), Some(0b0110));

        // Shift
        a <<= 2;
        assert_eq!(a.to_u64(), Some(0b101000));

        a >>= 1;
        assert_eq!(a.to_u64(), Some(0b10100));
    }
}
