//! 类型转换实现

use super::core::Mpz;

// 从标准整数类型的转换
impl From<i32> for Mpz {
    fn from(n: i32) -> Self {
        Self::from_i64(n as i64)
    }
}

impl From<u32> for Mpz {
    fn from(n: u32) -> Self {
        Self::from_u64(n as u64)
    }
}

impl From<i64> for Mpz {
    fn from(n: i64) -> Self {
        Self::from_i64(n)
    }
}

impl From<u64> for Mpz {
    fn from(n: u64) -> Self {
        Self::from_u64(n)
    }
}

impl From<i128> for Mpz {
    fn from(n: i128) -> Self {
        Self::from_i128(n)
    }
}

impl From<u128> for Mpz {
    fn from(n: u128) -> Self {
        Self::from_u128(n)
    }
}

impl From<usize> for Mpz {
    fn from(n: usize) -> Self {
        Self::from_u64(n as u64)
    }
}

impl From<isize> for Mpz {
    fn from(n: isize) -> Self {
        Self::from_i64(n as i64)
    }
}

// 尝试转换到标准整数类型
impl TryFrom<&Mpz> for i32 {
    type Error = crate::error::Error;

    fn try_from(value: &Mpz) -> Result<Self, Self::Error> {
        value
            .to_i64()
            .and_then(|n| {
                if n >= i32::MIN as i64 && n <= i32::MAX as i64 {
                    Some(n as i32)
                } else {
                    None
                }
            })
            .ok_or(crate::error::Error::Overflow)
    }
}

impl TryFrom<&Mpz> for u32 {
    type Error = crate::error::Error;

    fn try_from(value: &Mpz) -> Result<Self, Self::Error> {
        value
            .to_u64()
            .and_then(|n| {
                if n <= u32::MAX as u64 {
                    Some(n as u32)
                } else {
                    None
                }
            })
            .ok_or(crate::error::Error::Overflow)
    }
}

impl TryFrom<&Mpz> for i64 {
    type Error = crate::error::Error;

    fn try_from(value: &Mpz) -> Result<Self, Self::Error> {
        value.to_i64().ok_or(crate::error::Error::Overflow)
    }
}

impl TryFrom<&Mpz> for u64 {
    type Error = crate::error::Error;

    fn try_from(value: &Mpz) -> Result<Self, Self::Error> {
        value.to_u64().ok_or(crate::error::Error::Overflow)
    }
}

impl TryFrom<&Mpz> for i128 {
    type Error = crate::error::Error;

    fn try_from(value: &Mpz) -> Result<Self, Self::Error> {
        value.to_i128().ok_or(crate::error::Error::Overflow)
    }
}

impl TryFrom<&Mpz> for u128 {
    type Error = crate::error::Error;

    fn try_from(value: &Mpz) -> Result<Self, Self::Error> {
        value.to_u128().ok_or(crate::error::Error::Overflow)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_conversion() {
        let a: Mpz = 42i32.into();
        assert_eq!(a.to_i64(), Some(42));

        let b: Mpz = 42u64.into();
        assert_eq!(b.to_u64(), Some(42));
    }

    #[test]
    fn test_try_into_conversion() {
        let a = Mpz::from_i64(42);
        let result: Result<i32, _> = (&a).try_into();
        assert_eq!(result.unwrap(), 42i32);

        let b = Mpz::from_i64(i64::MAX);
        let result: Result<i32, _> = (&b).try_into();
        assert!(result.is_err());
    }
}
