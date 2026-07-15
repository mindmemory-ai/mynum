//! 复数算术运算模块

use super::core::Complex;
use crate::error::Result;
use crate::mpf::Mpf;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::time::Instant;

impl Complex {
    /// 模长 |z| = √(a² + b²)
    pub fn magnitude(&self) -> Result<Mpf> {
        let real_sq = self.real().mul(self.real());
        let imag_sq = self.imaginary().mul(self.imaginary());
        let sum = real_sq.add(&imag_sq);
        sum.sqrt()
    }

    /// 模长的平方 |z|² = a² + b²
    pub fn magnitude_squared(&self) -> Result<Mpf> {
        let real_sq = self.real().mul(self.real());
        let imag_sq = self.imaginary().mul(self.imaginary());
        Ok(real_sq.add(&imag_sq))
    }

    /// 幅角 arg(z) = atan2(b, a)
    pub fn argument(&self) -> Result<Mpf> {
        if self.is_zero() {
            return Err(crate::error::Error::DomainError(
                "Argument of zero is undefined".into(),
            ));
        }

        if self.is_real() {
            if self.real().is_negative() {
                return Ok(Mpf::pi(self.precision()));
            } else {
                return Ok(Mpf::new());
            }
        }

        if self.is_imaginary() {
            if self.imaginary().is_negative() {
                let pi_val = Mpf::pi(self.precision());
                let two = Mpf::from_i64(2, self.precision());
                return Ok(pi_val.div(&two)?.neg());
            } else {
                let pi_val = Mpf::pi(self.precision());
                let two = Mpf::from_i64(2, self.precision());
                return pi_val.div(&two);
            }
        }

        // 使用 atan2(b, a) 计算幅角
        // 注意：Mpf 没有 atan2 方法，我们需要手动计算
        let real_val = self.real();
        let imag_val = self.imaginary();

        if real_val.is_zero() {
            if imag_val.is_negative() {
                let pi_val = Mpf::pi(self.precision());
                let two = Mpf::from_i64(2, self.precision());
                Ok(pi_val.div(&two)?.neg())
            } else {
                let pi_val = Mpf::pi(self.precision());
                let two = Mpf::from_i64(2, self.precision());
                Ok(pi_val.div(&two)?)
            }
        } else {
            // 使用 atan(b/a) 然后根据象限调整
            let ratio = imag_val.div(real_val)?;
            let atan_result = ratio.atan()?;

            if real_val.is_negative() {
                if imag_val.is_negative() {
                    // 第三象限
                    Ok(atan_result.sub(&Mpf::pi(self.precision())))
                } else {
                    // 第二象限
                    Ok(atan_result.add(&Mpf::pi(self.precision())))
                }
            } else {
                // 第一或第四象限
                Ok(atan_result)
            }
        }
    }

    /// 共轭复数 a-bi
    pub fn conjugate(&self) -> Self {
        Self::from_real_imag(self.real().clone(), self.imaginary().neg())
    }

    /// 倒数 1/z
    pub fn reciprocal(&self) -> Result<Self> {
        if self.is_zero() {
            return Err(crate::error::Error::DomainError("Division by zero".into()));
        }

        let mag_sq = self.magnitude_squared()?;
        let conj = self.conjugate();

        Ok(Self::from_real_imag(
            conj.real().div(&mag_sq)?,
            conj.imaginary().div(&mag_sq)?,
        ))
    }

    /// 加法运算
    pub fn add(&self, other: &Complex) -> Result<Self> {
        let start = Instant::now();

        let result = Self::from_real_imag(
            self.real().add(other.real()),
            self.imaginary().add(other.imaginary()),
        );

        let _computation_time = start.elapsed().as_micros() as u64;

        // 检查缓存配置
        if crate::complex::ComplexConfig::is_caching_enabled() {
            // 这里可以添加缓存逻辑
        }

        Ok(result)
    }

    /// 减法运算
    pub fn sub(&self, other: &Complex) -> Result<Self> {
        let start = Instant::now();

        let result = Self::from_real_imag(
            self.real().sub(other.real()),
            self.imaginary().sub(other.imaginary()),
        );

        let _computation_time = start.elapsed().as_micros() as u64;

        Ok(result)
    }

    /// 乘法运算 (a+bi)(c+di) = (ac-bd) + (ad+bc)i
    pub fn mul(&self, other: &Complex) -> Result<Self> {
        let start = Instant::now();

        let ac = self.real().mul(other.real());
        let bd = self.imaginary().mul(other.imaginary());
        let ad = self.real().mul(other.imaginary());
        let bc = self.imaginary().mul(other.real());

        let real_part = ac.sub(&bd);
        let imag_part = ad.add(&bc);

        let result = Self::from_real_imag(real_part, imag_part);

        let _computation_time = start.elapsed().as_micros() as u64;

        Ok(result)
    }

    /// 除法运算 (a+bi)/(c+di) = [(ac+bd)/(c²+d²)] + [(bc-ad)/(c²+d²)]i
    pub fn div(&self, other: &Complex) -> Result<Self> {
        if other.is_zero() {
            return Err(crate::error::Error::DomainError("Division by zero".into()));
        }

        let start = Instant::now();

        let denom = other.magnitude_squared()?;

        let ac = self.real().mul(other.real());
        let bd = self.imaginary().mul(other.imaginary());
        let bc = self.imaginary().mul(other.real());
        let ad = self.real().mul(other.imaginary());

        let real_part = ac.add(&bd).div(&denom);
        let imag_part = bc.sub(&ad).div(&denom);

        let result = Self::from_real_imag(real_part?, imag_part?);

        let _computation_time = start.elapsed().as_micros() as u64;

        Ok(result)
    }

    /// 幂运算 z^n
    pub fn pow(&self, n: i64) -> Result<Self> {
        if n == 0 {
            return Ok(Self::from_real(Mpf::from_i64(1, self.precision())));
        }

        if n == 1 {
            return Ok(self.clone());
        }

        if n < 0 {
            let reciprocal = self.reciprocal()?;
            return reciprocal.pow(-n);
        }

        // 使用快速幂算法
        let mut result = Self::from_real(Mpf::from_i64(1, self.precision()));
        let mut base = self.clone();
        let mut exponent = n;

        while exponent > 0 {
            if exponent & 1 == 1 {
                result = result.mul(&base)?;
            }
            base = base.clone().mul(&base)?;
            exponent >>= 1;
        }

        Ok(result)
    }

    /// 平方 z²
    pub fn square(&self) -> Result<Self> {
        self.mul(self)
    }

    /// 立方 z³
    pub fn cube(&self) -> Result<Self> {
        self.square()?.mul(self)
    }

    /// 绝对值 |z|
    pub fn abs(&self) -> Result<Mpf> {
        self.magnitude()
    }

    /// 归一化（单位化）z/|z|
    pub fn normalize(&self) -> Result<Self> {
        let mag = self.magnitude()?;
        if mag.is_zero() {
            return Err(crate::error::Error::DomainError(
                "Cannot normalize zero".into(),
            ));
        }

        Ok(Self::from_real_imag(
            self.real().div(&mag)?,
            self.imaginary().div(&mag)?,
        ))
    }

    /// 投影到单位圆上
    pub fn project_to_unit_circle(&self) -> Result<Self> {
        if self.is_zero() {
            return Ok(Self::new());
        }

        let mag = self.magnitude()?;
        let one = Mpf::from_i64(1, self.precision());
        if mag == one {
            return Ok(self.clone());
        }

        self.normalize()
    }

    /// 计算两个复数的距离
    pub fn distance(&self, other: &Complex) -> Result<Mpf> {
        let diff = self.sub(other)?;
        diff.magnitude()
    }

    /// 计算两个复数的内积（实部乘积 + 虚部乘积）
    pub fn dot_product(&self, other: &Complex) -> Result<Mpf> {
        let real_product = self.real().mul(other.real());
        let imag_product = self.imaginary().mul(other.imaginary());
        Ok(real_product.add(&imag_product))
    }

    /// 计算两个复数的外积（实部乘积 - 虚部乘积）
    pub fn cross_product(&self, other: &Complex) -> Result<Mpf> {
        let real_product = self.real().mul(other.real());
        let imag_product = self.imaginary().mul(other.imaginary());
        Ok(real_product.sub(&imag_product))
    }

    /// 计算复数的平均值
    pub fn average(complexes: &[Complex]) -> Result<Self> {
        if complexes.is_empty() {
            return Err(crate::error::Error::InvalidInput("Empty array".into()));
        }

        let mut sum_real = Mpf::new();
        let mut sum_imag = Mpf::new();
        let mut count = 0;

        for complex in complexes {
            sum_real = sum_real.add(complex.real());
            sum_imag = sum_imag.add(complex.imaginary());
            count += 1;
        }

        let count_mpf = Mpf::from_i64(count as i64, complexes[0].precision());

        Ok(Self::from_real_imag(
            sum_real.div(&count_mpf)?,
            sum_imag.div(&count_mpf)?,
        ))
    }

    /// 计算复数的加权平均值
    pub fn weighted_average(complexes: &[Complex], weights: &[Mpf]) -> Result<Self> {
        if complexes.len() != weights.len() || complexes.is_empty() {
            return Err(crate::error::Error::InvalidInput(
                "Invalid input lengths".into(),
            ));
        }

        let mut sum_real = Mpf::new();
        let mut sum_imag = Mpf::new();
        let mut sum_weights = Mpf::new();

        for (complex, weight) in complexes.iter().zip(weights.iter()) {
            sum_real = sum_real.add(&complex.real().mul(weight));
            sum_imag = sum_imag.add(&complex.imaginary().mul(weight));
            sum_weights = sum_weights.add(weight);
        }

        if sum_weights.is_zero() {
            return Err(crate::error::Error::DomainError(
                "Sum of weights is zero".into(),
            ));
        }

        Ok(Self::from_real_imag(
            sum_real.div(&sum_weights)?,
            sum_imag.div(&sum_weights)?,
        ))
    }
}

// 实现标准运算符

impl Add for Complex {
    type Output = Result<Self>;

    fn add(self, other: Self) -> Self::Output {
        self.add(&other)
    }
}

impl Add<&Complex> for Complex {
    type Output = Result<Self>;

    fn add(self, other: &Complex) -> Self::Output {
        // 调用实例方法，避免递归
        let real = self.real().add(other.real());
        let imaginary = self.imaginary().add(other.imaginary());
        Ok(Self::from_real_imag(real, imaginary))
    }
}

impl Sub for Complex {
    type Output = Result<Self>;

    fn sub(self, other: Self) -> Self::Output {
        self.sub(&other)
    }
}

impl Sub<&Complex> for Complex {
    type Output = Result<Self>;

    fn sub(self, other: &Complex) -> Self::Output {
        // 调用实例方法，避免递归
        let real = self.real().sub(other.real());
        let imaginary = self.imaginary().sub(other.imaginary());
        Ok(Self::from_real_imag(real, imaginary))
    }
}

impl Mul for Complex {
    type Output = Result<Self>;

    fn mul(self, other: Self) -> Self::Output {
        self.mul(&other)
    }
}

impl Mul<&Complex> for Complex {
    type Output = Result<Self>;

    fn mul(self, other: &Complex) -> Self::Output {
        // 调用实例方法，避免递归
        let ac = self.real().mul(other.real());
        let bd = self.imaginary().mul(other.imaginary());
        let ad = self.real().mul(other.imaginary());
        let bc = self.imaginary().mul(other.real());

        let real_part = ac.sub(&bd);
        let imag_part = ad.add(&bc);

        Ok(Self::from_real_imag(real_part, imag_part))
    }
}

impl Div for Complex {
    type Output = Result<Self>;

    fn div(self, other: Self) -> Self::Output {
        self.div(&other)
    }
}

impl Div<&Complex> for Complex {
    type Output = Result<Self>;

    fn div(self, other: &Complex) -> Self::Output {
        // 调用实例方法，避免递归
        if other.is_zero() {
            return Err(crate::error::Error::DomainError("Division by zero".into()));
        }

        let denom = other.magnitude_squared()?;

        let ac = self.real().mul(other.real());
        let bd = self.imaginary().mul(other.imaginary());
        let bc = self.imaginary().mul(other.real());
        let ad = self.real().mul(other.imaginary());

        let real_part = ac.add(&bd).div(&denom);
        let imag_part = bc.sub(&ad).div(&denom);

        Ok(Self::from_real_imag(real_part?, imag_part?))
    }
}

impl Neg for Complex {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::from_real_imag(self.real().neg(), self.imaginary().neg())
    }
}

impl AddAssign for Complex {
    fn add_assign(&mut self, other: Self) {
        if let Ok(result) = self.clone().add(&other) {
            *self = result;
        }
    }
}

impl SubAssign for Complex {
    fn sub_assign(&mut self, other: Self) {
        if let Ok(result) = self.clone().sub(&other) {
            *self = result;
        }
    }
}

impl MulAssign for Complex {
    fn mul_assign(&mut self, other: Self) {
        if let Ok(result) = self.clone().mul(&other) {
            *self = result;
        }
    }
}

impl DivAssign for Complex {
    fn div_assign(&mut self, other: Self) {
        if let Ok(result) = self.clone().div(&other) {
            *self = result;
        }
    }
}

// 与实数的运算

impl Add<Mpf> for Complex {
    type Output = Result<Self>;

    fn add(self, other: Mpf) -> Self::Output {
        Ok(Self::from_real_imag(
            self.real().add(&other),
            self.imaginary().clone(),
        ))
    }
}

impl Sub<Mpf> for Complex {
    type Output = Result<Self>;

    fn sub(self, other: Mpf) -> Self::Output {
        Ok(Self::from_real_imag(
            self.real().sub(&other),
            self.imaginary().clone(),
        ))
    }
}

impl Mul<Mpf> for Complex {
    type Output = Result<Self>;

    fn mul(self, other: Mpf) -> Self::Output {
        Ok(Self::from_real_imag(
            self.real().mul(&other),
            self.imaginary().mul(&other),
        ))
    }
}

impl Div<Mpf> for Complex {
    type Output = Result<Self>;

    fn div(self, other: Mpf) -> Self::Output {
        if other.is_zero() {
            return Err(crate::error::Error::DomainError("Division by zero".into()));
        }

        Ok(Self::from_real_imag(
            self.real().div(&other)?,
            self.imaginary().div(&other)?,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complex_arithmetic() {
        let z1 = Complex::from_str("3+4i", 64).unwrap();
        let z2 = Complex::from_str("1+2i", 64).unwrap();

        // 测试加法
        let sum = z1.clone().add(&z2).unwrap();
        assert_eq!(sum.real().to_i64(), Some(4));
        assert_eq!(sum.imaginary().to_i64(), Some(6));

        // 测试减法
        let diff = z1.clone().sub(&z2).unwrap();
        assert_eq!(diff.real().to_i64(), Some(2));
        assert_eq!(diff.imaginary().to_i64(), Some(2));

        // 测试乘法
        let product = z1.clone().mul(&z2).unwrap();
        assert_eq!(product.real().to_i64(), Some(-5)); // 3*1 - 4*2 = -5
        assert_eq!(product.imaginary().to_i64(), Some(10)); // 3*2 + 4*1 = 10

        // 测试除法
        let quotient = z1.div(&z2).unwrap();
        // 验证结果
        assert!(!quotient.is_zero());
    }

    #[test]
    fn test_complex_properties() {
        let z = Complex::from_str("3+4i", 64).unwrap();

        // 测试模长
        let mag = z.magnitude().unwrap();
        assert_eq!(mag.to_i64(), Some(5)); // √(3² + 4²) = 5

        // 测试共轭
        let conj = z.conjugate();
        assert_eq!(conj.real().to_i64(), Some(3));
        assert_eq!(conj.imaginary().to_i64(), Some(-4));

        // 测试倒数
        let recip = z.reciprocal().unwrap();
        let product = z.mul(&recip).unwrap();
        assert!(product.is_one() || product.real().to_f64().unwrap().abs() - 1.0 < 1e-10);
    }

    #[test]
    fn test_complex_power() {
        let z = Complex::from_str("1+i", 64).unwrap();

        // 测试平方
        let square = z.square().unwrap();
        assert_eq!(square.real().to_i64(), Some(0)); // (1+i)² = 2i
        assert_eq!(square.imaginary().to_i64(), Some(2));

        // 测试立方
        let cube = z.cube().unwrap();
        assert_eq!(cube.real().to_i64(), Some(-2)); // (1+i)³ = -2 + 2i
        assert_eq!(cube.imaginary().to_i64(), Some(2));
    }
}
