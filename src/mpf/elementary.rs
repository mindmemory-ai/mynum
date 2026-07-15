//! Mpf 初等函数实现

use super::cordic::{Cordic, CordicConfig};
use super::core::Mpf;
use crate::error::{Error, Result};
use crate::mpz::Mpz;

impl Mpf {
    /// 平方根
    pub fn sqrt(&self) -> Result<Mpf> {
        if self.is_zero() {
            return Ok(Mpf::new());
        }

        if self.is_negative() {
            return Err(Error::DomainError(
                "Cannot take square root of negative number".into(),
            ));
        }

        // 使用牛顿法计算平方根
        // Start from max(self, 1) to ensure we start above sqrt,
        // which guarantees monotonic convergence (Newton's method for sqrt
        // decreases monotonically when starting above the root).
        let one = Mpf::from_mpz(Mpz::from_i64(1), self.precision());
        let mut x = if self.cmp(&one) == core::cmp::Ordering::Less {
            one.clone()
        } else {
            self.clone()
        };
        let two = Mpf::from_mpz(Mpz::from_i64(2), self.precision());

        for _ in 0..100 {
            let y = x.add(&self.div(&x)?).div(&two)?;
            // Converged when the value stops decreasing (Newton sqrt monotonic from above)
            if y.cmp(&x) != core::cmp::Ordering::Less {
                break;
            }
            x = y;
        }

        Ok(x)
    }

    /// 立方根
    pub fn cbrt(&self) -> Result<Mpf> {
        if self.is_zero() {
            return Ok(Mpf::new());
        }

        // 使用牛顿法计算立方根
        let mut x = self.clone();
        let three = Mpf::from_mpz(Mpz::from_i64(3), self.precision());

        for _ in 0..100 {
            let x_squared = x.mul(&x);
            let x_cubed = x_squared.mul(&x);

            if x_cubed.cmp(self) == core::cmp::Ordering::Equal {
                break;
            }

            let two_x = x.mul(&Mpf::from_mpz(Mpz::from_i64(2), self.precision()));
            let x_sq_div = self.div(&x_squared)?;
            let numerator = two_x.add(&x_sq_div);
            let new_x = numerator.div(&three)?;

            if new_x.cmp(&x) == core::cmp::Ordering::Equal {
                break;
            }

            x = new_x;
        }

        Ok(x)
    }

    /// 自然对数（高精度优化版）
    pub fn ln(&self) -> Result<Mpf> {
        if self.is_zero() {
            return Err(Error::DomainError("Cannot take logarithm of zero".into()));
        }

        if self.is_negative() {
            return Err(Error::DomainError(
                "Cannot take logarithm of negative number".into(),
            ));
        }

        let precision = self.precision();
        let one = Mpf::from_mpz(Mpz::from_i64(1), precision);

        // 如果 x = 1，返回 0
        if self.cmp(&one) == core::cmp::Ordering::Equal {
            return Ok(Mpf::new());
        }

        // 使用高精度常数ln(2)
        let ln2 = Self::ln_2_constant(precision)?;

        // 范围约简：将x转换到[1/√2, √2]范围内以获得更快收敛
        let mut x = self.clone();
        let mut exponent_adjustment = 0i64;

        let two = Mpf::from_mpz(Mpz::from_i64(2), precision);
        let sqrt2 = two.sqrt()?;
        let inv_sqrt2 = one.div(&sqrt2)?;

        // 约简到最优范围
        while x.cmp(&sqrt2) == core::cmp::Ordering::Greater {
            x = x.div(&two)?;
            exponent_adjustment += 1;
        }
        while x.cmp(&inv_sqrt2) == core::cmp::Ordering::Less {
            x = x.mul(&two);
            exponent_adjustment -= 1;
        }

        // 使用AGM (Arithmetic-Geometric Mean) 算法计算更高精度的对数
        // 对于范围内的值，使用改进的级数展开
        let ln_x = if (x.sub(&one))
            .abs()
            .cmp(&Mpf::from_str("0.1", 10).unwrap_or_else(|_| inv_sqrt2.clone()))
            == core::cmp::Ordering::Less
        {
            // 对于接近1的值，使用快速收敛的级数
            Self::ln_series_optimized(&x, precision)?
        } else {
            // 对于其他值，使用变换后的级数
            Self::ln_transformation(&x, precision)?
        };

        // 加上指数调整：ln(x * 2^n) = ln(x) + n * ln(2)
        let adjustment = ln2.mul(&Mpf::from_mpz(
            Mpz::from_i64(exponent_adjustment),
            precision,
        ));
        Ok(ln_x.add(&adjustment))
    }

    /// 高精度ln(2)常数
    fn ln_2_constant(_precision: usize) -> Result<Mpf> {
        // 使用高精度字符串表示
        let ln2_str = "0.69314718055994530941723212145817656807550013436025525412068000949339362196969471560586332699641868754200148102057068573368552023575813055703267075163507596193072757082837143519030703862389167347112335011536449795523912047517268157493206515552473413952588295045300709532636664265410423915781495204374043038550080194417064167151864471283996817178454695702627163106454615025720740248163777338963855069526066834113727387372292895649354702576265209885969320196505855476470330679365443254763274495125040606943814710468994650622016772042452780726884745251594634876608946";
        Mpf::from_str(ln2_str, 10)
    }

    /// 优化的级数展开（用于接近1的值）
    fn ln_series_optimized(x: &Mpf, precision: usize) -> Result<Mpf> {
        let one = Mpf::from_mpz(Mpz::from_i64(1), precision);
        let t = x.sub(&one);

        // 使用 ln(1+t) = t - t²/2 + t³/3 - t⁴/4 + ...
        // 但使用Horner方法和更高的项数
        let mut result = Mpf::new();
        let mut power = t.clone();

        for n in 1..=100 {
            let n_mpf = Mpf::from_mpz(Mpz::from_i64(n), precision);
            let term = power.div(&n_mpf)?;

            if n % 2 == 1 {
                result = result.add(&term);
            } else {
                result = result.sub(&term);
            }

            power = power.mul(&t);

            // 更严格的收敛检查
            let convergence_threshold = Mpf::from_str("1e-50", 10).unwrap_or_else(|_| {
                Mpf::from_mpz(Mpz::from_i64(1), precision)
                    .div(&Mpf::from_mpz(Mpz::from_i64(1000000000000i64), precision))
                    .unwrap_or_else(|_| Mpf::new())
            });

            if term.abs().cmp(&convergence_threshold) == core::cmp::Ordering::Less {
                break;
            }
        }

        Ok(result)
    }

    /// 使用变换的对数计算
    fn ln_transformation(x: &Mpf, precision: usize) -> Result<Mpf> {
        let one = Mpf::from_mpz(Mpz::from_i64(1), precision);

        // 使用变换 ln(x) = 2 * artanh((x-1)/(x+1))
        // artanh(t) = t + t³/3 + t⁵/5 + ...
        let numerator = x.sub(&one);
        let denominator = x.add(&one);
        let t = numerator.div(&denominator)?;

        let mut result = t.clone();
        let t_squared = t.mul(&t);
        let mut power = t.clone();

        for n in 1..=50 {
            power = power.mul(&t_squared);
            let term = power.div(&Mpf::from_mpz(Mpz::from_i64(2 * n + 1), precision))?;
            result = result.add(&term);

            // 检查收敛
            let convergence_threshold = Mpf::from_str("1e-40", 10).unwrap_or_else(|_| Mpf::new());
            if term.abs().cmp(&convergence_threshold) == core::cmp::Ordering::Less {
                break;
            }
        }

        // 乘以2得到ln(x)
        Ok(result.mul(&Mpf::from_mpz(Mpz::from_i64(2), precision)))
    }

    /// 指数函数
    pub fn exp(&self) -> Result<Mpf> {
        if self.is_zero() {
            return Ok(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));
        }

        let precision = self.precision();
        let one = Mpf::from_mpz(Mpz::from_i64(1), precision);

        // 对于大的负数，使用 e^x = 1/e^(-x)
        if self.is_negative() {
            let pos_x = self.neg();
            let exp_pos = pos_x.exp()?;
            return one.div(&exp_pos);
        }

        // 对于大的正数，使用 e^x = (e^(x/n))^n 来避免溢出
        let mut x = self.clone();
        let mut scaling_factor = 0;
        let max_safe = Mpf::from_mpz(Mpz::from_i64(10), precision);

        while x.cmp(&max_safe) == core::cmp::Ordering::Greater {
            x = x.div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?;
            scaling_factor += 1;
        }

        // 使用级数展开计算指数函数
        // e^x = 1 + x + x²/2! + x³/3! + ...
        let mut result = one.clone();
        let mut term = one.clone();

        for i in 1..=50 {
            term = term
                .mul(&x)
                .div(&Mpf::from_mpz(Mpz::from_i64(i), precision))?;
            result = result.add(&term);

            // 检查收敛
            if term.abs().cmp(
                &Mpf::from_str("1e-15", 10)
                    .unwrap_or_else(|_| Mpf::from_mpz(Mpz::from_i64(1), precision)),
            ) == core::cmp::Ordering::Less
            {
                break;
            }
        }

        // 应用缩放因子 (e^(x/2^n))^(2^n)
        for _ in 0..scaling_factor {
            result = result.mul(&result);
        }

        Ok(result)
    }

    /// 正弦函数（使用CORDIC算法）
    pub fn sin(&self) -> Result<Mpf> {
        if self.is_zero() {
            return Ok(Mpf::new());
        }

        // 对于小角度，使用泰勒展开 sin(x) ≈ x - x³/6 + x⁵/120
        let abs_x = self.abs();
        let small_threshold = Mpf::from_mpz(Mpz::from_i64(1), self.precision())
            .div(&Mpf::from_mpz(Mpz::from_i64(10), self.precision()))?;
        if abs_x.cmp(&small_threshold) == core::cmp::Ordering::Less {
            let x2 = self.mul(self);
            let x3 = x2.mul(self);
            let x5 = x2.mul(&x3);
            let six = Mpf::from_mpz(Mpz::from_i64(6), self.precision());
            let one_twenty = Mpf::from_mpz(Mpz::from_i64(120), self.precision());
            return Ok(self.sub(&x3.div(&six)?).add(&x5.div(&one_twenty)?));
        }

        // 使用CORDIC算法
        let config = CordicConfig::default();
        let cordic = Cordic::new(config)?;
        cordic.sin(self)
    }

    /// 余弦函数（使用CORDIC算法）
    pub fn cos(&self) -> Result<Mpf> {
        if self.is_zero() {
            return Ok(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));
        }

        // 对于小角度，使用 cos(x) ≈ 1 - x²/2
        let abs_x = self.abs();
        let small_threshold = Mpf::from_mpz(Mpz::from_i64(1), self.precision())
            .div(&Mpf::from_mpz(Mpz::from_i64(10), self.precision()))?;
        if abs_x.cmp(&small_threshold) == core::cmp::Ordering::Less {
            let x_squared = self.mul(self);
            let half_x_squared =
                x_squared.div(&Mpf::from_mpz(Mpz::from_i64(2), self.precision()))?;
            return Ok(Mpf::from_mpz(Mpz::from_i64(1), self.precision()).sub(&half_x_squared));
        }

        // 使用CORDIC算法
        let config = CordicConfig::default();
        let cordic = Cordic::new(config)?;
        cordic.cos(self)
    }

    /// 正切函数（使用CORDIC算法）
    pub fn tan(&self) -> Result<Mpf> {
        let sin_x = self.sin()?;
        let cos_x = self.cos()?;

        if cos_x.is_zero() {
            return Err(Error::DomainError("Tangent undefined at this point".into()));
        }

        sin_x.div(&cos_x)
    }

    /// 反正弦函数（使用CORDIC算法）
    pub fn asin(&self) -> Result<Mpf> {
        let abs_x = self.abs();
        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), self.precision()))
            == core::cmp::Ordering::Greater
        {
            return Err(Error::DomainError("Arcsin domain error".into()));
        }

        if self.is_zero() {
            return Ok(Mpf::new());
        }

        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), self.precision()))
            == core::cmp::Ordering::Equal
        {
            // asin(1) = π/2, asin(-1) = -π/2
            let pi_half = Mpf::pi(self.precision())
                .div(&Mpf::from_mpz(Mpz::from_i64(2), self.precision()))?;
            return Ok(if self.is_negative() {
                pi_half.neg()
            } else {
                pi_half
            });
        }

        // 对于小值，使用 asin(x) ≈ x
        let small_threshold = Mpf::from_mpz(Mpz::from_i64(1), self.precision())
            .div(&Mpf::from_mpz(Mpz::from_i64(10), self.precision()))?;
        if abs_x.cmp(&small_threshold) == core::cmp::Ordering::Less {
            return Ok(self.clone());
        }

        // 使用CORDIC算法
        let config = CordicConfig::default();
        let cordic = Cordic::new(config)?;
        cordic.asin(self)
    }

    /// 反余弦函数（使用CORDIC算法）
    pub fn acos(&self) -> Result<Mpf> {
        let abs_x = self.abs();
        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), self.precision()))
            == core::cmp::Ordering::Greater
        {
            return Err(Error::DomainError("Arccos domain error".into()));
        }

        if self.is_zero() {
            // acos(0) = π/2
            let pi = Mpf::pi(self.precision());
            let two = Mpf::from_mpz(Mpz::from_i64(2), self.precision());
            return pi.div(&two);
        }

        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), self.precision()))
            == core::cmp::Ordering::Equal
        {
            // acos(1) = 0, acos(-1) = π
            return Ok(if self.is_negative() {
                Mpf::pi(self.precision())
            } else {
                Mpf::new()
            });
        }

        // 使用CORDIC算法
        let config = CordicConfig::default();
        let cordic = Cordic::new(config)?;
        cordic.acos(self)
    }

    /// 反正切函数（使用CORDIC算法）
    pub fn atan(&self) -> Result<Mpf> {
        if self.is_zero() {
            return Ok(Mpf::new());
        }

        // 对于小值，使用 atan(x) ≈ x
        let abs_x = self.abs();
        let small_threshold = Mpf::from_mpz(Mpz::from_i64(1), self.precision())
            .div(&Mpf::from_mpz(Mpz::from_i64(10), self.precision()))?;
        if abs_x.cmp(&small_threshold) == core::cmp::Ordering::Less {
            return Ok(self.clone());
        }

        // 使用CORDIC算法
        let config = CordicConfig::default();
        let cordic = Cordic::new(config)?;
        cordic.atan(self)
    }

    /// 双曲正弦函数
    pub fn sinh(&self) -> Result<Mpf> {
        if self.is_zero() {
            return Ok(Mpf::new());
        }

        // 使用级数展开：sinh(x) = x + x³/3! + x⁵/5! + ...
        let mut result = self.clone();
        let x_squared = self.mul(self);
        let mut term = self.clone();
        let mut factorial = Mpz::from_i64(1);

        for i in 3..=15 {
            if i % 2 == 1 {
                factorial = factorial.mul(&Mpz::from_i64(i));
                term = term
                    .mul(&x_squared)
                    .div(&Mpf::from_mpz(factorial.clone(), self.precision()))?;
                result = result.add(&term);
            }
        }

        Ok(result)
    }

    /// 双曲余弦函数
    pub fn cosh(&self) -> Result<Mpf> {
        if self.is_zero() {
            return Ok(Mpf::from_mpz(Mpz::from_i64(1), self.precision()));
        }

        // 使用级数展开：cosh(x) = 1 + x²/2! + x⁴/4! + ...
        let mut result = Mpf::from_mpz(Mpz::from_i64(1), self.precision());
        let x_squared = self.mul(self);
        let mut term = Mpf::from_mpz(Mpz::from_i64(1), self.precision());
        let mut factorial = Mpz::from_i64(1);

        for i in 2..=14 {
            if i % 2 == 0 {
                factorial = factorial.mul(&Mpz::from_i64(i));
                term = term
                    .mul(&x_squared)
                    .div(&Mpf::from_mpz(factorial.clone(), self.precision()))?;
                result = result.add(&term);
            }
        }

        Ok(result)
    }

    /// 双曲正切函数
    pub fn tanh(&self) -> Result<Mpf> {
        let sinh_x = self.sinh()?;
        let cosh_x = self.cosh()?;

        if cosh_x.is_zero() {
            return Err(Error::DomainError(
                "Hyperbolic tangent undefined at this point".into(),
            ));
        }

        sinh_x.div(&cosh_x)
    }

    /// 反双曲正弦函数
    pub fn asinh(&self) -> Result<Mpf> {
        if self.is_zero() {
            return Ok(Mpf::new());
        }

        // asinh(x) = ln(x + √(x² + 1))
        let x_squared = self.mul(self);
        let one = Mpf::from_mpz(Mpz::from_i64(1), self.precision());
        let sqrt_term = x_squared.add(&one).sqrt()?;
        let argument = self.add(&sqrt_term);
        argument.ln()
    }

    /// 反双曲余弦函数
    pub fn acosh(&self) -> Result<Mpf> {
        let abs_x = self.abs();
        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), self.precision()))
            == core::cmp::Ordering::Less
        {
            return Err(Error::DomainError("Acosh domain error".into()));
        }

        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), self.precision()))
            == core::cmp::Ordering::Equal
        {
            return Ok(Mpf::new()); // acosh(1) = 0
        }

        // acosh(x) = ln(x + √(x² - 1))
        let x_squared = self.mul(self);
        let one = Mpf::from_mpz(Mpz::from_i64(1), self.precision());
        let sqrt_term = x_squared.sub(&one).sqrt()?;
        let argument = self.add(&sqrt_term);
        argument.ln()
    }

    /// 反双曲正切函数
    pub fn atanh(&self) -> Result<Mpf> {
        let abs_x = self.abs();
        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), self.precision()))
            == core::cmp::Ordering::Greater
        {
            return Err(Error::DomainError("Atanh domain error".into()));
        }

        if self.is_zero() {
            return Ok(Mpf::new());
        }

        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), self.precision()))
            == core::cmp::Ordering::Equal
        {
            return Err(Error::DomainError("Atanh undefined at x = ±1".into()));
        }

        // atanh(x) = (1/2) * ln((1 + x) / (1 - x))
        let one = Mpf::from_mpz(Mpz::from_i64(1), self.precision());
        let numerator = one.add(self);
        let denominator = one.sub(self);
        let ratio = numerator.div(&denominator)?;
        let ln_result = ratio.ln()?;
        let two = Mpf::from_mpz(Mpz::from_i64(2), self.precision());
        ln_result.div(&two)
    }

    /// 使用CORDIC算法计算三角函数（高级接口）
    pub fn sin_cordic(&self, config: CordicConfig) -> Result<Mpf> {
        let cordic = Cordic::new(config)?;
        cordic.sin(self)
    }

    /// 使用CORDIC算法计算余弦函数（高级接口）
    pub fn cos_cordic(&self, config: CordicConfig) -> Result<Mpf> {
        let cordic = Cordic::new(config)?;
        cordic.cos(self)
    }

    /// 使用CORDIC算法计算正切函数（高级接口）
    pub fn tan_cordic(&self, config: CordicConfig) -> Result<Mpf> {
        let cordic = Cordic::new(config)?;
        cordic.tan(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sqrt() {
        let four = Mpf::from_mpz(Mpz::from_i64(4), 64);
        let sqrt_four = four.sqrt().unwrap();
        assert_eq!(sqrt_four.to_i64(), Some(2));
    }

    #[test]
    fn test_cbrt() {
        let one = Mpf::from_mpz(Mpz::from_i64(1), 64);
        let cbrt_one = one.cbrt().unwrap();
        // 立方根1应该等于1
        assert_eq!(cbrt_one.to_i64(), Some(1));
    }

    #[test]
    fn test_exp() {
        let zero = Mpf::new();
        let exp_zero = zero.exp().unwrap();
        assert_eq!(exp_zero.to_i64(), Some(1));
    }
}
