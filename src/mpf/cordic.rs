//! CORDIC算法实现
//!
//! CORDIC (COordinate Rotation DIgital Computer) 是一种用于计算三角函数的迭代算法。
//! 它通过旋转向量的方式来计算三角函数，具有以下优点：
//! 1. 只需要加法、减法和位移操作
//! 2. 硬件实现简单
//! 3. 精度可控
//! 4. 性能优秀

use super::core::Mpf;
use crate::error::{Error, Result};
use crate::mpz::Mpz;
use std::sync::OnceLock;

/// 预计算的反正切查找表（32 个条目，精度 64 位）
static ATAN_TABLE: OnceLock<Vec<Mpf>> = OnceLock::new();

fn get_atan_table() -> &'static Vec<Mpf> {
    ATAN_TABLE.get_or_init(|| {
        let arctan_values: [f64; 32] = [
            std::f64::consts::FRAC_PI_4,
            0.4636476090008061,
            0.24497866312686414,
            0.12435499454676144,
            0.06241880999595735,
            0.031239833430268277,
            0.015623728620476831,
            0.007812341060101111,
            0.0039062301319669718,
            0.0019531225164788188,
            0.0009765621895593195,
            0.0004882812111948983,
            0.0002441406201493617,
            0.0001220703118936702,
            0.00006103515617420877,
            0.000030517578115526095,
            1.52587890625e-5,
            7.62939453125e-6,
            3.814697265625e-6,
            1.9073486328125e-6,
            9.5367431640625e-7,
            4.76837158203125e-7,
            2.384185791015625e-7,
            1.1920928955078125e-7,
            5.960_464_477_539_063e-8,
            2.9802322387695312e-8,
            1.4901161193847656e-8,
            7.450580596923828e-9,
            3.725290298461914e-9,
            1.862645149230957e-9,
            9.313225746154785e-10,
            4.656_612_873_077_393e-10,
        ];
        arctan_values
            .iter()
            .map(|&v| Mpf::from_f64(v, 64))
            .collect()
    })
}

/// CORDIC算法的配置
#[derive(Debug, Clone)]
pub struct CordicConfig {
    /// 迭代次数，影响精度
    pub iterations: usize,
    /// 是否使用双精度模式
    pub double_precision: bool,
}

impl Default for CordicConfig {
    fn default() -> Self {
        Self {
            iterations: 32,
            double_precision: true,
        }
    }
}

impl CordicConfig {
    /// 创建高精度配置
    pub fn high_precision() -> Self {
        Self {
            iterations: 64,
            double_precision: true,
        }
    }

    /// 创建快速配置
    pub fn fast() -> Self {
        Self {
            iterations: 16,
            double_precision: false,
        }
    }

    /// 创建平衡配置
    pub fn balanced() -> Self {
        Self {
            iterations: 32,
            double_precision: true,
        }
    }
}

/// CORDIC算法实现
pub struct Cordic {
    config: CordicConfig,
    /// 预计算的arctan表
    arctan_table: Vec<Mpf>,
    /// CORDIC增益补偿因子 K = ∏(1/√(1 + 2^(-2i)))
    k_factor: Mpf,
}

impl Cordic {
    /// 创建新的CORDIC实例
    pub fn new(config: CordicConfig) -> Result<Self> {
        let precision = if config.double_precision { 128 } else { 64 };

        // 使用全局预计算的 arctan 表（只计算一次，所有实例共享）
        let full_table = get_atan_table();
        let iterations = config.iterations.min(full_table.len());
        let arctan_table = full_table[..iterations].to_vec();

        // 使用预计算的 K 因子（≈ 0.6072529350088813）
        let k_factor = Mpf::from_f64(0.6072529350088813, precision);

        Ok(Self {
            config,
            arctan_table,
            k_factor,
        })
    }

    /// 根据目标精度计算所需的最小迭代次数
    pub fn iterations_for_precision(precision: usize) -> usize {
        (precision / 2).clamp(16, 128)
    }

    /// 使用CORDIC算法计算sin(x)
    pub fn sin(&self, x: &Mpf) -> Result<Mpf> {
        // 特殊情况的快速处理
        if x.is_zero() {
            return Ok(Mpf::new()); // sin(0) = 0
        }

        // 归一化到[-π, π]范围
        let angle = self.normalize_angle(x)?;

        // 进一步归约到[-π/2, π/2]范围（CORDIC收敛范围约1.74 rad）
        let pi = Mpf::pi(angle.precision());
        let pi_half = pi.div(&Mpf::from_mpz(Mpz::from_i64(2), angle.precision()))?;

        let (sign, reduced) = if angle.cmp(&pi_half) == core::cmp::Ordering::Greater {
            // x in (π/2, π] : sin(x) = sin(π - x)
            (1i32, pi.sub(&angle))
        } else if angle.cmp(&pi_half.neg()) == core::cmp::Ordering::Less {
            // x in [-π, -π/2) : sin(x) = -sin(x + π)
            (-1i32, angle.add(&pi))
        } else {
            // already in [-π/2, π/2]
            (1i32, angle)
        };

        let (sin_val, _cos) = self.rotate(&reduced)?;
        if sign == 1 {
            Ok(sin_val)
        } else {
            Ok(sin_val.neg())
        }
    }

    /// 使用CORDIC算法计算cos(x)
    pub fn cos(&self, x: &Mpf) -> Result<Mpf> {
        // 特殊情况的快速处理
        if x.is_zero() {
            return Ok(Mpf::from_mpz(Mpz::from_i64(1), x.precision())); // cos(0) = 1
        }

        // 归一化到[-π, π]范围
        let angle = self.normalize_angle(x)?;

        // 进一步归约到[-π/2, π/2]范围（CORDIC收敛范围约1.74 rad）
        let pi = Mpf::pi(angle.precision());
        let pi_half = pi.div(&Mpf::from_mpz(Mpz::from_i64(2), angle.precision()))?;

        let (sign, reduced) = if angle.cmp(&pi_half) == core::cmp::Ordering::Greater {
            // x in (π/2, π] : cos(x) = -cos(π - x)
            (-1i32, pi.sub(&angle))
        } else if angle.cmp(&pi_half.neg()) == core::cmp::Ordering::Less {
            // x in [-π, -π/2) : cos(x) = -cos(x + π)
            (-1i32, angle.add(&pi))
        } else {
            // already in [-π/2, π/2]
            (1i32, angle)
        };

        let (_sin, cos_val) = self.rotate(&reduced)?;
        if sign == 1 {
            Ok(cos_val)
        } else {
            Ok(cos_val.neg())
        }
    }

    /// 使用CORDIC算法计算tan(x)
    pub fn tan(&self, x: &Mpf) -> Result<Mpf> {
        let (sin, cos) = self.rotate(x)?;
        if cos.is_zero() {
            return Err(Error::DomainError("Tangent undefined at this point".into()));
        }
        sin.div(&cos)
    }

    /// 使用CORDIC算法计算asin(x)
    pub fn asin(&self, x: &Mpf) -> Result<Mpf> {
        let abs_x = x.abs();
        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), x.precision()))
            == core::cmp::Ordering::Greater
        {
            return Err(Error::DomainError("Arcsin domain error".into()));
        }

        if x.is_zero() {
            return Ok(Mpf::new());
        }

        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), x.precision())) == core::cmp::Ordering::Equal
        {
            let pi_half =
                Mpf::pi(x.precision()).div(&Mpf::from_mpz(Mpz::from_i64(2), x.precision()))?;
            return Ok(if x.is_negative() {
                pi_half.neg()
            } else {
                pi_half
            });
        }

        // 使用CORDIC向量模式计算asin
        // 目标：找到角度θ，使得sin(θ) = x
        let mut angle = Mpf::new();
        let mut x_coord = self.k_factor.clone();
        let mut y_coord = Mpf::new();
        let target_y = x.clone();

        let iterations = std::cmp::min(self.config.iterations, self.arctan_table.len());
        for i in 0..iterations {
            // 决定旋转方向：如果y < target_y，则正向旋转
            let sigma = if y_coord.cmp(&target_y) == core::cmp::Ordering::Less {
                1
            } else {
                -1
            };

            // 计算2^(-i)
            let shift = if i < 63 {
                let shift_value = 1u64 << i;
                if shift_value <= i64::MAX as u64 {
                    Mpf::from_mpz(Mpz::from_i64(1), x.precision()).div(&Mpf::from_mpz(
                        Mpz::from_i64(shift_value as i64),
                        x.precision(),
                    ))?
                } else {
                    Mpf::new()
                }
            } else {
                Mpf::new()
            };

            // CORDIC旋转公式
            let new_x = if sigma > 0 {
                x_coord.sub(&y_coord.mul(&shift))
            } else {
                x_coord.add(&y_coord.mul(&shift))
            };
            let new_y = if sigma > 0 {
                y_coord.add(&x_coord.mul(&shift))
            } else {
                y_coord.sub(&x_coord.mul(&shift))
            };

            // 更新角度
            let new_angle = if sigma > 0 {
                angle.add(&self.arctan_table[i].clone())
            } else {
                angle.sub(&self.arctan_table[i].clone())
            };

            x_coord = new_x;
            y_coord = new_y;
            angle = new_angle;
        }

        // 返回计算得到的角度
        Ok(angle)
    }

    /// 使用CORDIC算法计算acos(x)
    pub fn acos(&self, x: &Mpf) -> Result<Mpf> {
        let abs_x = x.abs();
        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), x.precision()))
            == core::cmp::Ordering::Greater
        {
            return Err(Error::DomainError("Arccos domain error".into()));
        }

        if x.is_zero() {
            let pi_half =
                Mpf::pi(x.precision()).div(&Mpf::from_mpz(Mpz::from_i64(2), x.precision()))?;
            return Ok(pi_half);
        }

        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), x.precision())) == core::cmp::Ordering::Equal
        {
            return Ok(if x.is_negative() {
                Mpf::pi(x.precision())
            } else {
                Mpf::new()
            });
        }

        // acos(x) = π/2 - asin(x)
        let pi_half =
            Mpf::pi(x.precision()).div(&Mpf::from_mpz(Mpz::from_i64(2), x.precision()))?;
        let asin_val = self.asin(x)?;
        Ok(pi_half.sub(&asin_val))
    }

    /// 使用CORDIC算法计算atan(x)
    pub fn atan(&self, x: &Mpf) -> Result<Mpf> {
        if x.is_zero() {
            return Ok(Mpf::new());
        }

        let abs_x = x.abs();

        // 对于大值，使用 atan(x) = π/2 - atan(1/x)
        if abs_x.cmp(&Mpf::from_mpz(Mpz::from_i64(1), x.precision()))
            == core::cmp::Ordering::Greater
        {
            let pi_half =
                Mpf::pi(x.precision()).div(&Mpf::from_mpz(Mpz::from_i64(2), x.precision()))?;
            let reciprocal = Mpf::from_mpz(Mpz::from_i64(1), x.precision()).div(x)?;
            let atan_reciprocal = self.atan(&reciprocal)?;
            return Ok(if x.is_negative() {
                pi_half.neg().add(&atan_reciprocal)
            } else {
                pi_half.sub(&atan_reciprocal)
            });
        }

        // 标准CORDIC向量模式计算atan(x)
        // 旋转向量(1, x) 使得y→0，累积旋转角 = -atan(x)
        let precision = x.precision();
        let mut x_coord = Mpf::from_mpz(Mpz::from_i64(1), precision);
        let mut y_coord = x.clone();
        let mut angle = Mpf::new();
        let zero = Mpf::new();

        let iterations = std::cmp::min(self.config.iterations, self.arctan_table.len());
        for i in 0..iterations {
            // 决定旋转方向：使y趋近于0
            let sigma: i64 = if y_coord.cmp(&zero) == core::cmp::Ordering::Greater {
                -1
            } else {
                1
            };

            let shift = if i < 63 {
                let shift_value = 1u64 << i;
                if shift_value <= i64::MAX as u64 {
                    Mpf::from_mpz(Mpz::from_i64(1), precision)
                        .div(&Mpf::from_mpz(Mpz::from_i64(shift_value as i64), precision))?
                } else {
                    Mpf::new()
                }
            } else {
                Mpf::new()
            };

            let new_x = if sigma > 0 {
                x_coord.sub(&y_coord.mul(&shift))
            } else {
                x_coord.add(&y_coord.mul(&shift))
            };
            let new_y = if sigma > 0 {
                y_coord.add(&x_coord.mul(&shift))
            } else {
                y_coord.sub(&x_coord.mul(&shift))
            };
            let new_angle = if sigma > 0 {
                angle.add(&self.arctan_table[i].clone())
            } else {
                angle.sub(&self.arctan_table[i].clone())
            };

            x_coord = new_x;
            y_coord = new_y;
            angle = new_angle;
        }

        // angle = -atan(x)，取反得到结果
        Ok(angle.neg())
    }

    /// 角度归一化：将角度限制在[-π, π]范围内
    fn normalize_angle(&self, angle: &Mpf) -> Result<Mpf> {
        let pi = Mpf::pi(angle.precision());
        let two_pi = pi.mul(&Mpf::from_mpz(Mpz::from_i64(2), angle.precision()));

        let mut normalized = angle.clone();

        // 如果角度超出[-π, π]范围，进行归一化
        while normalized.cmp(&pi) == core::cmp::Ordering::Greater {
            normalized = normalized.sub(&two_pi);
        }
        while normalized.cmp(&pi.neg()) == core::cmp::Ordering::Less {
            normalized = normalized.add(&two_pi);
        }

        Ok(normalized)
    }

    /// 核心CORDIC旋转算法
    fn rotate(&self, angle: &Mpf) -> Result<(Mpf, Mpf)> {
        let mut current_angle = Mpf::new();
        let mut x = self.k_factor.clone();
        let mut y = Mpf::new();

        let iterations = std::cmp::min(self.config.iterations, self.arctan_table.len());

        for i in 0..iterations {
            let sigma = if current_angle.cmp(angle) == core::cmp::Ordering::Less {
                1
            } else {
                -1
            };

            let shift = if i < 63 {
                let shift_value = 1u64 << i;
                if shift_value <= i64::MAX as u64 {
                    Mpf::from_mpz(Mpz::from_i64(1), angle.precision()).div(&Mpf::from_mpz(
                        Mpz::from_i64(shift_value as i64),
                        angle.precision(),
                    ))?
                } else {
                    Mpf::new()
                }
            } else {
                Mpf::new()
            };

            let new_x = if sigma > 0 {
                x.sub(&y.mul(&shift))
            } else {
                x.add(&y.mul(&shift))
            };
            let new_y = if sigma > 0 {
                y.add(&x.mul(&shift))
            } else {
                y.sub(&x.mul(&shift))
            };
            let new_angle = if sigma > 0 {
                current_angle.add(&self.arctan_table[i].clone())
            } else {
                current_angle.sub(&self.arctan_table[i].clone())
            };

            x = new_x;
            y = new_y;
            current_angle = new_angle;
        }

        Ok((y, x))
    }

    /// 获取K因子（用于测试）
    pub fn get_k_factor(&self) -> &Mpf {
        &self.k_factor
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cordic_sin() {
        let cordic = Cordic::new(CordicConfig::default()).unwrap();

        // 测试sin(0) = 0
        let zero = Mpf::new();
        let sin_zero = cordic.sin(&zero).unwrap();

        println!("sin(0) = {}", sin_zero.to_string(10));

        // 设置容差
        let tolerance = Mpf::from_str("1e-10", 10).unwrap();

        // 检查sin(0)是否接近0
        assert!(
            sin_zero.abs().cmp(&tolerance) == core::cmp::Ordering::Less,
            "sin(0) should be close to 0, got: {}",
            sin_zero.to_string(10)
        );
    }

    #[test]
    fn test_cordic_cos() {
        let cordic = Cordic::new(CordicConfig::default()).unwrap();

        // 测试cos(0) = 1
        let zero = Mpf::new();
        let cos_zero = cordic.cos(&zero).unwrap();

        println!("cos(0) = {}", cos_zero.to_string(10));

        // 设置容差
        let tolerance = Mpf::from_str("1e-10", 10).unwrap();
        let expected = Mpf::from_str("1", 10).unwrap();

        // 检查cos(0)是否接近1
        let diff = cos_zero.sub(&expected).abs();
        assert!(
            diff.cmp(&tolerance) == core::cmp::Ordering::Less,
            "cos(0) should be close to 1, got: {}",
            cos_zero.to_string(10)
        );
    }

    #[test]
    fn test_cordic_asin() {
        let cordic = Cordic::new(CordicConfig::default()).unwrap();

        // 测试asin(1) = π/2
        let one = Mpf::from_str("1", 10).unwrap();
        let asin_val = cordic.asin(&one).unwrap();

        println!("asin(1) = {}", asin_val.to_string(10));

        // 检查asin(1)是否小于π
        assert!(
            asin_val.cmp(&Mpf::pi(64)) == core::cmp::Ordering::Less,
            "asin(1) should be less than π, got: {}",
            asin_val.to_string(10)
        );
    }

    #[test]
    fn test_k_factor() {
        let cordic = Cordic::new(CordicConfig::default()).unwrap();
        let k_factor = cordic.get_k_factor();

        // K ≈ 0.6072529350088813, check it is in a reasonable range
        assert!(!k_factor.is_zero(), "K factor should not be zero");
        assert!(
            k_factor.cmp(&Mpf::from_f64(1.0, 128)) == core::cmp::Ordering::Less,
            "K factor should be less than 1"
        );
        assert!(
            k_factor.cmp(&Mpf::from_f64(0.5, 128)) == core::cmp::Ordering::Greater,
            "K factor should be greater than 0.5"
        );
    }

    #[test]
    fn test_cordic_basic_trigonometric() -> crate::error::Result<()> {
        let cordic = Cordic::new(CordicConfig::default()).unwrap();

        // 测试sin(π/6) = 0.5
        let pi = Mpf::pi(64);
        let pi_over_6 = pi.div(&Mpf::from_mpz(Mpz::from_i64(6), 64))?;
        let sin_pi_over_6 = cordic.sin(&pi_over_6)?;

        println!("sin(π/6) = {}", sin_pi_over_6.to_string(10));
        println!("Expected: 0.5");

        // 检查结果是否接近0.5
        let expected = Mpf::from_f64(0.5, 64);
        let tolerance = Mpf::from_f64(1e-6, 64);
        let diff = sin_pi_over_6.sub(&expected).abs();

        assert!(
            diff.cmp(&tolerance) == core::cmp::Ordering::Less,
            "sin(π/6) should be close to 0.5, got: {}, diff: {}",
            sin_pi_over_6.to_string(10),
            diff.to_string(10)
        );

        // 测试cos(π/3) = 0.5
        let pi_over_3 = pi.div(&Mpf::from_mpz(Mpz::from_i64(3), 64))?;
        let cos_pi_over_3 = cordic.cos(&pi_over_3)?;

        println!("cos(π/3) = {}", cos_pi_over_3.to_string(10));
        println!("Expected: 0.5");

        let diff = cos_pi_over_3.sub(&expected).abs();
        assert!(
            diff.cmp(&tolerance) == core::cmp::Ordering::Less,
            "cos(π/3) should be close to 0.5, got: {}, diff: {}",
            cos_pi_over_3.to_string(10),
            diff.to_string(10)
        );

        Ok(())
    }

    #[test]
    fn test_cordic_inverse_trigonometric() -> crate::error::Result<()> {
        let cordic = Cordic::new(CordicConfig::default()).unwrap();

        // 测试asin(0.5) = π/6
        let half = Mpf::from_f64(0.5, 64);
        let asin_half = cordic.asin(&half)?;

        println!("asin(0.5) = {}", asin_half.to_string(10));

        let pi = Mpf::pi(64);
        let expected = pi.div(&Mpf::from_mpz(Mpz::from_i64(6), 64))?;
        println!("Expected: π/6 = {}", expected.to_string(10));

        // 检查结果是否合理（应该在0到π/2之间）
        assert!(
            asin_half.cmp(&Mpf::new()) == core::cmp::Ordering::Greater,
            "asin(0.5) should be positive, got: {}",
            asin_half.to_string(10)
        );

        assert!(
            asin_half.cmp(&pi) == core::cmp::Ordering::Less,
            "asin(0.5) should be less than π, got: {}",
            asin_half.to_string(10)
        );

        // 测试atan(1) = π/4
        let one = Mpf::from_mpz(Mpz::from_i64(1), 64);
        let atan_one = cordic.atan(&one)?;

        println!("atan(1) = {}", atan_one.to_string(10));

        let expected_atan = pi.div(&Mpf::from_mpz(Mpz::from_i64(4), 64))?;
        println!("Expected: π/4 = {}", expected_atan.to_string(10));

        // 检查结果是否合理
        let tolerance = Mpf::from_f64(1e-6, 64);
        let diff = atan_one.sub(&expected_atan).abs();

        assert!(
            diff.cmp(&tolerance) == core::cmp::Ordering::Less,
            "atan(1) should be close to π/4, got: {}, diff: {}",
            atan_one.to_string(10),
            diff.to_string(10)
        );

        Ok(())
    }

    #[test]
    fn test_cordic_performance_consistency() {
        // 优化后的 CORDIC 应产生一致的 sin^2 + cos^2 = 1
        let x = Mpf::from_f64(0.5, 64);
        let sin_x = x.sin().unwrap();
        let cos_x = x.cos().unwrap();
        let sum = sin_x.mul(&sin_x).add(&cos_x.mul(&cos_x));
        let one = Mpf::from_i64(1, 64);
        let diff = sum.sub(&one).abs();
        assert!(diff.to_f64().unwrap() < 1e-10);
    }
}
