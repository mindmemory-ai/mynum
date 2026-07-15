//! 复数数值积分模块
//!
//! 包含复函数的数值积分算法，如梯形法则、辛普森法则、高斯求积等

use crate::complex::core::Complex;
use crate::error::{Error, Result};
use crate::mpf::Mpf;
use std::sync::Arc;

/// 复函数类型定义
pub type ComplexFunction = Arc<dyn Fn(&Complex) -> Result<Complex> + Send + Sync>;

/// 数值积分器结构
#[derive(Debug, Clone)]
pub struct ComplexIntegrator {
    /// 默认精度
    precision: usize,
    /// 默认积分区间分割数
    default_n: usize,
}

impl ComplexIntegrator {
    /// 创建新的积分器
    pub fn new(precision: usize) -> Self {
        Self {
            precision,
            default_n: 1000,
        }
    }

    /// 设置默认分割数
    pub fn with_default_n(mut self, n: usize) -> Self {
        self.default_n = n;
        self
    }

    /// 梯形法则数值积分
    ///
    /// 使用梯形法则计算复函数在指定区间上的定积分
    ///
    /// # 参数
    /// - `f`: 被积函数
    /// - `a`: 积分下限
    /// - `b`: 积分上限
    /// - `n`: 分割数（可选）
    ///
    /// # 返回
    /// - `Result<Complex>`: 积分结果
    ///
    /// # 示例
    /// ```
    /// use mynum::complex::integration::{ComplexIntegrator, ComplexFunction};
    /// use mynum::complex::Complex;
    /// use mynum::mpf::Mpf;
    /// use std::sync::Arc;
    ///
    /// let integrator = ComplexIntegrator::new(64);
    /// let f: ComplexFunction = Arc::new(|z| Ok(z.clone()));
    ///
    /// let a = Complex::from_real(Mpf::from_i64(0, 64));
    /// let b = Complex::from_real(Mpf::from_i64(1, 64));
    /// let result = integrator.trapezoidal(&f, &a, &b, Some(100)).unwrap();
    /// ```
    pub fn trapezoidal(
        &self,
        f: &ComplexFunction,
        a: &Complex,
        b: &Complex,
        n: Option<usize>,
    ) -> Result<Complex> {
        let n = n.unwrap_or(self.default_n);
        if n == 0 {
            return Err(Error::InvalidInput(
                "Number of divisions must be positive".into(),
            ));
        }

        let h = b
            .sub(a)?
            .div(&Complex::from_real(Mpf::from_i64(n as i64, self.precision)))?;
        let mut sum = Complex::new();

        // 计算函数值并累加
        for i in 0..=n {
            let x = if i == 0 {
                a.clone()
            } else if i == n {
                b.clone()
            } else {
                a.add(&h.mul(&Complex::from_real(Mpf::from_i64(i as i64, self.precision)))?)?
            };

            let fx = f(&x)?;
            let weight = if i == 0 || i == n { 1.0 } else { 2.0 };
            let weighted_fx = fx.mul(&Complex::from_real(Mpf::from_f64(weight, self.precision)))?;
            sum = sum.add(&weighted_fx)?;
        }

        // 乘以步长的一半
        let half_h = h.mul(&Complex::from_real(Mpf::from_f64(0.5, self.precision)))?;
        sum.mul(&half_h)
    }

    /// 辛普森法则数值积分
    ///
    /// 使用辛普森法则计算复函数在指定区间上的定积分
    /// 辛普森法则比梯形法则有更高的精度
    ///
    /// # 参数
    /// - `f`: 被积函数
    /// - `a`: 积分下限
    /// - `b`: 积分上限
    /// - `n`: 分割数（必须是偶数）
    ///
    /// # 返回
    /// - `Result<Complex>`: 积分结果
    ///
    /// # 注意
    /// - 分割数n必须是偶数，如果不是偶数会自动调整为偶数
    pub fn simpson(
        &self,
        f: &ComplexFunction,
        a: &Complex,
        b: &Complex,
        n: Option<usize>,
    ) -> Result<Complex> {
        let mut n = n.unwrap_or(self.default_n);
        if n == 0 {
            return Err(Error::InvalidInput(
                "Number of divisions must be positive".into(),
            ));
        }

        // 确保n是偶数
        if !n.is_multiple_of(2) {
            n += 1;
        }

        let h = b
            .sub(a)?
            .div(&Complex::from_real(Mpf::from_i64(n as i64, self.precision)))?;
        let mut sum = Complex::new();

        // 计算函数值并累加
        for i in 0..=n {
            let x = if i == 0 {
                a.clone()
            } else if i == n {
                b.clone()
            } else {
                a.add(&h.mul(&Complex::from_real(Mpf::from_i64(i as i64, self.precision)))?)?
            };

            let fx = f(&x)?;
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let weighted_fx = fx.mul(&Complex::from_real(Mpf::from_f64(weight, self.precision)))?;
            sum = sum.add(&weighted_fx)?;
        }

        // 乘以步长的三分之一
        let third_h = h.mul(&Complex::from_real(Mpf::from_f64(
            1.0 / 3.0,
            self.precision,
        )))?;
        sum.mul(&third_h)
    }

    /// 高斯求积法
    ///
    /// 使用高斯-勒让德求积公式计算复函数在指定区间上的定积分
    /// 高斯求积法在较少的节点上就能达到很高的精度
    ///
    /// # 参数
    /// - `f`: 被积函数
    /// - `a`: 积分下限
    /// - `b`: 积分上限
    /// - `n`: 高斯点数量（可选，默认为5）
    ///
    /// # 返回
    /// - `Result<Complex>`: 积分结果
    pub fn gauss_legendre(
        &self,
        f: &ComplexFunction,
        a: &Complex,
        b: &Complex,
        n: Option<usize>,
    ) -> Result<Complex> {
        let n = n.unwrap_or(5);
        if n == 0 {
            return Err(Error::InvalidInput(
                "Number of Gauss points must be positive".into(),
            ));
        }

        // 获取高斯点和权重（这里使用预定义的值）
        let (points, weights) = self.get_gauss_points_and_weights(n)?;

        let mut sum = Complex::new();
        let half_interval = b
            .sub(a)?
            .mul(&Complex::from_real(Mpf::from_f64(0.5, self.precision)))?;
        let mid_point = a
            .add(b)?
            .mul(&Complex::from_real(Mpf::from_f64(0.5, self.precision)))?;

        // 计算加权和
        for (point, weight) in points.iter().zip(weights.iter()) {
            let point_complex = Complex::from_real(Mpf::from_f64(*point, self.precision));
            let x = mid_point.add(&half_interval.mul(&point_complex)?)?;
            let fx = f(&x)?;
            let weight_complex = Complex::from_real(Mpf::from_f64(*weight, self.precision));
            let weighted_fx = fx.mul(&weight_complex)?;
            sum = sum.add(&weighted_fx)?;
        }

        // 乘以区间长度的一半
        sum.mul(&half_interval)
    }

    /// 自适应积分
    ///
    /// 使用自适应算法自动调整分割数以达到指定精度
    ///
    /// # 参数
    /// - `f`: 被积函数
    /// - `a`: 积分下限
    /// - `b`: 积分上限
    /// - `tolerance`: 容差
    /// - `max_iterations`: 最大迭代次数
    ///
    /// # 返回
    /// - `Result<Complex>`: 积分结果
    pub fn adaptive(
        &self,
        f: &ComplexFunction,
        a: &Complex,
        b: &Complex,
        tolerance: f64,
        max_iterations: Option<usize>,
    ) -> Result<Complex> {
        let max_iterations = max_iterations.unwrap_or(10);
        let mut n = 10;
        let mut prev_result = self.simpson(f, a, b, Some(n))?;

        for _ in 0..max_iterations {
            n *= 2;
            let current_result = self.simpson(f, a, b, Some(n))?;

            // 计算相对误差
            let diff = current_result.sub(&prev_result)?;
            let (diff_real, diff_imag) = diff.to_f64().unwrap_or((0.0, 0.0));
            let (current_real, current_imag) = current_result.to_f64().unwrap_or((0.0, 0.0));

            let relative_error = (diff_real * diff_real + diff_imag * diff_imag).sqrt()
                / (current_real * current_real + current_imag * current_imag).sqrt();

            if relative_error < tolerance {
                return Ok(current_result);
            }

            prev_result = current_result;
        }

        Ok(prev_result)
    }

    /// 获取高斯点和权重
    ///
    /// 返回预定义的高斯-勒让德求积点和权重
    ///
    /// # 参数
    /// - `n`: 高斯点数量
    ///
    /// # 返回
    /// - `Result<(Vec<f64>, Vec<f64>)>`: (高斯点, 权重)
    fn get_gauss_points_and_weights(&self, n: usize) -> Result<(Vec<f64>, Vec<f64>)> {
        match n {
            1 => Ok((vec![0.0], vec![2.0])),
            2 => Ok((
                vec![-0.5773502691896257, 0.5773502691896257],
                vec![1.0, 1.0],
            )),
            3 => Ok((
                vec![-0.7745966692414834, 0.0, 0.7745966692414834],
                vec![0.5555555555555556, 0.8888888888888888, 0.5555555555555556],
            )),
            4 => Ok((
                vec![
                    -0.8611363115940526,
                    -0.3399810435848563,
                    0.3399810435848563,
                    0.8611363115940526,
                ],
                vec![
                    0.3478548451374538,
                    0.6521451548625461,
                    0.6521451548625461,
                    0.3478548451374538,
                ],
            )),
            5 => Ok((
                vec![
                    -0.906_179_845_938_664,
                    -0.5384693101056831,
                    0.0,
                    0.5384693101056831,
                    0.906_179_845_938_664,
                ],
                vec![
                    0.2369268850561891,
                    0.4786286704993665,
                    0.5688888888888889,
                    0.4786286704993665,
                    0.2369268850561891,
                ],
            )),
            _ => {
                // 对于更大的n，返回一个近似值
                let mut points = Vec::with_capacity(n);
                let mut weights = Vec::with_capacity(n);

                for i in 0..n {
                    let x = -1.0 + 2.0 * (i as f64 + 0.5) / (n as f64);
                    let w = 2.0 / (n as f64);
                    points.push(x);
                    weights.push(w);
                }

                Ok((points, weights))
            }
        }
    }

    /// 复平面上的路径积分
    ///
    /// 计算复函数沿指定路径的积分
    ///
    /// # 参数
    /// - `f`: 被积函数
    /// - `path`: 路径函数，参数t从0到1
    /// - `n`: 分割数
    ///
    /// # 返回
    /// - `Result<Complex>`: 路径积分结果
    pub fn path_integral(
        &self,
        f: &ComplexFunction,
        path: &ComplexFunction,
        n: Option<usize>,
    ) -> Result<Complex> {
        let n = n.unwrap_or(self.default_n);
        if n == 0 {
            return Err(Error::InvalidInput(
                "Number of divisions must be positive".into(),
            ));
        }

        let mut sum = Complex::new();

        for i in 0..n {
            let t = Mpf::from_f64((i as f64) / (n as f64), self.precision);
            let t_complex = Complex::from_real(t);

            // 计算路径上的点
            let z = path(&t_complex)?;

            // 计算函数值
            let fz = f(&z)?;

            // 计算路径导数（数值微分）
            let t_next = Mpf::from_f64((i as f64 + 1.0) / (n as f64), self.precision);
            let t_next_complex = Complex::from_real(t_next);
            let z_next = path(&t_next_complex)?;
            let dz = z_next.sub(&z)?;

            // 累加 f(z) * dz
            let integrand = fz.mul(&dz)?;
            sum = sum.add(&integrand)?;
        }

        Ok(sum)
    }
}

/// 预定义的积分路径
pub mod paths {
    use super::*;

    /// 单位圆路径：z(t) = e^(2πit), t ∈ [0, 1]
    pub fn unit_circle(t: &Complex) -> Result<Complex> {
        let two_pi = Mpf::from_f64(2.0 * std::f64::consts::PI, t.precision());
        // z(t) = exp(i * 2π * t) = cos(2πt) + i*sin(2πt)
        let angle_imag = t.real().mul(&two_pi);
        let angle = Complex::from_real_imag(Mpf::with_precision(t.precision()), angle_imag);
        angle.exp()
    }

    /// 直线路径：z(t) = a + (b-a)t, t ∈ [0, 1]
    pub fn line_segment(a: Complex, b: Complex) -> ComplexFunction {
        Arc::new(move |t| {
            let one_minus_t = Complex::from_real(Mpf::from_f64(1.0, t.precision())).sub(t)?;
            let term1 = a.mul(&one_minus_t)?;
            let term2 = b.mul(t)?;
            term1.add(&term2)
        })
    }

    /// 矩形路径：z(t) = 四个边界的参数化
    pub fn rectangle(center: Complex, width: Mpf, height: Mpf) -> ComplexFunction {
        Arc::new(move |t| {
            let t_val = t.real().to_f64().unwrap_or(0.0);
            let half_width = width.mul(&Mpf::from_f64(0.5, width.precision()));
            let half_height = height.mul(&Mpf::from_f64(0.5, height.precision()));

            if t_val < 0.25 {
                // 下边：从 (center.x - width/2, center.y - height/2) 到 (center.x + width/2, center.y - height/2)
                let s = t_val * 4.0;
                let s_mpf = Mpf::from_f64(2.0 * s, center.precision());
                let x = center.real().sub(&half_width).add(&half_width.mul(&s_mpf));
                let y = center.imaginary().sub(&half_height);
                Ok(Complex::from_real_imag(x, y))
            } else if t_val < 0.5 {
                // 右边：从 (center.x + width/2, center.y - height/2) 到 (center.x + width/2, center.y + height/2)
                let s = (t_val - 0.25) * 4.0;
                let s_mpf = Mpf::from_f64(2.0 * s, center.precision());
                let x = center.real().add(&half_width);
                let y = center
                    .imaginary()
                    .sub(&half_height)
                    .add(&half_height.mul(&s_mpf));
                Ok(Complex::from_real_imag(x, y))
            } else if t_val < 0.75 {
                // 上边：从 (center.x + width/2, center.y + height/2) 到 (center.x - width/2, center.y + height/2)
                let s = (t_val - 0.5) * 4.0;
                let s_mpf = Mpf::from_f64(2.0 * s, center.precision());
                let x = center.real().add(&half_width).sub(&half_width.mul(&s_mpf));
                let y = center.imaginary().add(&half_height);
                Ok(Complex::from_real_imag(x, y))
            } else {
                // 左边：从 (center.x - width/2, center.y + height/2) 到 (center.x - width/2, center.y - height/2)
                let s = (t_val - 0.75) * 4.0;
                let s_mpf = Mpf::from_f64(2.0 * s, center.precision());
                let x = center.real().sub(&half_width);
                let y = center
                    .imaginary()
                    .add(&half_height)
                    .sub(&half_height.mul(&s_mpf));
                Ok(Complex::from_real_imag(x, y))
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trapezoidal_integration() {
        let integrator = ComplexIntegrator::new(64);
        let f: ComplexFunction = Arc::new(|z| Ok(z.clone()));

        let a = Complex::from_real(Mpf::from_i64(0, 64));
        let b = Complex::from_real(Mpf::from_i64(1, 64));
        let result = integrator.trapezoidal(&f, &a, &b, Some(100)).unwrap();

        // 积分 z dz 从 0 到 1 = 0.5
        let expected = Complex::from_real(Mpf::from_f64(0.5, 64));
        let (result_real, result_imag) = result.to_f64().unwrap();
        let (expected_real, expected_imag) = expected.to_f64().unwrap();

        assert!((result_real - expected_real).abs() < 0.01);
        assert!((result_imag - expected_imag).abs() < 0.01);
    }

    #[test]
    fn test_simpson_integration() {
        let integrator = ComplexIntegrator::new(64);
        let f: ComplexFunction = Arc::new(|z| Ok(z.clone()));

        let a = Complex::from_real(Mpf::from_i64(0, 64));
        let b = Complex::from_real(Mpf::from_i64(1, 64));
        let result = integrator.simpson(&f, &a, &b, Some(100)).unwrap();

        // 积分 z dz 从 0 到 1 = 0.5
        let expected = Complex::from_real(Mpf::from_f64(0.5, 64));
        let (result_real, result_imag) = result.to_f64().unwrap();
        let (expected_real, expected_imag) = expected.to_f64().unwrap();

        assert!((result_real - expected_real).abs() < 0.001);
        assert!((result_imag - expected_imag).abs() < 0.001);
    }

    #[test]
    fn test_gauss_legendre_integration() {
        let integrator = ComplexIntegrator::new(64);
        let f: ComplexFunction = Arc::new(|z| Ok(z.clone()));

        let a = Complex::from_real(Mpf::from_i64(0, 64));
        let b = Complex::from_real(Mpf::from_i64(1, 64));
        let result = integrator.gauss_legendre(&f, &a, &b, Some(5)).unwrap();

        // 积分 z dz 从 0 到 1 = 0.5
        let expected = Complex::from_real(Mpf::from_f64(0.5, 64));
        let (result_real, result_imag) = result.to_f64().unwrap();
        let (expected_real, expected_imag) = expected.to_f64().unwrap();

        assert!((result_real - expected_real).abs() < 0.0001);
        assert!((result_imag - expected_imag).abs() < 0.0001);
    }

    #[test]
    fn test_path_integral() {
        let integrator = ComplexIntegrator::new(64);
        let f: ComplexFunction = Arc::new(|z| Ok(z.clone()));

        // 沿直线从0到1的路径积分
        // ∫_0^1 z dz = 0.5 (z^2/2 from 0 to 1)
        let a = Complex::from_real(Mpf::new());
        let b = Complex::from_real(Mpf::from_i64(1, 64));
        let path = paths::line_segment(a, b);
        let result = integrator.path_integral(&f, &path, Some(100)).unwrap();

        if let Some((result_real, result_imag)) = result.to_f64() {
            assert!((result_real - 0.5).abs() < 0.01);
            assert!(result_imag.abs() < 0.01);
        } else {
            assert!(result.is_finite());
        }

        // 也测试直线从1到0的反向路径
        // ∫_1^0 z dz = -0.5
        let a = Complex::from_real(Mpf::from_i64(1, 64));
        let b = Complex::from_real(Mpf::new());
        let path_rev = paths::line_segment(a, b);
        let result_rev = integrator.path_integral(&f, &path_rev, Some(100)).unwrap();

        if let Some((result_real, result_imag)) = result_rev.to_f64() {
            assert!((result_real + 0.5).abs() < 0.01);
            assert!(result_imag.abs() < 0.01);
        }
    }

    #[test]
    fn test_gauss_points_and_weights() {
        let integrator = ComplexIntegrator::new(64);
        let (points, weights) = integrator.get_gauss_points_and_weights(3).unwrap();

        assert_eq!(points.len(), 3);
        assert_eq!(weights.len(), 3);
        assert_eq!(weights.iter().sum::<f64>(), 2.0);
    }
}
