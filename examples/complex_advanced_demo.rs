//! 复数高级功能演示
//!
//! 展示复数多项式运算和数值积分功能

use mynum::complex::{
    integration::{paths, ComplexFunction, ComplexIntegrator},
    polynomial::ComplexPolynomial,
    Complex,
};
use mynum::mpf::Mpf;
use std::sync::Arc;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== 复数高级功能演示 ===\n");

    // 1. 复数多项式运算演示
    println!("1. 复数多项式运算");
    println!("==================");

    // 创建多项式 P(z) = 1 + 2z + z²
    let coeffs = vec![
        Complex::from_real(Mpf::from_i64(1, 64)),
        Complex::from_real(Mpf::from_i64(2, 64)),
        Complex::from_real(Mpf::from_i64(1, 64)),
    ];
    let poly1 = ComplexPolynomial::from_coefficients(coeffs)?;
    println!("多项式 P(z) = {}", poly1);
    println!("次数: {}", poly1.degree());

    // 创建多项式 Q(z) = 1 + z
    let coeffs2 = vec![
        Complex::from_real(Mpf::from_i64(1, 64)),
        Complex::from_real(Mpf::from_i64(1, 64)),
    ];
    let poly2 = ComplexPolynomial::from_coefficients(coeffs2)?;
    println!("多项式 Q(z) = {}", poly2);

    // 多项式加法
    let sum = poly1.add(&poly2)?;
    println!("P(z) + Q(z) = {}", sum);

    // 多项式乘法
    let product = poly1.mul(&poly2)?;
    println!("P(z) × Q(z) = {}", product);

    // 多项式求值
    let z = Complex::from_real(Mpf::from_i64(3, 64));
    let value = poly1.evaluate(&z)?;
    println!("P({}) = {}", z, value);

    // 多项式导数
    let derivative = poly1.derivative()?;
    println!("P'(z) = {}", derivative);

    // 多项式积分
    let constant = Complex::from_real(Mpf::from_i64(0, 64));
    let integral = poly1.integral(&constant)?;
    println!("∫P(z)dz = {}", integral);

    println!();

    // 2. 复数数值积分演示
    println!("2. 复数数值积分");
    println!("================");

    let integrator = ComplexIntegrator::new(64);

    // 定义被积函数 f(z) = z
    let f: ComplexFunction = Arc::new(|z| Ok(z.clone()));

    // 积分区间 [0, 1]
    let a = Complex::from_real(Mpf::from_i64(0, 64));
    let b = Complex::from_real(Mpf::from_i64(1, 64));

    // 梯形法则
    let trapezoidal_result = integrator.trapezoidal(&f, &a, &b, Some(100))?;
    println!("梯形法则积分 ∫z dz [0,1] = {}", trapezoidal_result);

    // 辛普森法则
    let simpson_result = integrator.simpson(&f, &a, &b, Some(100))?;
    println!("辛普森法则积分 ∫z dz [0,1] = {}", simpson_result);

    // 高斯求积
    let gauss_result = integrator.gauss_legendre(&f, &a, &b, Some(5))?;
    println!("高斯求积积分 ∫z dz [0,1] = {}", gauss_result);

    // 自适应积分
    let adaptive_result = integrator.adaptive(&f, &a, &b, 1e-6, Some(10))?;
    println!("自适应积分 ∫z dz [0,1] = {}", adaptive_result);

    println!();

    // 3. 路径积分演示
    println!("3. 路径积分");
    println!("============");

    // 沿单位圆的路径积分
    let unit_circle_path: ComplexFunction = Arc::new(|t| paths::unit_circle(t));
    let circle_result = integrator.path_integral(&f, &unit_circle_path, Some(100))?;
    println!("沿单位圆的路径积分 ∫z dz = {}", circle_result);

    // 沿直线段的路径积分
    let line_path = paths::line_segment(a.clone(), b.clone());
    let line_result = integrator.path_integral(&f, &line_path, Some(100))?;
    println!("沿直线段的路径积分 ∫z dz [0,1] = {}", line_result);

    // 沿矩形的路径积分
    let center = Complex::from_real(Mpf::from_i64(0, 64));
    let width = Mpf::from_i64(2, 64);
    let height = Mpf::from_i64(2, 64);
    let rect_path = paths::rectangle(center, width, height);
    let rect_result = integrator.path_integral(&f, &rect_path, Some(100))?;
    println!("沿矩形的路径积分 ∫z dz = {}", rect_result);

    println!();

    // 4. 更复杂的函数积分
    println!("4. 复杂函数积分");
    println!("================");

    // 定义函数 f(z) = z²
    let f_square: ComplexFunction = Arc::new(|z| Ok(z.mul(z)?));

    // 积分 ∫z² dz [0,1] = 1/3
    let square_result = integrator.simpson(&f_square, &a, &b, Some(100))?;
    println!("积分 ∫z² dz [0,1] = {}", square_result);

    // 定义函数 f(z) = e^z
    let f_exp: ComplexFunction = Arc::new(|z| Ok(z.exp()?));

    // 积分 ∫e^z dz [0,1] = e - 1
    let exp_result = integrator.simpson(&f_exp, &a, &b, Some(100))?;
    println!("积分 ∫e^z dz [0,1] = {}", exp_result);

    println!();
    println!("演示完成！");

    Ok(())
}
