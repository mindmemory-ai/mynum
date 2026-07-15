//! 复数基本功能演示
//!
//! 展示新的复数模块的基本功能

use mynum::{Complex, Mpf};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== 复数基本功能演示 ===\n");

    // 创建复数
    println!("1. 创建复数:");
    let z1 = Complex::from_str("3+4i", 64)?;
    let z2 = Complex::from_str("1+2i", 64)?;
    println!("   z1 = {}", z1);
    println!("   z2 = {}\n", z2);

    // 基本运算
    println!("2. 基本运算:");
    let sum = z1.add(&z2)?;
    let diff = z1.sub(&z2)?;
    let product = z1.mul(&z2)?;
    let quotient = z1.div(&z2)?;

    println!("   z1 + z2 = {}", sum);
    println!("   z1 - z2 = {}", diff);
    println!("   z1 * z2 = {}", product);
    println!("   z1 / z2 = {}\n", quotient);

    // 复数属性
    println!("3. 复数属性:");
    match z1.magnitude() {
        Ok(mag) => println!("   z1 的模长: {}", mag),
        Err(_) => println!("   z1 的模长: 计算错误"),
    }
    match z1.argument() {
        Ok(arg) => println!("   z1 的幅角: {}", arg),
        Err(_) => println!("   z1 的幅角: 计算错误"),
    }
    println!("   z1 的共轭: {}", z1.conjugate());
    match z1.reciprocal() {
        Ok(rec) => println!("   z1 的倒数: {}", rec),
        Err(_) => println!("   z1 的倒数: 计算错误"),
    }
    println!();

    // 幂运算
    println!("4. 幂运算:");
    match z1.pow(2) {
        Ok(square) => println!("   z1² = {}", square),
        Err(_) => println!("   z1² = 计算错误"),
    }
    match z1.pow(3) {
        Ok(cube) => println!("   z1³ = {}", cube),
        Err(_) => println!("   z1³ = 计算错误"),
    }
    match z1.pow(3) {
        Ok(power_3) => println!("   z1^3 = {}", power_3),
        Err(_) => println!("   z1^3 = 计算错误"),
    }
    println!();

    // 常量
    println!("5. 复数常量:");
    println!("   0 = {}", Complex::new());
    println!("   1 = {}", Complex::from_real(Mpf::from_i64(1, 64)));
    println!("   i = {}", Complex::from_imag(Mpf::from_i64(1, 64)));
    println!("   -1 = {}", Complex::from_real(Mpf::from_i64(-1, 64)));
    println!("   -i = {}", Complex::from_imag(Mpf::from_i64(-1, 64)));
    println!();

    // 特殊值
    println!("6. 特殊值:");
    let pi_complex = Complex::from_real(Mpf::from_f64(std::f64::consts::PI, 64));
    let e_complex = Complex::from_real(Mpf::from_f64(std::f64::consts::E, 64));
    println!("   π ≈ {}", pi_complex);
    println!("   e ≈ {}\n", e_complex);

    // 从极坐标创建
    println!("7. 极坐标:");
    match z1.magnitude() {
        Ok(r) => match z1.argument() {
            Ok(theta) => {
                println!("   z1 的极坐标: r = {}, θ = {}", r, theta);

                match Complex::from_polar(r, theta, 64) {
                    Ok(z1_from_polar) => {
                        println!("   从极坐标重建: {}", z1_from_polar);
                        match z1.sub(&z1_from_polar) {
                            Ok(error) => println!("   重建误差: {}", error),
                            Err(_) => println!("   重建误差: 计算错误"),
                        }
                    }
                    Err(_) => println!("   从极坐标重建: 计算错误"),
                }
            }
            Err(_) => println!("   z1 的幅角: 计算错误"),
        },
        Err(_) => println!("   z1 的模长: 计算错误"),
    }
    println!();

    // 统计函数 - 简化版本
    println!("8. 统计函数:");
    let complexes = vec![z1.clone(), z2.clone()];
    let avg_real = complexes
        .iter()
        .map(|z| z.real().clone())
        .fold(Mpf::new(), |acc, x| acc.add(&x))
        .div(&Mpf::from_i64(complexes.len() as i64, 64))?;
    let avg_imag = complexes
        .iter()
        .map(|z| z.imaginary().clone())
        .fold(Mpf::new(), |acc, x| acc.add(&x))
        .div(&Mpf::from_i64(complexes.len() as i64, 64))?;
    let avg = Complex::from_real_imag(avg_real, avg_imag);
    println!("   平均值: {}", avg);
    println!();

    // 距离和内积 - 简化版本
    println!("9. 几何运算:");
    match z1.distance(&z2) {
        Ok(distance) => println!("   z1 和 z2 的距离: {}", distance),
        Err(_) => println!("   z1 和 z2 的距离: 计算错误"),
    }

    // 简化的内积计算
    let dot_product = z1
        .real()
        .mul(z2.real())
        .add(&z1.imaginary().mul(z2.imaginary()));
    println!("   z1 和 z2 的内积: {}", dot_product);
    println!();

    // 归一化
    println!("10. 归一化:");
    match z1.normalize() {
        Ok(normalized) => {
            println!("   z1 归一化后: {}", normalized);
            match normalized.magnitude() {
                Ok(mag) => println!("   归一化后的模长: {}", mag),
                Err(_) => println!("   归一化后的模长: 计算错误"),
            }
        }
        Err(_) => println!("   z1 归一化: 计算错误"),
    }
    println!();

    println!("=== 演示完成 ===");
    Ok(())
}
