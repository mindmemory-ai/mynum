//! 复数运算演示
//!
//! 演示MyNum库中的高精度复数运算功能

use mynum::{Complex, Mpf, Mpz};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== MyNum 复数运算演示 ===\n");

    // 基础复数操作
    basic_complex_operations()?;

    // 复数解析
    complex_parsing()?;

    // 复数算术运算
    complex_arithmetic()?;

    // 复数几何运算
    complex_geometry()?;

    // 复数函数
    complex_functions()?;

    // 欧拉公式验证
    euler_formula_verification()?;

    // 复数方程求解
    complex_equations()?;

    // 性能测试
    performance_tests()?;

    // 边界测试
    boundary_tests()?;

    // 信号处理示例
    signal_processing_example()?;

    // Julia集示例
    julia_set_example()?;

    println!("\n=== 演示完成 ===");
    Ok(())
}

fn basic_complex_operations() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔢 基础复数操作");
    println!("────────────────");

    // 创建复数
    let z1 = Complex::new();
    println!("零复数: {}", z1);

    let z2 = Complex::from_i64(3, 4, 64);
    println!("从整数创建: {}", z2);

    let real = Mpf::from_mpz(Mpz::from_i64(5), 64);
    let imag = Mpf::from_mpz(Mpz::from_i64(-2), 64);
    let z3 = Complex::from_real_imag(real, imag);
    println!("从实部虚部创建: {}", z3);

    let z4 = Complex::from_real(Mpf::from_mpz(Mpz::from_i64(7), 64));
    println!("纯实数: {}", z4);

    let z5 = Complex::from_imag(Mpf::from_mpz(Mpz::from_i64(3), 64));
    println!("纯虚数: {}", z5);

    // 复数属性
    println!("\n复数属性:");
    println!("{} 是否为零: {}", z2, z2.is_zero());
    println!("{} 是否为实数: {}", z2, z2.is_real());
    println!("{} 是否为纯虚数: {}", z2, z2.is_imaginary());
    println!("{} 是否为实数: {}", z4, z4.is_real());
    println!("{} 是否为纯虚数: {}", z5, z5.is_imaginary());

    println!();
    Ok(())
}

fn complex_parsing() -> Result<(), Box<dyn std::error::Error>> {
    println!("📝 复数解析");
    println!("────────────");

    let test_strings = [
        "3+4i", "5-2i", "7", "3i", "-2+5i", "6-i", "i", "-i", "-3-4i",
    ];

    for s in test_strings.iter() {
        match Complex::from_str(s, 64) {
            Ok(z) => println!("解析 '{}' -> {}", s, z),
            Err(e) => println!("解析 '{}' 失败: {:?}", s, e),
        }
    }

    println!();
    Ok(())
}

fn complex_arithmetic() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧮 复数算术运算");
    println!("────────────────");

    let z1 = Complex::from_i64(1, 2, 64);
    let z2 = Complex::from_i64(3, 4, 64);

    println!("z1 = {}", z1);
    println!("z2 = {}", z2);

    // 加法
    let sum = z1.add(&z2)?;
    println!("z1 + z2 = {}", sum);

    // 减法
    let diff = z1.sub(&z2)?;
    println!("z1 - z2 = {}", diff);

    // 乘法
    let product = z1.mul(&z2)?;
    println!("z1 * z2 = {}", product);

    // 除法
    let quotient = z1.div(&z2)?;
    println!("z1 / z2 = {}", quotient);

    // 验证交换律
    let ab = z1.mul(&z2)?;
    let ba = z2.mul(&z1)?;
    println!("z1*z2 = {}, z2*z1 = {} (交换律)", ab, ba);

    // 验证结合律
    let z3 = Complex::from_i64(5, 6, 64);
    let left = z1.add(&z2)?.add(&z3)?;
    let right = z1.add(&z2.add(&z3)?)?;
    println!("(z1 + z2) + z3 = {}", left);
    println!("z1 + (z2 + z3) = {}", right);

    println!();
    Ok(())
}

fn complex_geometry() -> Result<(), Box<dyn std::error::Error>> {
    println!("📐 复数几何运算");
    println!("────────────────");

    let z1 = Complex::from_i64(3, 4, 64);
    println!("z1 = {}", z1);

    // 模长
    match z1.magnitude() {
        Ok(mag) => println!("|z1| = {}", mag),
        Err(_) => println!("|z1| = 计算错误"),
    }

    // 幅角
    match z1.argument() {
        Ok(arg) => println!("arg(z1) = {}", arg),
        Err(_) => println!("arg(z1) = 计算错误"),
    }

    // 共轭复数
    let conjugate = z1.conjugate();
    println!("z1* = {}", conjugate);

    // 验证 |z|² = z * z*
    match z1.magnitude() {
        Ok(mag) => {
            let mag_squared = mag.mul(&mag);
            let z_conj = z1.mul(&conjugate)?;
            println!("|z1|² = {}", mag_squared);
            println!("z1 * z1* = {}", z_conj);
        }
        Err(_) => println!("模长计算错误"),
    }

    println!();
    Ok(())
}

fn complex_functions() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔬 复数函数");
    println!("────────────");

    let z = Complex::from_i64(2, 3, 64);
    println!("z = {}", z);

    // 幂运算
    for n in 1..=5 {
        match z.pow(n) {
            Ok(result) => println!("z^{} = {}", n, result),
            Err(_) => println!("z^{} = 计算错误", n),
        }
    }

    // 平方和立方
    match z.square() {
        Ok(square) => println!("z² = {}", square),
        Err(_) => println!("z² = 计算错误"),
    }

    match z.cube() {
        Ok(cube) => println!("z³ = {}", cube),
        Err(_) => println!("z³ = 计算错误"),
    }

    // 绝对值
    match z.abs() {
        Ok(abs_z) => println!("|z| = {}", abs_z),
        Err(_) => println!("|z| = 计算错误"),
    }

    println!();
    Ok(())
}

fn euler_formula_verification() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧮 欧拉公式验证");
    println!("────────────────");

    // 欧拉公式: e^(iπ) + 1 = 0
    // 注意：由于 Complex 没有 exp 方法，我们只能演示概念
    let pi = Mpf::from_f64(std::f64::consts::PI, 64);
    let i_pi = Complex::from_real_imag(Mpf::new(), pi);

    println!("iπ = {}", i_pi);
    println!("理论上 e^(iπ) = -1");
    println!("所以 e^(iπ) + 1 = 0");

    // 欧拉恒等式: e^(ix) = cos(x) + i*sin(x)
    let x = Mpf::from_f64(1.0, 64);
    let ix = Complex::from_real_imag(Mpf::new(), x);

    println!("ix = {}", ix);
    println!("理论上 e^(ix) = cos(x) + i*sin(x)");

    println!();
    Ok(())
}

fn complex_equations() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔍 复数方程求解");
    println!("────────────────");

    // 解方程 z² + 1 = 0
    println!("解方程: z² + 1 = 0");
    println!("理论解: z = ±i");

    let z1 = Complex::from_i64(0, 1, 64);
    let z2 = Complex::from_i64(0, -1, 64);
    println!("z₁ = {}", z1);
    println!("z₂ = {}", z2);

    // 验证解
    for z in [z1, z2] {
        let z_squared = z.mul(&z)?;
        let result = z_squared.add(&Complex::from_real(Mpf::from_i64(1, 64)))?;
        println!("验证 z = {}: z² + 1 = {}", z, result);
    }

    // 解方程 z² - 2z + 2 = 0
    println!("\n解方程: z² - 2z + 2 = 0");
    println!("理论解: z = 1 ± i");

    let z1 = Complex::from_real_imag(Mpf::from_i64(1, 64), Mpf::from_i64(1, 64));
    let z2 = Complex::from_real_imag(Mpf::from_i64(1, 64), Mpf::from_i64(-1, 64));
    println!("z₁ = {}", z1);
    println!("z₂ = {}", z2);

    // 验证
    for z in [z1, z2] {
        let z_squared = z.mul(&z)?;
        let neg_2z = z.mul(&Complex::from_real(Mpf::from_i64(-2, 64)))?;
        let result = z_squared
            .add(&neg_2z)?
            .add(&Complex::from_real(Mpf::from_i64(2, 64)))?;
        println!("验证 z = {}: z² - 2z + 2 = {}", z, result);
    }

    println!();
    Ok(())
}

fn performance_tests() -> Result<(), Box<dyn std::error::Error>> {
    println!("⚡ 性能测试");
    println!("────────────");

    let test_sizes = [100, 500, 1000];

    for size in test_sizes.iter() {
        println!("测试大小: {} 位", size);

        let z = Complex::from_i64(1, 1, *size);

        // 幂运算性能
        let start = Instant::now();
        let _ = z.pow(10);
        let elapsed = start.elapsed();
        println!("  10次幂运算: {:?}", elapsed);

        // 乘法性能
        let start = Instant::now();
        let _ = z.mul(&z);
        let elapsed = start.elapsed();
        println!("  乘法运算: {:?}", elapsed);

        println!();
    }

    Ok(())
}

fn boundary_tests() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔬 边界测试");
    println!("────────────");

    // 零复数测试
    let zero = Complex::new();
    println!("零复数: {}", zero);
    println!("  是否为零: {}", zero.is_zero());
    println!("  是否为实数: {}", zero.is_real());
    println!("  是否为纯虚数: {}", zero.is_imaginary());

    // 无穷大测试
    let large_real = Mpf::from_mpz(
        Mpz::from_str("9999999999999999999999999999999999999999", 10)?,
        64,
    );
    let large_complex = Complex::from_real(large_real);
    println!("大数复数: {}", large_complex);

    // 精度测试
    let high_precision = Complex::from_i64(1, 1, 256);
    println!("高精度复数: {}", high_precision);

    println!();
    Ok(())
}

fn signal_processing_example() -> Result<(), Box<dyn std::error::Error>> {
    println!("📡 信号处理示例");
    println!("────────────────");

    // 模拟信号
    let signal = Complex::from_i64(1, 0, 64);
    let frequency = Complex::from_i64(0, 1, 64);

    // 信号调制
    let result = signal.mul(&frequency)?;
    println!("信号 {} 与频率 {} 的乘积: {}", signal, frequency, result);

    // 阻抗计算
    let resistance = Mpf::from_i64(100, 64);
    let reactance = Mpf::from_i64(50, 64);
    let impedance = Complex::from_real_imag(resistance, reactance);
    println!("阻抗 Z = R + jX = {}", impedance);

    match impedance.magnitude() {
        Ok(mag) => {
            println!("阻抗模长: {}", mag);
        }
        Err(_) => println!("阻抗模长计算错误"),
    }

    println!();
    Ok(())
}

fn julia_set_example() -> Result<(), Box<dyn std::error::Error>> {
    println!("🎨 Julia集示例");
    println!("──────────────");

    let c = Complex::from_real_imag(Mpf::from_f64(-0.7, 64), Mpf::from_f64(0.27, 64));
    println!("参数 c = {}", c);

    let mut z = Complex::from_real_imag(Mpf::from_f64(0.0, 64), Mpf::from_f64(0.0, 64));

    println!("迭代过程:");
    for i in 1..=5 {
        z = z.mul(&z)?.add(&c)?;
        println!("  z_{} = {}", i, z);

        match z.magnitude() {
            Ok(mag) => {
                if mag > Mpf::from_i64(2, 64) {
                    println!("    发散!");
                    break;
                }
            }
            Err(_) => println!("    模长计算错误"),
        }
    }

    println!();
    Ok(())
}
