//! 特殊数学函数演示
//!
//! 演示MyNum库中的Gamma函数、贝塞尔函数、误差函数等高级数学函数

use mynum::mpf::special::{bessel_j0, bessel_j1, erf, erfc, gamma};
use mynum::{Mpf, Mpz};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== MyNum 特殊函数演示 ===\n");

    // 测试Gamma函数
    test_gamma_function()?;

    // 测试贝塞尔函数
    test_bessel_functions()?;

    // 测试误差函数
    test_error_functions()?;

    // 性能测试
    performance_tests()?;

    println!("\n=== 演示完成 ===");
    Ok(())
}

fn test_gamma_function() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔬 测试Gamma函数");
    println!("────────────────");

    let precision = 64;

    // 测试整数值
    let test_values = [1, 2, 3, 4, 5];
    let expected_results = [1, 1, 2, 6, 24]; // Γ(n) = (n-1)!

    for (i, &n) in test_values.iter().enumerate() {
        let x = Mpf::from_mpz(Mpz::from_i64(n), precision);
        match gamma(&x) {
            Ok(result) => {
                let result_str = result.to_string(10);
                println!("Γ({}) = {} (期望: {})", n, result_str, expected_results[i]);
            }
            Err(e) => println!("Γ({}) 计算失败: {:?}", n, e),
        }
    }

    // 测试半整数值
    println!("\n半整数值:");
    let half = Mpf::from_mpz(Mpz::from_i64(1), precision)
        .div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?;
    match gamma(&half) {
        Ok(result) => {
            let result_str = result.to_string(10);
            println!("Γ(1/2) = {} (期望: √π ≈ 1.7724)", result_str);
        }
        Err(e) => println!("Γ(1/2) 计算失败: {:?}", e),
    }

    // 测试对数Gamma函数
    println!("\n对数Gamma函数:");
    let large_x = Mpf::from_mpz(Mpz::from_i64(10), precision);
    // 使用ln(gamma(x))代替log_gamma
    match gamma(&large_x) {
        Ok(gamma_val) => match gamma_val.ln() {
            Ok(result) => {
                let result_str = result.to_string(10);
                println!("ln(Γ(10)) = {}", result_str);
            }
            Err(e) => println!("ln(Γ(10)) 计算失败: {:?}", e),
        },
        Err(e) => println!("Γ(10) 计算失败: {:?}", e),
    }

    println!();
    Ok(())
}

fn test_bessel_functions() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔬 测试贝塞尔函数");
    println!("────────────────");

    let precision = 64;

    // 测试J₀(x)函数
    println!("J₀(x) 函数:");
    let test_values = [0.0, 1.0, 2.0, 3.0, 5.0];

    for &x_val in test_values.iter() {
        let x = Mpf::from_f64(x_val, precision);
        match bessel_j0(&x) {
            Ok(result) => {
                let result_str = result.to_string(10);
                println!("J₀({}) = {}", x_val, result_str);
            }
            Err(e) => println!("J₀({}) 计算失败: {:?}", x_val, e),
        }
    }

    // 测试J₁(x)函数
    println!("\nJ₁(x) 函数:");
    for &x_val in test_values.iter() {
        let x = Mpf::from_f64(x_val, precision);
        match bessel_j1(&x) {
            Ok(result) => {
                let result_str = result.to_string(10);
                println!("J₁({}) = {}", x_val, result_str);
            }
            Err(e) => println!("J₁({}) 计算失败: {:?}", x_val, e),
        }
    }

    println!();
    Ok(())
}

fn test_error_functions() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔬 测试误差函数");
    println!("────────────────");

    let precision = 64;

    // 测试erf函数
    println!("erf(x) 函数:");
    let test_values = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0];

    for &x_val in test_values.iter() {
        let x = Mpf::from_f64(x_val, precision);
        match erf(&x) {
            Ok(result) => {
                let result_str = result.to_string(10);
                println!("erf({}) = {}", x_val, result_str);
            }
            Err(e) => println!("erf({}) 计算失败: {:?}", x_val, e),
        }
    }

    // 测试erfc函数
    println!("\nerfc(x) 函数:");
    for &x_val in test_values.iter() {
        let x = Mpf::from_f64(x_val, precision);
        match erfc(&x) {
            Ok(result) => {
                let result_str = result.to_string(10);
                println!("erfc({}) = {}", x_val, result_str);
            }
            Err(e) => println!("erfc({}) 计算失败: {:?}", x_val, e),
        }
    }

    // 验证erf(x) + erfc(x) = 1
    println!("\n验证 erf(x) + erfc(x) = 1:");
    let test_x = Mpf::from_f64(1.5, precision);
    if let (Ok(erf_val), Ok(erfc_val)) = (erf(&test_x), erfc(&test_x)) {
        let sum = erf_val.add(&erfc_val);
        println!(
            "erf(1.5) + erfc(1.5) = {} + {} = {}",
            erf_val.to_string(10),
            erfc_val.to_string(10),
            sum.to_string(10)
        );
    }

    println!();
    Ok(())
}

fn performance_tests() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔬 性能测试");
    println!("──────────");

    let precision = 64;
    let test_x = Mpf::from_f64(2.5, precision);
    let iterations = 100;

    // Gamma函数性能测试
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = gamma(&test_x);
    }
    let gamma_time = start.elapsed();
    println!("Gamma函数 ({}次): {:?}", iterations, gamma_time);

    // 贝塞尔函数性能测试
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = bessel_j0(&test_x);
    }
    let bessel_time = start.elapsed();
    println!("J₀函数 ({}次): {:?}", iterations, bessel_time);

    // 误差函数性能测试
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = erf(&test_x);
    }
    let erf_time = start.elapsed();
    println!("erf函数 ({}次): {:?}", iterations, erf_time);

    println!();
    Ok(())
}
