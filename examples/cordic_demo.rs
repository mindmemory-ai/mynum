//! CORDIC算法演示程序

use mynum::mpf::{Cordic, CordicConfig};
use mynum::{Mpf, Mpz};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== MyNum CORDIC算法演示 ===\n");

    let precision = 64;

    // 1. 基础CORDIC测试
    println!("1. 基础CORDIC测试:");

    let config = CordicConfig::default();
    let cordic = Cordic::new(config)?;

    let test_angles = vec![
        ("0", Mpf::new()),
        (
            "0.1",
            Mpf::from_mpz(Mpz::from_i64(1), precision)
                .div(&Mpf::from_mpz(Mpz::from_i64(10), precision))?,
        ),
        (
            "0.5",
            Mpf::from_mpz(Mpz::from_i64(1), precision)
                .div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?,
        ),
        ("1.0", Mpf::from_mpz(Mpz::from_i64(1), precision)),
    ];

    for (name, angle) in test_angles {
        println!("测试角度: {} = {}", name, angle.to_string(10));

        let sin_val = cordic.sin(&angle)?;
        let cos_val = cordic.cos(&angle)?;
        let tan_val = cordic.tan(&angle);

        println!("  CORDIC sin({}) = {}", name, sin_val.to_string(10));
        println!("  CORDIC cos({}) = {}", name, cos_val.to_string(10));
        match tan_val {
            Ok(tan) => println!("  CORDIC tan({}) = {}", name, tan.to_string(10)),
            Err(e) => println!("  CORDIC tan({}) = 错误: {}", name, e),
        }
        println!();
    }

    // 2. 反三角函数测试
    println!("2. 反三角函数测试:");

    let test_values = vec![
        ("0", Mpf::new()),
        (
            "0.1",
            Mpf::from_mpz(Mpz::from_i64(1), precision)
                .div(&Mpf::from_mpz(Mpz::from_i64(10), precision))?,
        ),
        (
            "0.5",
            Mpf::from_mpz(Mpz::from_i64(1), precision)
                .div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?,
        ),
        ("1.0", Mpf::from_mpz(Mpz::from_i64(1), precision)),
    ];

    for (name, value) in test_values {
        println!("测试值: {} = {}", name, value.to_string(10));

        let asin_val = cordic.asin(&value)?;
        let acos_val = cordic.acos(&value)?;
        let atan_val = cordic.atan(&value)?;

        println!("  CORDIC asin({}) = {}", name, asin_val.to_string(10));
        println!("  CORDIC acos({}) = {}", name, acos_val.to_string(10));
        println!("  CORDIC atan({}) = {}", name, atan_val.to_string(10));
        println!();
    }

    // 3. 性能对比测试
    println!("3. 性能对比测试:");

    let test_angle = Mpf::from_mpz(Mpz::from_i64(1), precision)
        .div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?;

    // 测试CORDIC性能
    let start = std::time::Instant::now();
    for _ in 0..100 {
        let _sin = cordic.sin(&test_angle)?;
    }
    let cordic_time = start.elapsed();

    // 测试传统算法性能
    let start = std::time::Instant::now();
    for _ in 0..100 {
        let _sin = test_angle.sin()?;
    }
    let traditional_time = start.elapsed();

    println!("角度: {} rad", test_angle.to_string(10));
    println!("  CORDIC sin (100次): {:?}", cordic_time);
    println!("  传统算法 sin (100次): {:?}", traditional_time);
    println!(
        "  性能提升: {:.2}x",
        traditional_time.as_nanos() as f64 / cordic_time.as_nanos() as f64
    );
    println!();

    // 4. 精度对比测试
    println!("4. 精度对比测试:");

    let test_angle = Mpf::from_mpz(Mpz::from_i64(1), precision)
        .div(&Mpf::from_mpz(Mpz::from_i64(4), precision))?;

    let cordic_sin = cordic.sin(&test_angle)?;
    let traditional_sin = test_angle.sin()?;

    println!("角度: {} rad", test_angle.to_string(10));
    println!("  CORDIC sin: {}", cordic_sin.to_string(10));
    println!("  传统算法 sin: {}", traditional_sin.to_string(10));
    println!("  差异: {}", cordic_sin.sub(&traditional_sin).to_string(10));
    println!();

    // 5. 不同配置的CORDIC测试
    println!("5. 不同配置的CORDIC测试:");

    let configs = vec![
        (
            "低精度",
            CordicConfig {
                iterations: 16,
                double_precision: false,
            },
        ),
        (
            "标准精度",
            CordicConfig {
                iterations: 32,
                double_precision: false,
            },
        ),
        (
            "高精度",
            CordicConfig {
                iterations: 64,
                double_precision: false,
            },
        ),
    ];

    let test_angle = Mpf::from_mpz(Mpz::from_i64(1), precision)
        .div(&Mpf::from_mpz(Mpz::from_i64(3), precision))?;

    for (name, config) in configs {
        let cordic = Cordic::new(config)?;
        let start = std::time::Instant::now();
        let sin_val = cordic.sin(&test_angle)?;
        let time = start.elapsed();

        println!(
            "{}: sin({}) = {} (耗时: {:?})",
            name,
            test_angle.to_string(10),
            sin_val.to_string(10),
            time
        );
    }
    println!();

    // 6. 三角函数恒等式验证
    println!("6. 三角函数恒等式验证:");

    let test_angles = vec![
        (
            "0.1",
            Mpf::from_mpz(Mpz::from_i64(1), precision)
                .div(&Mpf::from_mpz(Mpz::from_i64(10), precision))?,
        ),
        (
            "0.5",
            Mpf::from_mpz(Mpz::from_i64(1), precision)
                .div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?,
        ),
        ("1.0", Mpf::from_mpz(Mpz::from_i64(1), precision)),
    ];

    for (name, angle) in test_angles {
        let sin_x = cordic.sin(&angle)?;
        let cos_x = cordic.cos(&angle)?;
        let sin_squared = sin_x.mul(&sin_x);
        let cos_squared = cos_x.mul(&cos_x);
        let sum = sin_squared.add(&cos_squared);

        println!("测试角度: {} rad", name);
        println!("  CORDIC sin({}) = {}", name, sin_x.to_string(10));
        println!("  CORDIC cos({}) = {}", name, cos_x.to_string(10));
        println!("  sin²({}) + cos²({}) = {}", name, name, sum.to_string(10));
        println!("  期望值: 1");
        println!(
            "  误差: {}",
            sum.sub(&Mpf::from_mpz(Mpz::from_i64(1), precision))
                .to_string(10)
        );
        println!();
    }

    // 7. 边界条件测试
    println!("7. 边界条件测试:");

    // 测试接近0的值
    let small_angle = Mpf::from_mpz(Mpz::from_i64(1), precision)
        .div(&Mpf::from_mpz(Mpz::from_i64(1000), precision))?;
    let sin_small = cordic.sin(&small_angle)?;
    println!("CORDIC sin(0.001) = {}", sin_small.to_string(10));
    println!("期望接近: 0.001");
    println!("误差: {}", sin_small.sub(&small_angle).to_string(10));
    println!();

    // 8. 错误处理测试
    println!("8. 错误处理测试:");

    // 测试超出定义域的值
    let invalid_asin = Mpf::from_mpz(Mpz::from_i64(2), precision);
    match cordic.asin(&invalid_asin) {
        Ok(_) => println!("CORDIC asin(2) 意外成功"),
        Err(e) => println!("CORDIC asin(2) 正确返回错误: {}", e),
    }

    let invalid_acos = Mpf::from_mpz(Mpz::from_i64(2), precision);
    match cordic.acos(&invalid_acos) {
        Ok(_) => println!("CORDIC acos(2) 意外成功"),
        Err(e) => println!("CORDIC acos(2) 正确返回错误: {}", e),
    }

    println!("\n=== CORDIC算法演示完成 ===");
    println!("CORDIC算法提供了高性能的三角函数计算！");

    Ok(())
}
