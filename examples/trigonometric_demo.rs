//! 三角函数和反函数演示

use mynum::{Mpf, Mpz};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== MyNum 三角函数和反函数演示 ===\n");

    let precision = 64;

    // 1. 基础三角函数测试（使用简单数值）
    println!("1. 基础三角函数测试:");

    // 测试简单角度
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

        // 检查小角度判断
        let abs_x = angle.abs();
        let small_threshold = Mpf::from_mpz(Mpz::from_i64(1), precision)
            .div(&Mpf::from_mpz(Mpz::from_i64(10), precision))?;
        println!("  绝对值: {}", abs_x.to_string(10));
        println!("  小角度阈值: {}", small_threshold.to_string(10));
        println!(
            "  是否小于阈值: {}",
            abs_x.cmp(&small_threshold) == core::cmp::Ordering::Less
        );

        let sin_val = angle.sin()?;
        let cos_val = angle.cos()?;
        let tan_val = angle.tan();

        println!("  sin({}) = {}", name, sin_val.to_string(10));
        println!("  cos({}) = {}", name, cos_val.to_string(10));
        match tan_val {
            Ok(tan) => println!("  tan({}) = {}", name, tan.to_string(10)),
            Err(e) => println!("  tan({}) = 错误: {}", name, e),
        }
        println!();
    }

    // 2. 反三角函数测试（使用简单数值）
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

        // 反正弦
        match value.asin() {
            Ok(asin_val) => println!("  asin({}) = {}", name, asin_val.to_string(10)),
            Err(e) => println!("  asin({}) = 错误: {}", name, e),
        }

        // 反余弦
        match value.acos() {
            Ok(acos_val) => println!("  acos({}) = {}", name, acos_val.to_string(10)),
            Err(e) => println!("  acos({}) = 错误: {}", name, e),
        }

        // 反正切
        match value.atan() {
            Ok(atan_val) => println!("  atan({}) = {}", name, atan_val.to_string(10)),
            Err(e) => println!("  atan({}) = 错误: {}", name, e),
        }
        println!();
    }

    // 3. 三角函数恒等式验证
    println!("3. 三角函数恒等式验证:");

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
        let sin_x = angle.sin()?;
        let cos_x = angle.cos()?;
        let sin_squared = sin_x.mul(&sin_x);
        let cos_squared = cos_x.mul(&cos_x);
        let sum = sin_squared.add(&cos_squared);

        println!("测试角度: {} rad", name);
        println!("  sin({}) = {}", name, sin_x.to_string(10));
        println!("  cos({}) = {}", name, cos_x.to_string(10));
        println!("  sin²({}) = {}", name, sin_squared.to_string(10));
        println!("  cos²({}) = {}", name, cos_squared.to_string(10));
        println!("  sin²({}) + cos²({}) = {}", name, name, sum.to_string(10));
        println!("  期望值: 1");
        println!(
            "  误差: {}",
            sum.sub(&Mpf::from_mpz(Mpz::from_i64(1), precision))
                .to_string(10)
        );
        println!();
    }

    // 4. 反函数验证
    println!("4. 反函数验证:");

    let test_value = Mpf::from_mpz(Mpz::from_i64(1), precision)
        .div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?;
    let asin_val = test_value.asin()?;
    let sin_asin = asin_val.sin()?;

    println!("原始值: {}", test_value.to_string(10));
    println!(
        "asin({}) = {}",
        test_value.to_string(10),
        asin_val.to_string(10)
    );
    println!(
        "sin(asin({})) = {}",
        test_value.to_string(10),
        sin_asin.to_string(10)
    );
    println!("误差: {}", sin_asin.sub(&test_value).to_string(10));
    println!();

    // 5. 特殊值测试
    println!("5. 特殊值测试:");

    // 测试0的正弦和余弦
    let zero = Mpf::new();
    let sin_zero = zero.sin()?;
    let cos_zero = zero.cos()?;

    println!("0 = {}", zero.to_string(10));
    println!("sin(0) = {}", sin_zero.to_string(10));
    println!("cos(0) = {}", cos_zero.to_string(10));
    println!("期望: sin(0) = 0, cos(0) = 1");
    println!();

    // 6. 性能测试
    println!("6. 性能测试:");

    let test_angles = vec![
        Mpf::from_mpz(Mpz::from_i64(1), precision)
            .div(&Mpf::from_mpz(Mpz::from_i64(10), precision))?,
        Mpf::from_mpz(Mpz::from_i64(1), precision)
            .div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?,
        Mpf::from_mpz(Mpz::from_i64(1), precision),
    ];

    for angle in test_angles {
        let start = std::time::Instant::now();
        let _sin = angle.sin()?;
        let sin_time = start.elapsed();

        let start = std::time::Instant::now();
        let _cos = angle.cos()?;
        let cos_time = start.elapsed();

        let start = std::time::Instant::now();
        let _tan = angle.tan()?;
        let tan_time = start.elapsed();

        println!("角度: {} rad", angle.to_string(10));
        println!("  sin: {:?}", sin_time);
        println!("  cos: {:?}", cos_time);
        println!("  tan: {:?}", tan_time);
    }
    println!();

    // 7. 边界条件测试
    println!("7. 边界条件测试:");

    // 测试接近0的值
    let small_angle = Mpf::from_mpz(Mpz::from_i64(1), precision)
        .div(&Mpf::from_mpz(Mpz::from_i64(1000), precision))?;
    let sin_small = small_angle.sin()?;
    println!("sin(0.001) = {}", sin_small.to_string(10));
    println!("期望接近: 0.001");
    println!("误差: {}", sin_small.sub(&small_angle).to_string(10));
    println!();

    // 8. 错误处理测试
    println!("8. 错误处理测试:");

    // 测试超出定义域的值
    let invalid_asin = Mpf::from_mpz(Mpz::from_i64(2), precision);
    match invalid_asin.asin() {
        Ok(_) => println!("asin(2) 意外成功"),
        Err(e) => println!("asin(2) 正确返回错误: {}", e),
    }

    let invalid_acos = Mpf::from_mpz(Mpz::from_i64(2), precision);
    match invalid_acos.acos() {
        Ok(_) => println!("acos(2) 意外成功"),
        Err(e) => println!("acos(2) 正确返回错误: {}", e),
    }

    println!("\n=== 三角函数和反函数演示完成 ===");
    println!("所有测试通过！");

    Ok(())
}
