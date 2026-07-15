//! 双曲函数和反双曲函数演示程序

use mynum::{Mpf, Mpz};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== MyNum 双曲函数和反双曲函数演示 ===\n");

    let precision = 64;

    // 1. 基础双曲函数测试
    println!("1. 基础双曲函数测试:");

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

        let sinh_val = value.sinh()?;
        let cosh_val = value.cosh()?;
        let tanh_val = value.tanh()?;

        println!("  sinh({}) = {}", name, sinh_val.to_string(10));
        println!("  cosh({}) = {}", name, cosh_val.to_string(10));
        println!("  tanh({}) = {}", name, tanh_val.to_string(10));
        println!();
    }

    // 2. 反双曲函数测试
    println!("2. 反双曲函数测试:");

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

        let asinh_val = value.asinh()?;
        let acosh_val = if value.cmp(&Mpf::from_mpz(Mpz::from_i64(1), precision))
            == core::cmp::Ordering::Greater
        {
            value.acosh()
        } else {
            Err(mynum::Error::DomainError("Acosh domain error".into()))
        };
        let atanh_val = if value.cmp(&Mpf::from_mpz(Mpz::from_i64(1), precision))
            == core::cmp::Ordering::Less
        {
            value.atanh()
        } else {
            Err(mynum::Error::DomainError("Atanh domain error".into()))
        };

        println!("  asinh({}) = {}", name, asinh_val.to_string(10));
        match acosh_val {
            Ok(acosh) => println!("  acosh({}) = {}", name, acosh.to_string(10)),
            Err(e) => println!("  acosh({}) = 错误: {}", name, e),
        }
        match atanh_val {
            Ok(atanh) => println!("  atanh({}) = {}", name, atanh.to_string(10)),
            Err(e) => println!("  atanh({}) = 错误: {}", name, e),
        }
        println!();
    }

    // 3. 双曲函数恒等式验证
    println!("3. 双曲函数恒等式验证:");

    let test_values = vec![
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
        let sinh_x = value.sinh()?;
        let cosh_x = value.cosh()?;
        let sinh_squared = sinh_x.mul(&sinh_x);
        let cosh_squared = cosh_x.mul(&cosh_x);
        let difference = cosh_squared.sub(&sinh_squared);

        println!("测试值: {} = {}", name, value.to_string(10));
        println!("  sinh({}) = {}", name, sinh_x.to_string(10));
        println!("  cosh({}) = {}", name, cosh_x.to_string(10));
        println!(
            "  cosh²({}) - sinh²({}) = {}",
            name,
            name,
            difference.to_string(10)
        );
        println!("  期望值: 1");
        println!(
            "  误差: {}",
            difference
                .sub(&Mpf::from_mpz(Mpz::from_i64(1), precision))
                .to_string(10)
        );
        println!();
    }

    // 4. 反函数验证
    println!("4. 反函数验证:");

    let test_value = Mpf::from_mpz(Mpz::from_i64(1), precision)
        .div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?;
    let sinh_val = test_value.sinh()?;
    let asinh_sinh = sinh_val.asinh()?;

    println!("原始值: {}", test_value.to_string(10));
    println!(
        "sinh({}) = {}",
        test_value.to_string(10),
        sinh_val.to_string(10)
    );
    println!(
        "asinh(sinh({})) = {}",
        test_value.to_string(10),
        asinh_sinh.to_string(10)
    );
    println!("误差: {}", asinh_sinh.sub(&test_value).to_string(10));
    println!();

    // 5. 性能测试
    println!("5. 性能测试:");

    let test_value = Mpf::from_mpz(Mpz::from_i64(1), precision)
        .div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?;

    let start = std::time::Instant::now();
    for _ in 0..100 {
        let _sinh = test_value.sinh()?;
    }
    let sinh_time = start.elapsed();

    let start = std::time::Instant::now();
    for _ in 0..100 {
        let _cosh = test_value.cosh()?;
    }
    let cosh_time = start.elapsed();

    let start = std::time::Instant::now();
    for _ in 0..100 {
        let _tanh = test_value.tanh()?;
    }
    let tanh_time = start.elapsed();

    println!(
        "测试值: {} = {}",
        test_value.to_string(10),
        test_value.to_string(10)
    );
    println!("  sinh (100次): {:?}", sinh_time);
    println!("  cosh (100次): {:?}", cosh_time);
    println!("  tanh (100次): {:?}", tanh_time);
    println!();

    // 6. 边界条件测试
    println!("6. 边界条件测试:");

    // 测试接近0的值
    let small_value = Mpf::from_mpz(Mpz::from_i64(1), precision)
        .div(&Mpf::from_mpz(Mpz::from_i64(1000), precision))?;
    let sinh_small = small_value.sinh()?;
    println!("sinh(0.001) = {}", sinh_small.to_string(10));
    println!("期望接近: 0.001");
    println!("误差: {}", sinh_small.sub(&small_value).to_string(10));
    println!();

    // 7. 错误处理测试
    println!("7. 错误处理测试:");

    // 测试acosh的定义域
    let invalid_acosh = Mpf::from_mpz(Mpz::from_i64(1), precision)
        .div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?;
    match invalid_acosh.acosh() {
        Ok(_) => println!("acosh(0.5) 意外成功"),
        Err(e) => println!("acosh(0.5) 正确返回错误: {}", e),
    }

    // 测试atanh的定义域
    let invalid_atanh = Mpf::from_mpz(Mpz::from_i64(2), precision);
    match invalid_atanh.atanh() {
        Ok(_) => println!("atanh(2) 意外成功"),
        Err(e) => println!("atanh(2) 正确返回错误: {}", e),
    }

    println!("\n=== 双曲函数和反双曲函数演示完成 ===");
    println!("双曲函数提供了完整的双曲三角函数支持！");

    Ok(())
}
