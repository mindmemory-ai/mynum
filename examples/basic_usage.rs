//! 基本使用示例

use mynum::{Mpz, MpzMultiplicationConfig, MultiplicationBackend};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== MyNum 大数运算库示例 ===\n");

    // 基本运算
    println!("1. 基本运算演示:");
    let a = Mpz::from_str("123456789012345678901234567890", 10)?;
    let b = Mpz::from_str("987654321098765432109876543210", 10)?;

    println!("a = {}", a);
    println!("b = {}", b);

    let sum = a.add(&b);
    println!("a + b = {}", sum);

    let diff = a.sub(&b);
    println!("a - b = {}", diff);

    let product = a.mul(&b);
    println!("a * b = {}", product);

    // 幂运算
    println!("\n2. 幂运算演示:");
    let base = Mpz::from_i64(2);
    let power = base.pow_u32(100);
    println!("2^100 = {}", power);

    // 不同乘法后端对比
    println!("\n3. 乘法后端对比:");
    let x = Mpz::from_str("12345678901234567890", 10)?;
    let y = Mpz::from_str("98765432109876543210", 10)?;

    println!("x = {}", x);
    println!("y = {}", y);

    // 基础乘法
    let result_basic = x.mul_with_backend(&y, MultiplicationBackend::Schoolbook);
    println!("基础乘法结果: {}", result_basic);

    // Karatsuba乘法
    let result_karatsuba = x.mul_with_backend(&y, MultiplicationBackend::Karatsuba);
    println!("Karatsuba乘法结果: {}", result_karatsuba);

    // 自适应乘法
    let result_adaptive = x.mul_with_backend(&y, MultiplicationBackend::Adaptive);
    println!("自适应乘法结果: {}", result_adaptive);

    // 验证结果一致性
    assert_eq!(result_basic, result_karatsuba);
    assert_eq!(result_karatsuba, result_adaptive);
    println!("✓ 所有乘法算法结果一致");

    // 位运算演示
    println!("\n4. 位运算演示:");
    let num1 = Mpz::from_u64(0b1010_1100);
    let num2 = Mpz::from_u64(0b1100_0011);

    println!("num1 = {} (二进制: {})", num1, num1.to_string(2));
    println!("num2 = {} (二进制: {})", num2, num2.to_string(2));

    let and_result = num1.bitwise_and(&num2);
    println!(
        "num1 & num2 = {} (二进制: {})",
        and_result,
        and_result.to_string(2)
    );

    let or_result = num1.bitwise_or(&num2);
    println!(
        "num1 | num2 = {} (二进制: {})",
        or_result,
        or_result.to_string(2)
    );

    let xor_result = num1.bitwise_xor(&num2);
    println!(
        "num1 ^ num2 = {} (二进制: {})",
        xor_result,
        xor_result.to_string(2)
    );

    // 数论函数演示
    println!("\n5. 数论函数演示:");
    let m = Mpz::from_i64(48);
    let n = Mpz::from_i64(18);

    let gcd = m.gcd(&n);
    println!("gcd({}, {}) = {}", m, n, gcd);

    let lcm = m.lcm(&n);
    println!("lcm({}, {}) = {}", m, n, lcm);

    let sqrt_val = Mpz::from_i64(16).sqrt();
    println!("sqrt(16) = {}", sqrt_val);

    // 阶乘演示
    let factorial_10 = Mpz::factorial(10);
    println!("10! = {}", factorial_10);

    // 斐波那契数列
    let fib_20 = Mpz::fibonacci(20);
    println!("fib(20) = {}", fib_20);

    // 配置演示
    println!("\n6. 配置演示:");
    println!(
        "当前乘法后端: {:?}",
        MpzMultiplicationConfig::get_global_backend()
    );

    // 设置为Karatsuba算法
    MpzMultiplicationConfig::set_global_backend(MultiplicationBackend::Karatsuba)?;
    println!(
        "设置后的乘法后端: {:?}",
        MpzMultiplicationConfig::get_global_backend()
    );

    // 获取算法阈值
    let thresholds = MpzMultiplicationConfig::get_thresholds();
    println!("算法切换阈值:");
    println!(
        "  Schoolbook -> Karatsuba: {} limbs",
        thresholds.schoolbook_to_karatsuba
    );
    println!(
        "  Karatsuba -> Toom-Cook 3: {} limbs",
        thresholds.karatsuba_to_toom3
    );
    println!(
        "  Toom-Cook 3 -> Toom-Cook 4: {} limbs",
        thresholds.toom3_to_toom4
    );

    println!("\n=== 演示完成 ===");

    Ok(())
}
