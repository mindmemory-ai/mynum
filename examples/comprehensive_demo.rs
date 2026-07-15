//! MyNum 全面功能演示
//!
//! 这个示例展示了 MyNum 库的所有主要功能，包括：
//! - 基础算术运算
//! - 高级数论函数
//! - 随机数生成
//! - 位运算
//! - 乘法后端
//! - 性能测试

use mynum::Mpz;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== MyNum 全面功能演示 ===\n");

    // 1. 基础算术运算
    println!("1. 基础算术运算:");
    let a = Mpz::from_str("12345678901234567890", 10)?;
    let b = Mpz::from_str("98765432109876543210", 10)?;

    println!("a = {}", a.to_string(10));
    println!("b = {}", b.to_string(10));
    println!("a + b = {}", a.add(&b).to_string(10));
    println!("a - b = {}", a.sub(&b).to_string(10));
    println!("a * b = {}", a.mul(&b).to_string(10));
    println!("a / b = {}", a.div(&b)?.to_string(10));
    println!("a % b = {}", a.rem(&b)?.to_string(10));
    println!("a mod b = {}", a.mod_(&b)?.to_string(10));

    // 2. 幂运算
    println!("\n2. 幂运算:");
    let base = Mpz::from_i64(2);
    let exp = 100u32;
    let start = Instant::now();
    let result = base.pow_u32(exp);
    let duration = start.elapsed();
    println!(
        "2^100 = {}...",
        result.to_string(10).chars().take(20).collect::<String>()
    );
    println!("计算耗时: {:?}", duration);

    // 3. 数论函数
    println!("\n3. 数论函数:");

    // GCD 和 LCM
    let x = Mpz::from_i64(48);
    let y = Mpz::from_i64(18);
    println!(
        "gcd({}, {}) = {}",
        x.to_string(10),
        y.to_string(10),
        x.gcd(&y).to_string(10)
    );
    println!(
        "lcm({}, {}) = {}",
        x.to_string(10),
        y.to_string(10),
        x.lcm(&y).to_string(10)
    );

    // 扩展欧几里得算法
    let (gcd, s, t) = x.extended_gcd(&y);
    println!(
        "扩展欧几里得: {} * {} + {} * {} = {}",
        x.to_string(10),
        s.to_string(10),
        y.to_string(10),
        t.to_string(10),
        gcd.to_string(10)
    );

    // 模逆元
    let a_inv = Mpz::from_i64(3);
    let modulus = Mpz::from_i64(11);
    let inverse = a_inv.mod_inverse(&modulus)?;
    println!(
        "{}^(-1) mod {} = {}",
        a_inv.to_string(10),
        modulus.to_string(10),
        inverse.to_string(10)
    );

    // 4. 素数测试和生成
    println!("\n4. 素数测试和生成:");

    // 测试小素数
    let test_primes = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31];
    for &p in &test_primes {
        let mpz_p = Mpz::from_i64(p);
        println!(
            "{}: {}",
            p,
            if mpz_p.is_probably_prime(25) {
                "素数"
            } else {
                "合数"
            }
        );
    }

    // 生成随机素数
    println!("正在生成32位随机素数...");
    let start = Instant::now();
    let random_prime = Mpz::random_prime(32)?;
    let duration = start.elapsed();
    println!(
        "32位随机素数: {} (耗时: {:?})",
        random_prime.to_string(10),
        duration
    );

    // 5. 二次剩余
    println!("\n5. 二次剩余:");
    let a = Mpz::from_i64(5);
    let p = Mpz::from_i64(11);
    let sqrt_mod = a.sqrt_mod(&p)?;
    println!(
        "sqrt({}) mod {} = {:?}",
        a.to_string(10),
        p.to_string(10),
        sqrt_mod.iter().map(|x| x.to_string(10)).collect::<Vec<_>>()
    );

    // 6. 位运算
    println!("\n6. 位运算:");
    let num1 = Mpz::from_u64(0b1010);
    let num2 = Mpz::from_u64(0b1100);

    println!(
        "{} & {} = {}",
        num1.to_string(2),
        num2.to_string(2),
        num1.bitwise_and(&num2).to_string(2)
    );
    println!(
        "{} | {} = {}",
        num1.to_string(2),
        num2.to_string(2),
        num1.bitwise_or(&num2).to_string(2)
    );
    println!(
        "{} ^ {} = {}",
        num1.to_string(2),
        num2.to_string(2),
        num1.bitwise_xor(&num2).to_string(2)
    );
    println!(
        "~{} = {}",
        num1.to_string(2),
        num1.bitwise_not().to_string(2)
    );

    // 位操作
    let mut test_num = Mpz::from_u64(0b1010);
    println!(
        "原始数: {} (二进制: {})",
        test_num.to_string(10),
        test_num.to_string(2)
    );

    test_num.set_bit(0);
    println!(
        "设置第0位: {} (二进制: {})",
        test_num.to_string(10),
        test_num.to_string(2)
    );

    test_num.clear_bit(3);
    println!(
        "清除第3位: {} (二进制: {})",
        test_num.to_string(10),
        test_num.to_string(2)
    );

    test_num.flip_bit(2);
    println!(
        "翻转第2位: {} (二进制: {})",
        test_num.to_string(10),
        test_num.to_string(2)
    );

    // 7. 乘法后端性能测试
    println!("\n7. 乘法后端性能测试:");
    let large_a = Mpz::random_bits(256)?;
    let large_b = Mpz::random_bits(256)?;

    // Schoolbook 乘法
    let start = Instant::now();
    let _result1 = large_a.mul_schoolbook(&large_b);
    let duration1 = start.elapsed();
    println!("Schoolbook 乘法 (256位): {:?}", duration1);

    // Karatsuba 乘法
    let start = Instant::now();
    let _result2 = large_a.mul_karatsuba(&large_b);
    let duration2 = start.elapsed();
    println!("Karatsuba 乘法 (256位): {:?}", duration2);

    // 自适应乘法
    let start = Instant::now();
    let _result3 = large_a.mul_with_backend(&large_b, mynum::MultiplicationBackend::Adaptive);
    let duration3 = start.elapsed();
    println!("自适应乘法 (256位): {:?}", duration3);

    // 8. 随机数生成
    println!("\n8. 随机数生成:");

    // 随机位
    let random_bits = Mpz::random_bits(64)?;
    println!(
        "64位随机数: {} (二进制: {})",
        random_bits.to_string(10),
        random_bits.to_string(2)
    );

    // 随机范围
    let min = Mpz::from_i64(1000);
    let max = Mpz::from_i64(10000);
    let random_range = Mpz::random_range(&min, &max)?;
    println!("随机数 [1000, 10000): {}", random_range.to_string(10));

    // 密码学安全随机数
    let secure_random = Mpz::cryptographically_secure_random(128)?;
    println!(
        "128位密码学安全随机数: {}...",
        secure_random
            .to_string(16)
            .chars()
            .take(32)
            .collect::<String>()
    );

    // 9. 数学函数
    println!("\n9. 数学函数:");

    // 阶乘
    let factorial_10 = Mpz::factorial(10);
    println!("10! = {}", factorial_10.to_string(10));

    // 斐波那契数
    let fib_30 = Mpz::fibonacci(30);
    println!("F(30) = {}", fib_30.to_string(10));

    // 卢卡斯数
    let lucas_20 = Mpz::lucas(20);
    println!("L(20) = {}", lucas_20.to_string(10));

    // 二项式系数
    let binomial = Mpz::binomial(20, 10);
    println!("C(20,10) = {}", binomial.to_string(10));

    // 10. 平方根和立方根
    println!("\n10. 平方根和立方根:");
    let perfect_square = Mpz::from_i64(1000000);
    let sqrt_result = perfect_square.sqrt();
    println!(
        "sqrt({}) = {}",
        perfect_square.to_string(10),
        sqrt_result.to_string(10)
    );

    let perfect_cube = Mpz::from_i64(1000000);
    let cbrt_result = perfect_cube.cbrt();
    println!(
        "cbrt({}) = {}",
        perfect_cube.to_string(10),
        cbrt_result.to_string(10)
    );

    // 11. 大数运算性能测试
    println!("\n11. 大数运算性能测试:");

    // 1024位大数运算
    let huge_a = Mpz::random_bits(1024)?;
    let huge_b = Mpz::random_bits(1024)?;

    let start = Instant::now();
    let _huge_sum = huge_a.add(&huge_b);
    let duration = start.elapsed();
    println!("1024位加法: {:?}", duration);

    let start = Instant::now();
    let _huge_product = huge_a.mul(&huge_b);
    let duration = start.elapsed();
    println!("1024位乘法: {:?}", duration);

    // 12. 内存使用情况
    println!("\n12. 内存使用情况:");
    let test_number = Mpz::random_bits(2048)?;
    println!("2048位数字的limb数量: {}", test_number.limb_count());
    println!("2048位数字的位长度: {}", test_number.bit_length());

    println!("\n=== 全面功能演示完成 ===");
    println!("所有功能测试通过！");

    Ok(())
}
