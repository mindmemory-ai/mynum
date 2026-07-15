//! 随机数生成演示

use mynum::Mpz;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== MyNum 随机数生成演示 ===\n");

    // 1. 生成指定位数的随机数
    println!("1. 指定位数随机数生成:");
    let random_64bit = Mpz::random_bits(64)?;
    let random_128bit = Mpz::random_bits(128)?;
    let random_256bit = Mpz::random_bits(256)?;

    println!("64位随机数: {}", random_64bit.to_string(16));
    println!("128位随机数: {}", random_128bit.to_string(16));
    println!("256位随机数: {}", random_256bit.to_string(16));

    // 2. 范围随机数生成
    println!("\n2. 范围随机数生成:");
    let min = Mpz::from_i64(1000);
    let max = Mpz::from_i64(10000);

    for i in 0..5 {
        let random_in_range = Mpz::random_range(&min, &max)?;
        println!("随机数 {}: {}", i + 1, random_in_range.to_string(10));
    }

    // 3. 随机奇数生成
    println!("\n3. 随机奇数生成:");
    let random_odd_32 = Mpz::random_odd_bits(32)?;
    let random_odd_64 = Mpz::random_odd_bits(64)?;

    println!(
        "32位随机奇数: {} (二进制: {})",
        random_odd_32.to_string(10),
        random_odd_32.to_string(2)
    );
    println!(
        "64位随机奇数: {} (二进制: {})",
        random_odd_64.to_string(10),
        random_odd_64.to_string(2)
    );

    // 验证确实是奇数
    assert!(random_odd_32.test_bit(0));
    assert!(random_odd_64.test_bit(0));
    println!("✓ 验证：两个数都是奇数");

    // 4. 密码学安全随机数
    println!("\n4. 密码学安全随机数:");
    let secure_random = Mpz::cryptographically_secure_random(128)?;
    println!("128位安全随机数: {}", secure_random.to_string(16));

    // 5. 随机素数生成（小位数用于演示）
    println!("\n5. 随机素数生成:");
    println!("正在生成16位随机素数...");
    let prime_16 = Mpz::cryptographically_secure_prime(16)?;
    println!(
        "16位随机素数: {} (二进制: {})",
        prime_16.to_string(10),
        prime_16.to_string(2)
    );

    // 验证确实是素数
    assert!(prime_16.is_probably_prime(25));
    println!("✓ 验证：确实是素数");

    // 6. 统计测试
    println!("\n6. 随机数统计测试:");
    let mut numbers = Vec::new();
    for _ in 0..1000 {
        numbers.push(Mpz::random_below(&Mpz::from_i64(100))?);
    }

    // 计算平均值（简化版本）
    let sum: u64 = numbers.iter().map(|n| n.to_u64().unwrap_or(0)).sum();
    let average = sum as f64 / numbers.len() as f64;

    println!("生成1000个[0,100)范围内的随机数");
    println!("平均值: {:.2} (期望值: 50.0)", average);
    println!("分布偏差: {:.2}%", ((average - 50.0) / 50.0 * 100.0).abs());

    // 7. 性能测试
    println!("\n7. 性能测试:");
    use std::time::Instant;

    let start = Instant::now();
    for _ in 0..1000 {
        let _ = Mpz::random_bits(64)?;
    }
    let duration = start.elapsed();
    println!("生成1000个64位随机数耗时: {:?}", duration);
    println!("平均每个随机数耗时: {:?}", duration / 1000);

    println!("\n=== 随机数演示完成 ===");

    Ok(())
}
