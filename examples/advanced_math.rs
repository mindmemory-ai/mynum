//! 高级数学功能演示

use mynum::Mpz;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== MyNum 高级数学功能演示 ===\n");

    // 1. 扩展欧几里得算法
    println!("1. 扩展欧几里得算法:");
    let a = Mpz::from_i64(48);
    let b = Mpz::from_i64(18);
    let (gcd, x, y) = a.extended_gcd(&b);

    println!(
        "gcd({}, {}) = {}",
        a.to_string(10),
        b.to_string(10),
        gcd.to_string(10)
    );
    println!("x = {}, y = {}", x.to_string(10), y.to_string(10));
    println!(
        "验证: {} * {} + {} * {} = {}",
        a.to_string(10),
        x.to_string(10),
        b.to_string(10),
        y.to_string(10),
        a.mul(&x).add(&b.mul(&y)).to_string(10)
    );

    // 2. 模逆元
    println!("\n2. 模逆元:");
    let a = Mpz::from_i64(3);
    let modulus = Mpz::from_i64(11);
    let inverse = a.mod_inverse(&modulus)?;

    println!(
        "{}^(-1) mod {} = {}",
        a.to_string(10),
        modulus.to_string(10),
        inverse.to_string(10)
    );
    let product = a.mul(&inverse);
    let remainder = product.rem(&modulus)?;
    println!(
        "验证: {} * {} ≡ {} (mod {})",
        a.to_string(10),
        inverse.to_string(10),
        remainder.to_string(10),
        modulus.to_string(10)
    );

    // 3. 雅可比符号
    println!("\n3. 雅可比符号:");
    let a = Mpz::from_i64(2);
    let n = Mpz::from_i64(15);
    let jacobi = a.jacobi(&n)?;
    println!("({}/{}) = {}", a.to_string(10), n.to_string(10), jacobi);

    // 4. 勒让德符号
    println!("\n4. 勒让德符号:");
    let a = Mpz::from_i64(2);
    let p = Mpz::from_i64(7);
    let legendre = a.legendre(&p)?;
    println!("({}/{}) = {}", a.to_string(10), p.to_string(10), legendre);

    // 5. 素数测试
    println!("\n5. 素数测试:");
    let test_numbers = vec![2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 17, 19, 23, 29, 31];

    for &num in &test_numbers {
        let mpz_num = Mpz::from_i64(num);
        let is_prime = mpz_num.is_probably_prime(25);
        println!("{}: {}", num, if is_prime { "素数" } else { "合数" });
    }

    // 6. 大素数测试
    println!("\n6. 大素数测试:");
    let large_prime = Mpz::from_str("1000000007", 10)?;
    println!("1000000007 是素数: {}", large_prime.is_probably_prime(25));

    let large_composite = Mpz::from_str("1000000008", 10)?;
    println!(
        "1000000008 是素数: {}",
        large_composite.is_probably_prime(25)
    );

    // 7. 随机素数生成
    println!("\n7. 随机素数生成:");
    println!("正在生成16位随机素数...");
    let random_prime = Mpz::cryptographically_secure_prime(16)?;
    println!(
        "16位随机素数: {} (二进制: {})",
        random_prime.to_string(10),
        random_prime.to_string(2)
    );
    println!("验证是素数: {}", random_prime.is_probably_prime(25));

    // 8. 数学函数
    println!("\n8. 数学函数:");

    // 阶乘
    let factorial_10 = Mpz::factorial(10);
    println!("10! = {}", factorial_10.to_string(10));

    // 斐波那契数
    let fib_20 = Mpz::fibonacci(20);
    println!("F(20) = {}", fib_20.to_string(10));

    // 卢卡斯数
    let lucas_10 = Mpz::lucas(10);
    println!("L(10) = {}", lucas_10.to_string(10));

    // 二项式系数
    let binomial = Mpz::binomial(10, 5);
    println!("C(10,5) = {}", binomial.to_string(10));

    // 9. 平方根和立方根
    println!("\n9. 平方根和立方根:");
    let perfect_square = Mpz::from_i64(256);
    let sqrt_result = perfect_square.sqrt();
    println!(
        "sqrt({}) = {}",
        perfect_square.to_string(10),
        sqrt_result.to_string(10)
    );

    let perfect_cube = Mpz::from_i64(125);
    let cbrt_result = perfect_cube.cbrt();
    println!(
        "cbrt({}) = {}",
        perfect_cube.to_string(10),
        cbrt_result.to_string(10)
    );

    // 10. 性能测试
    println!("\n10. 性能测试:");
    use std::time::Instant;

    // 大数乘法性能
    let large_a = Mpz::random_bits(512)?;
    let large_b = Mpz::random_bits(512)?;

    let start = Instant::now();
    let _result = large_a.mul(&large_b);
    let duration = start.elapsed();
    println!("512位大数乘法耗时: {:?}", duration);

    // 素数测试性能
    let test_prime = Mpz::cryptographically_secure_prime(64)?;
    let start = Instant::now();
    let _is_prime = test_prime.is_probably_prime(25);
    let duration = start.elapsed();
    println!("64位素数测试耗时: {:?}", duration);

    println!("\n=== 高级数学功能演示完成 ===");

    Ok(())
}
