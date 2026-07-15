//! 基本功能演示
//!
//! 展示 MyNum 库的核心功能：大整数运算、高精度浮点数和复数运算

use mynum::{Complex, Mpf, Mpz};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== MyNum 高精度数值计算库演示 ===\n");

    // 1. 大整数运算演示
    println!("1. 大整数运算演示");
    println!("   -----------------");

    let a = Mpz::from_i64(123456789);
    let b = Mpz::from_i64(987654321);

    println!("   a = {}", a);
    println!("   b = {}", b);
    println!("   a + b = {}", a.add(&b));
    println!("   a * b = {}", a.mul(&b));
    println!("   a / b = {}", a.div(&b).unwrap());
    println!("   a % b = {}", a.rem(&b).unwrap());

    // 大整数运算
    let large_num = Mpz::from_str("123456789012345678901234567890", 10).unwrap();
    println!("   大整数: {}", large_num);
    println!("   位数: {}", large_num.bit_length());
    println!();

    // 2. 高精度浮点数演示
    println!("2. 高精度浮点数演示");
    println!("   -----------------");

    let pi = Mpf::from_f64(std::f64::consts::PI, 100);
    let e = Mpf::from_f64(std::f64::consts::E, 100);

    println!("   π ≈ {}", pi);
    println!("   e ≈ {}", e);
    println!("   π * e ≈ {}", pi.mul(&e));
    println!("   π / e ≈ {}", pi.div(&e).unwrap());

    // 高精度计算
    let precise_pi = Mpf::from_str("3.141592653589793238462643383279", 10).unwrap();
    println!("   高精度 π: {}", precise_pi);
    println!();

    // 3. 复数运算演示
    println!("3. 复数运算演示");
    println!("   -----------------");

    let z1 = Complex::from_real_imag(Mpf::from_f64(1.0, 64), Mpf::from_f64(2.0, 64));
    let z2 = Complex::from_real_imag(Mpf::from_f64(3.0, 64), Mpf::from_f64(4.0, 64));

    println!("   z1 = {}", z1);
    println!("   z2 = {}", z2);
    println!("   z1 + z2 = {}", z1.add(&z2).unwrap());
    println!("   z1 * z2 = {}", z1.mul(&z2).unwrap());

    // 复数的模和幅角
    match z1.magnitude() {
        Ok(magnitude) => println!("   |z1| = {}", magnitude),
        Err(_) => println!("   |z1| = 计算错误"),
    }
    match z1.argument() {
        Ok(argument) => println!("   arg(z1) = {}", argument),
        Err(_) => println!("   arg(z1) = 计算错误"),
    }
    println!();

    // 4. 数论函数演示
    println!("4. 数论函数演示");
    println!("   -----------------");

    let n = Mpz::from_i64(1234567);
    println!("   数字: {}", n);
    println!("   是否为素数: {}", n.is_probably_prime(25));
    println!("   平方根: {}", n.sqrt());

    let m = Mpz::from_i64(7654321);
    let gcd = Mpz::gcd(&n, &m);
    let lcm = Mpz::lcm(&n, &m);
    println!("   另一个数字: {}", m);
    println!("   GCD({}, {}) = {}", n, m, gcd);
    println!("   LCM({}, {}) = {}", n, m, lcm);
    println!();

    // 5. 性能演示
    println!("5. 性能演示");
    println!("   -----------------");

    let start = std::time::Instant::now();
    let factorial_100 = Mpz::factorial(100);
    let elapsed = start.elapsed();
    println!("   100! 计算时间: {:?}", elapsed);
    println!("   100! 位数: {}", factorial_100.bit_length());

    let start = std::time::Instant::now();
    let fib_100 = Mpz::fibonacci(100);
    let elapsed = start.elapsed();
    println!("   第100个斐波那契数计算时间: {:?}", elapsed);
    println!("   第100个斐波那契数位数: {}", fib_100.bit_length());
    println!();

    println!("=== 演示完成 ===");
    println!("MyNum 库提供了高精度数值计算的所有核心功能！");

    Ok(())
}
