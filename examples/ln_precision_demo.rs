use mynum::{Mpf, Mpz};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== 高精度ln函数演示 ===\n");

    // 测试不同精度下的ln函数
    let precisions = [64, 128, 256];

    for &precision in &precisions {
        println!("精度: {} bits", precision);
        println!("{}", "-".repeat(50));

        // 测试ln(2)
        let two = Mpf::from_mpz(Mpz::from_i64(2), precision);
        let ln2 = two.ln()?;
        println!("ln(2) = {}", ln2.to_string(10));

        // 测试ln(e) ≈ 1
        let e = Mpf::from_str("2.718281828459045", 10)?;
        let ln_e = e.ln()?;
        println!("ln(e) = {}", ln_e.to_string(10));

        // 测试ln(10)
        let ten = Mpf::from_mpz(Mpz::from_i64(10), precision);
        let ln10 = ten.ln()?;
        println!("ln(10) = {}", ln10.to_string(10));

        // 测试ln(0.5) = -ln(2)
        let half = Mpf::from_mpz(Mpz::from_i64(1), precision)
            .div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?;
        let ln_half = half.ln()?;
        println!("ln(0.5) = {}", ln_half.to_string(10));

        // 测试ln(100) = 2*ln(10)
        let hundred = Mpf::from_mpz(Mpz::from_i64(100), precision);
        let ln100 = hundred.ln()?;
        println!("ln(100) = {}", ln100.to_string(10));

        println!();
    }

    // 验证对数恒等式
    println!("=== 验证对数恒等式 ===\n");

    let precision = 128;
    let a = Mpf::from_mpz(Mpz::from_i64(3), precision);
    let b = Mpf::from_mpz(Mpz::from_i64(4), precision);

    // ln(a*b) = ln(a) + ln(b)
    let product = a.mul(&b);
    let ln_product = product.ln()?;
    let ln_a = a.ln()?;
    let ln_b = b.ln()?;
    let ln_sum = ln_a.add(&ln_b);

    println!("验证: ln(a*b) = ln(a) + ln(b)");
    println!("a = {}, b = {}", a.to_string(10), b.to_string(10));
    println!(
        "ln(a*b) = ln({}) = {}",
        product.to_string(10),
        ln_product.to_string(10)
    );
    println!(
        "ln(a) + ln(b) = {} + {} = {}",
        ln_a.to_string(10),
        ln_b.to_string(10),
        ln_sum.to_string(10)
    );

    let diff = ln_product.sub(&ln_sum).abs();
    println!("差值: {}", diff.to_string(10));
    println!();

    // ln(a^n) = n*ln(a)
    let n = 3;
    let a_pow_n = a.pow(n as u32)?;
    let ln_a_pow_n = a_pow_n.ln()?;
    let n_ln_a = ln_a.mul(&Mpf::from_mpz(Mpz::from_i64(n), precision));

    println!("验证: ln(a^n) = n*ln(a)");
    println!("a = {}, n = {}", a.to_string(10), n);
    println!(
        "ln(a^n) = ln({}) = {}",
        a_pow_n.to_string(10),
        ln_a_pow_n.to_string(10)
    );
    println!(
        "n*ln(a) = {}*{} = {}",
        n,
        ln_a.to_string(10),
        n_ln_a.to_string(10)
    );

    let diff2 = ln_a_pow_n.sub(&n_ln_a).abs();
    println!("差值: {}", diff2.to_string(10));
    println!();

    // 性能测试
    println!("=== 性能测试 ===\n");

    let start = std::time::Instant::now();
    let test_values = [1.1, 1.5, 2.0, std::f64::consts::PI, 10.0, 100.0, 1000.0];

    for &val in &test_values {
        let x = Mpf::from_str(&val.to_string(), 10)?;
        let ln_x = x.ln()?;
        println!("ln({}) = {}", val, ln_x.to_string(10));
    }

    let duration = start.elapsed();
    println!("\n计算 {} 个ln值耗时: {:?}", test_values.len(), duration);

    Ok(())
}
