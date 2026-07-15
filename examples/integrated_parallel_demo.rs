//! 集成并行乘法演示
//!
//! 展示如何通过配置启用并行乘法，以及统一接口的使用

use mynum::{Mpz, MpzMultiplicationConfig, MultiplicationBackend};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🚀 MyNum 集成并行乘法演示");
    println!("{}", "=".repeat(50));

    // 创建两个中等大小的整数进行测试
    let a = Mpz::from_str("123456789012345678901234567890", 10)?;
    let b = Mpz::from_str("987654321098765432109876543210", 10)?;

    println!("测试数 A: {}", a);
    println!("测试数 B: {}", b);
    println!("A 的位数: {}", a.bit_length());
    println!("B 的位数: {}", b.bit_length());
    println!();

    // 1. 默认配置（自适应算法）
    println!("📊 1. 默认配置（自适应算法）");
    let start = std::time::Instant::now();
    let result_default = a.mul(&b);
    let duration_default = start.elapsed();
    println!("   结果: {}", result_default);
    println!("   耗时: {:?}", duration_default);
    println!();

    // 2. 启用并行计算
    println!("⚡ 2. 启用并行计算");
    MpzMultiplicationConfig::enable_parallel(None)?;
    MpzMultiplicationConfig::set_parallel_threshold(100)?; // 设置较低的阈值以便测试

    let start = std::time::Instant::now();
    let result_parallel = a.mul(&b);
    let duration_parallel = start.elapsed();
    println!("   结果: {}", result_parallel);
    println!("   耗时: {:?}", duration_parallel);
    println!(
        "   性能提升: {:.2}x",
        duration_default.as_nanos() as f64 / duration_parallel.as_nanos() as f64
    );
    println!();

    // 3. 验证结果一致性
    println!("✅ 3. 结果验证");
    if result_default == result_parallel {
        println!("   ✅ 并行结果与串行结果完全一致");
    } else {
        println!("   ❌ 结果不一致！");
        return Err("并行计算结果错误".into());
    }
    println!();

    // 4. 测试不同乘法后端
    println!("🔧 4. 测试不同乘法后端");

    // Schoolbook 算法
    let start = std::time::Instant::now();
    let _result_schoolbook = a.mul_with_backend(&b, MultiplicationBackend::Schoolbook);
    let duration_schoolbook = start.elapsed();
    println!("   Schoolbook: {:?}", duration_schoolbook);

    // Karatsuba 算法
    let start = std::time::Instant::now();
    let _result_karatsuba = a.mul_with_backend(&b, MultiplicationBackend::Karatsuba);
    let duration_karatsuba = start.elapsed();
    println!("   Karatsuba: {:?}", duration_karatsuba);

    // 并行算法
    let start = std::time::Instant::now();
    let _result_parallel_backend = a.mul_with_backend(&b, MultiplicationBackend::Parallel);
    let duration_parallel_backend = start.elapsed();
    println!("   Parallel: {:?}", duration_parallel_backend);

    // 自适应算法
    let start = std::time::Instant::now();
    let _result_adaptive = a.mul_with_backend(&b, MultiplicationBackend::Adaptive);
    let duration_adaptive = start.elapsed();
    println!("   Adaptive: {:?}", duration_adaptive);
    println!();

    // 5. 性能对比
    println!("📈 5. 性能对比分析");
    let algorithms = [
        ("Schoolbook", duration_schoolbook),
        ("Karatsuba", duration_karatsuba),
        ("Parallel", duration_parallel_backend),
        ("Adaptive", duration_adaptive),
    ];

    let fastest = algorithms
        .iter()
        .min_by_key(|(_, duration)| duration)
        .unwrap();
    println!("   最快算法: {} ({:?})", fastest.0, fastest.1);

    for (name, duration) in algorithms {
        let ratio = duration.as_nanos() as f64 / fastest.1.as_nanos() as f64;
        println!("   {}: {:.2}x", name, ratio);
    }
    println!();

    // 6. 大数测试
    println!("🔢 6. 大数测试");
    let big_a = Mpz::from_str(
        "123456789012345678901234567890123456789012345678901234567890",
        10,
    )?;
    let big_b = Mpz::from_str(
        "987654321098765432109876543210987654321098765432109876543210",
        10,
    )?;

    println!("   大数 A 位数: {}", big_a.bit_length());
    println!("   大数 B 位数: {}", big_b.bit_length());

    let start = std::time::Instant::now();
    let big_result = big_a.mul(&b);
    let big_duration = start.elapsed();
    println!("   大数乘法结果位数: {}", big_result.bit_length());
    println!("   大数乘法耗时: {:?}", big_duration);
    println!();

    // 7. 配置管理演示
    println!("⚙️ 7. 配置管理演示");
    println!(
        "   当前并行状态: {}",
        MpzMultiplicationConfig::is_parallel_enabled()
    );
    println!(
        "   当前并行阈值: {} 位",
        MpzMultiplicationConfig::get_parallel_threshold()
    );
    println!(
        "   当前全局后端: {:?}",
        MpzMultiplicationConfig::get_global_backend()
    );

    // 禁用并行计算
    MpzMultiplicationConfig::disable_parallel();
    println!(
        "   禁用并行计算后: {}",
        MpzMultiplicationConfig::is_parallel_enabled()
    );

    // 重新启用
    MpzMultiplicationConfig::enable_parallel(None)?;
    println!(
        "   重新启用并行计算: {}",
        MpzMultiplicationConfig::is_parallel_enabled()
    );
    println!();

    // 8. 错误处理测试
    println!("⚠️ 8. 错误处理测试");
    let tiny_a = Mpz::from_i64(123);
    let tiny_b = Mpz::from_i64(456);

    // 即使启用了并行计算，小数仍然使用串行算法
    let start = std::time::Instant::now();
    let tiny_result = tiny_a.mul(&tiny_b);
    let tiny_duration = start.elapsed();
    println!("   小数乘法: {} * {} = {}", tiny_a, tiny_b, tiny_result);
    println!("   小数乘法耗时: {:?}", tiny_duration);
    println!("   注意：小数自动使用串行算法，不触发并行计算");
    println!();

    println!("🎉 演示完成！");
    println!("{}", "=".repeat(50));
    println!("💡 关键特性:");
    println!("   • 统一的乘法接口 (a.mul(&b))");
    println!("   • 自动算法选择（基于配置和操作数大小）");
    println!("   • 并行计算作为后端无缝集成");
    println!("   • 智能阈值控制（避免小数的并行开销）");
    println!("   • 失败时自动回退到串行算法");

    Ok(())
}
