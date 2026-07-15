//! 并行FFT乘法演示
//!
//! 展示并行FFT大数乘法的性能和功能

use mynum::mpz::{ParallelConfig, ParallelMultiplier};
use mynum::{Mpz, MpzMultiplicationConfig, MultiplicationBackend};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🚀 MyNum 并行FFT乘法演示");
    println!("{}", "=".repeat(60));

    // 1. 测试并行FFT乘法
    println!("\n📊 1. 并行FFT乘法测试");
    println!("{}", "-".repeat(40));

    // 创建大数进行测试
    let a = Mpz::from_str("12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890", 10)?;
    let b = Mpz::from_str("98765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210", 10)?;

    println!("测试数 A 位数: {}", a.bit_length());
    println!("测试数 B 位数: {}", b.bit_length());

    // 启用并行计算
    MpzMultiplicationConfig::enable_parallel(None)?;
    MpzMultiplicationConfig::set_parallel_threshold(100)?;

    // 测试不同算法
    let start = Instant::now();
    let result_fft = a.mul_with_backend(&b, MultiplicationBackend::FFT);
    let duration_fft = start.elapsed();
    println!("   FFT (串行): {:?}", duration_fft);

    let start = Instant::now();
    let _result_parallel_fft = a.mul_with_backend(&b, MultiplicationBackend::ParallelFFT);
    let duration_parallel_fft = start.elapsed();
    println!("   Parallel FFT: {:?}", duration_parallel_fft);

    if duration_parallel_fft < duration_fft {
        let speedup = duration_fft.as_nanos() as f64 / duration_parallel_fft.as_nanos() as f64;
        println!("   🚀 并行FFT加速比: {:.2}x", speedup);
    } else {
        let slowdown = duration_parallel_fft.as_nanos() as f64 / duration_fft.as_nanos() as f64;
        println!("   ⚠️ 并行FFT减速: {:.2}x", slowdown);
    }

    // 2. 测试直接使用并行乘法器
    println!("\n⚡ 2. 直接并行乘法器测试");
    println!("{}", "-".repeat(40));

    let parallel_config = ParallelConfig {
        num_threads: 4,
        min_task_size: 100,
        enabled: true,
    };

    let parallel_multiplier = ParallelMultiplier::new(parallel_config);

    let start = Instant::now();
    let result_direct_fft = parallel_multiplier.parallel_fft_multiply(&a, &b)?;
    let duration_direct_fft = start.elapsed();
    println!("   直接并行FFT: {:?}", duration_direct_fft);

    // 验证结果一致性
    if result_fft == result_direct_fft {
        println!("   ✅ 直接并行FFT结果与串行FFT结果完全一致");
    } else {
        println!("   ❌ 结果不一致！");
        return Err("并行FFT计算结果错误".into());
    }

    // 3. 测试不同大小的数字
    println!("\n🔢 3. 不同大小数字测试");
    println!("{}", "-".repeat(40));

    let test_cases = vec![
        (
            "中等数字 (512位)",
            Mpz::from_str(&"1".repeat(512), 2)?,
            Mpz::from_str(&"1".repeat(512), 2)?,
        ),
        (
            "大数字 (1024位)",
            Mpz::from_str(&"1".repeat(1024), 2)?,
            Mpz::from_str(&"1".repeat(1024), 2)?,
        ),
        (
            "超大数字 (2048位)",
            Mpz::from_str(&"1".repeat(2048), 2)?,
            Mpz::from_str(&"1".repeat(2048), 2)?,
        ),
    ];

    for (name, a, b) in test_cases {
        println!("\n   📊 测试: {}", name);
        println!("   A 位数: {}, B 位数: {}", a.bit_length(), b.bit_length());

        let start = Instant::now();
        let result = parallel_multiplier.parallel_fft_multiply(&a, &b)?;
        let duration = start.elapsed();
        println!("   并行FFT耗时: {:?}", duration);
        println!("   结果位数: {}", result.bit_length());
    }

    // 4. 测试分块FFT策略
    println!("\n🧩 4. 分块FFT策略测试");
    println!("{}", "-".repeat(40));

    // 创建超大数字来测试分块策略
    let huge_a = Mpz::from_str(&"1".repeat(16384), 2)?;
    let huge_b = Mpz::from_str(&"1".repeat(16384), 2)?;

    println!("   超大数 A 位数: {}", huge_a.bit_length());
    println!("   超大数 B 位数: {}", huge_b.bit_length());

    let start = Instant::now();
    let huge_result = parallel_multiplier.parallel_fft_multiply(&huge_a, &huge_b)?;
    let huge_duration = start.elapsed();
    println!("   分块FFT耗时: {:?}", huge_duration);
    println!("   结果位数: {}", huge_result.bit_length());

    // 5. 性能对比总结
    println!("\n📈 5. 性能对比总结");
    println!("{}", "-".repeat(40));

    let algorithms = [
        ("FFT (串行)", duration_fft),
        ("Parallel FFT", duration_parallel_fft),
        ("Direct Parallel FFT", duration_direct_fft),
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

    // 6. 配置建议
    println!("\n💡 6. 配置建议");
    println!("{}", "-".repeat(40));
    println!("   • 小数字 (< 1000位): 使用串行算法");
    println!("   • 中等数字 (1000-4000位): 使用并行Karatsuba");
    println!("   • 大数字 (4000-16000位): 使用并行FFT");
    println!("   • 超大数字 (> 16000位): 使用分块并行FFT");
    println!("   • 线程数: 建议设置为CPU核心数");
    println!("   • 最小任务大小: 根据硬件性能调整");

    println!("\n🎉 并行FFT乘法演示完成！");
    println!("{}", "=".repeat(60));

    Ok(())
}
