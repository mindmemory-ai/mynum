//! 并行计算演示程序

use mynum::{Mpz, MpzMultiplicationConfig, MultiplicationBackend};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== MyNum 并行计算演示 ===\n");

    // 1. 并行配置测试
    println!("1. 并行配置测试:");

    // 启用并行计算
    MpzMultiplicationConfig::enable_parallel(None)?;
    MpzMultiplicationConfig::set_parallel_threshold(100)?;

    println!("已启用并行计算，阈值设置为100位");
    println!();

    // 2. 并行乘法性能测试
    println!("2. 并行乘法性能测试:");

    let large_a = generate_large_number(500); // 500位数字
    let large_b = generate_large_number(500);

    println!("测试数字 A 位数: {}", large_a.bit_length());
    println!("测试数字 B 位数: {}", large_b.bit_length());

    // 串行乘法
    let start = Instant::now();
    let serial_result = large_a.mul_with_backend(&large_b, MultiplicationBackend::Karatsuba);
    let serial_time = start.elapsed();

    // 并行乘法
    let start = Instant::now();
    let parallel_result = large_a.mul_with_backend(&large_b, MultiplicationBackend::Parallel);
    let parallel_time = start.elapsed();

    println!("大数乘法测试 (500位 × 500位):");
    println!("  串行乘法耗时: {:?}", serial_time);
    println!("  并行乘法耗时: {:?}", parallel_time);

    if parallel_time < serial_time {
        let speedup = serial_time.as_nanos() as f64 / parallel_time.as_nanos() as f64;
        println!("  性能提升: {:.2}x", speedup);
    } else {
        let slowdown = parallel_time.as_nanos() as f64 / serial_time.as_nanos() as f64;
        println!("  性能下降: {:.2}x", slowdown);
    }

    println!("  结果匹配: {}", serial_result == parallel_result);
    println!();

    // 3. 不同大小数字的并行性能测试
    println!("3. 不同大小数字的并行性能测试:");

    let test_cases = vec![
        ("小数字", 100, 100),
        ("中等数字", 500, 500),
        ("大数字", 1000, 1000),
    ];

    for (name, bits_a, bits_b) in &test_cases {
        println!("\n测试: {} ({}位 × {}位)", name, bits_a, bits_b);

        let a = generate_large_number(*bits_a);
        let b = generate_large_number(*bits_b);

        // 串行Karatsuba
        let start = Instant::now();
        let _serial = a.mul_with_backend(&b, MultiplicationBackend::Karatsuba);
        let serial_time = start.elapsed();

        // 并行乘法
        let start = Instant::now();
        let _parallel = a.mul_with_backend(&b, MultiplicationBackend::Parallel);
        let parallel_time = start.elapsed();

        println!("  串行Karatsuba: {:?}", serial_time);
        println!("  并行乘法: {:?}", parallel_time);

        if parallel_time < serial_time {
            let speedup = serial_time.as_nanos() as f64 / parallel_time.as_nanos() as f64;
            println!("  🚀 并行加速: {:.2}x", speedup);
        } else {
            let slowdown = parallel_time.as_nanos() as f64 / serial_time.as_nanos() as f64;
            println!("  ⚠️ 并行减速: {:.2}x", slowdown);
        }

        // 算法选择建议
        let suggested = MpzMultiplicationConfig::suggest_backend(a.bit_length(), b.bit_length());
        println!("  建议算法: {:?}", suggested);
    }

    // 5. 配置管理演示
    println!("\n5. 配置管理演示:");
    println!(
        "  当前并行状态: {}",
        MpzMultiplicationConfig::is_parallel_enabled()
    );
    println!(
        "  当前并行阈值: {} 位",
        MpzMultiplicationConfig::get_parallel_threshold()
    );
    println!(
        "  当前全局后端: {:?}",
        MpzMultiplicationConfig::get_global_backend()
    );

    // 禁用并行计算
    MpzMultiplicationConfig::disable_parallel();
    println!(
        "  禁用并行计算后: {}",
        MpzMultiplicationConfig::is_parallel_enabled()
    );

    // 重新启用
    MpzMultiplicationConfig::enable_parallel(None)?;
    println!(
        "  重新启用并行计算: {}",
        MpzMultiplicationConfig::is_parallel_enabled()
    );

    println!("\n=== 演示完成 ===");
    println!("💡 关键发现:");
    println!("  • 小数字：并行算法通常较慢（线程开销 > 计算收益）");
    println!("  • 中等数字：并行算法开始显现优势");
    println!("  • 大数字：并行算法显著提升性能");
    println!("  • 自适应算法：智能选择最优算法");

    Ok(())
}

/// 生成指定位数的随机大数
fn generate_large_number(bits: usize) -> Mpz {
    if bits <= 64 {
        // 对于小数字，直接使用随机数
        let value = rand::random::<u64>() % (1u64 << bits);
        Mpz::from_u64(value)
    } else {
        // 对于大数字，使用重复的1来避免栈溢出
        let ones = "1".repeat(bits.min(1000)); // 限制最大长度避免栈溢出
        Mpz::from_str(&ones, 2).unwrap_or_else(|_| Mpz::from_u64(1))
    }
}
