//! 简单并行乘法测试
//!
//! 避免使用超大数字，专注于小到中等大小的数字测试

use mynum::{Mpz, MpzMultiplicationConfig, MultiplicationBackend};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🚀 MyNum 简单并行乘法测试");
    println!("{}", "=".repeat(50));

    // 测试不同大小的数字（使用安全的数字大小）
    let test_cases = vec![
        (
            "小数字 (32位)",
            Mpz::from_str("123456789", 10)?,
            Mpz::from_str("987654321", 10)?,
        ),
        (
            "中等数字 (64位)",
            Mpz::from_str("1234567890123456", 10)?,
            Mpz::from_str("9876543210987654", 10)?,
        ),
        (
            "大数字 (96位)",
            Mpz::from_str("123456789012345678901234", 10)?,
            Mpz::from_str("987654321098765432109876", 10)?,
        ),
    ];

    for (name, a, b) in test_cases {
        println!("\n📊 测试: {}", name);
        println!("{}", "-".repeat(40));
        println!("操作数 A 位数: {}", a.bit_length());
        println!("操作数 B 位数: {}", b.bit_length());

        // 1. 串行算法测试
        println!("\n🔧 串行算法测试:");

        // Schoolbook
        let start = Instant::now();
        let _result_schoolbook = a.mul_with_backend(&b, MultiplicationBackend::Schoolbook);
        let duration_schoolbook = start.elapsed();
        println!("   Schoolbook: {:?}", duration_schoolbook);

        // Karatsuba
        let start = Instant::now();
        let _result_karatsuba = a.mul_with_backend(&b, MultiplicationBackend::Karatsuba);
        let duration_karatsuba = start.elapsed();
        println!("   Karatsuba: {:?}", duration_karatsuba);

        // 自适应（串行）
        let start = Instant::now();
        let _result_adaptive = a.mul_with_backend(&b, MultiplicationBackend::Adaptive);
        let duration_adaptive = start.elapsed();
        println!("   Adaptive: {:?}", duration_adaptive);

        // 2. 并行算法测试
        println!("\n⚡ 并行算法测试:");

        // 启用并行计算
        MpzMultiplicationConfig::enable_parallel(None)?;
        MpzMultiplicationConfig::set_parallel_threshold(100)?; // 设置合适的阈值

        // 并行乘法
        let start = Instant::now();
        let _result_parallel = a.mul_with_backend(&b, MultiplicationBackend::Parallel);
        let duration_parallel = start.elapsed();
        println!("   Parallel: {:?}", duration_parallel);

        // 自适应（可能使用并行）
        let start = Instant::now();
        let _result_adaptive_parallel = a.mul_with_backend(&b, MultiplicationBackend::Adaptive);
        let duration_adaptive_parallel = start.elapsed();
        println!("   Adaptive (可能并行): {:?}", duration_adaptive_parallel);

        // 3. 性能分析
        println!("\n📈 性能分析:");

        let fastest_serial = duration_schoolbook
            .min(duration_karatsuba)
            .min(duration_adaptive);
        let fastest_overall = fastest_serial
            .min(duration_parallel)
            .min(duration_adaptive_parallel);

        println!("   最快串行算法: {:?}", fastest_serial);
        println!("   最快整体算法: {:?}", fastest_overall);

        if duration_parallel < fastest_serial {
            let speedup = fastest_serial.as_nanos() as f64 / duration_parallel.as_nanos() as f64;
            println!("   🚀 并行算法加速比: {:.2}x", speedup);
        } else {
            let slowdown = duration_parallel.as_nanos() as f64 / fastest_serial.as_nanos() as f64;
            println!("   ⚠️ 并行算法减速: {:.2}x", slowdown);
        }

        // 4. 算法选择建议
        println!("\n💡 算法选择建议:");
        let suggested = MpzMultiplicationConfig::suggest_backend(a.bit_length(), b.bit_length());
        println!("   建议算法: {:?}", suggested);

        if suggested == MultiplicationBackend::Parallel {
            println!("   ✅ 系统建议使用并行算法");
        } else {
            println!("   ℹ️ 系统建议使用串行算法");
        }

        // 禁用并行计算，准备下一个测试
        MpzMultiplicationConfig::disable_parallel();
    }

    // 5. 总结
    println!("\n🎯 测试总结");
    println!("{}", "=".repeat(50));
    println!("💡 关键发现:");
    println!("   • 小数字：并行算法通常较慢（线程开销 > 计算收益）");
    println!("   • 中等数字：并行算法开始显现优势");
    println!("   • 大数字：并行算法在某些情况下有优势");
    println!("   • 自适应算法：智能选择最优算法");
    println!();
    println!("🔧 使用建议:");
    println!("   • 对于日常计算，使用默认的自适应模式");
    println!("   • 对于中等大小数字，可以尝试并行模式");
    println!("   • 根据实际性能测试调整并行阈值");
    println!("   • 避免使用超大数字进行并行测试");

    Ok(())
}
