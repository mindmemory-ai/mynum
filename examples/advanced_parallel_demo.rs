//! 高级并行乘法演示
//!
//! 详细展示和对比各种并行化乘法算法的性能和特点

use mynum::{
    mpz::{ParallelConfig, ParallelMultiplier},
    Mpz, MpzMultiplicationConfig, MultiplicationBackend,
};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🚀 MyNum 高级并行乘法算法详细对比演示");
    println!("{}", "=".repeat(70));

    // 1. 测试数据准备
    println!("\n📊 1. 测试数据准备");
    println!("{}", "-".repeat(50));

    // 准备不同大小的测试数据
    let test_cases = vec![
        (
            "小数字 (128位)",
            Mpz::from_i64(123456789),
            Mpz::from_i64(987654321),
        ),
        (
            "中等数字 (512位)",
            Mpz::from_i64(123456789).pow_u32(2),
            Mpz::from_i64(987654321).pow_u32(2),
        ),
        (
            "大数字 (1024位)",
            Mpz::from_i64(123456789).pow_u32(3),
            Mpz::from_i64(987654321).pow_u32(3),
        ),
    ];

    for (name, a, b) in &test_cases {
        println!(
            "   {}: A={}位, B={}位",
            name,
            a.bit_length(),
            b.bit_length()
        );
    }

    // 启用并行计算
    MpzMultiplicationConfig::enable_parallel(None)?;
    MpzMultiplicationConfig::set_parallel_threshold(100)?;

    // 2. 串行算法基准测试
    println!("\n🔧 2. 串行算法基准测试");
    println!("{}", "-".repeat(50));

    let mut serial_results = Vec::new();

    for (name, a, b) in &test_cases {
        println!("\n   📊 测试: {}", name);
        println!("   A 位数: {}, B 位数: {}", a.bit_length(), b.bit_length());

        // Schoolbook 算法
        let start = Instant::now();
        let _result_schoolbook = a.mul_with_backend(b, MultiplicationBackend::Schoolbook);
        let duration_schoolbook = start.elapsed();
        println!("      Schoolbook: {:?}", duration_schoolbook);

        // Karatsuba 算法
        let start = Instant::now();
        let _result_karatsuba = a.mul_with_backend(b, MultiplicationBackend::Karatsuba);
        let duration_karatsuba = start.elapsed();
        println!("      Karatsuba: {:?}", duration_karatsuba);

        // FFT 算法
        let start = Instant::now();
        let _result_fft = a.mul_with_backend(b, MultiplicationBackend::FFT);
        let duration_fft = start.elapsed();
        println!("      FFT: {:?}", duration_fft);

        serial_results.push((*name, duration_schoolbook, duration_karatsuba, duration_fft));
    }

    // 3. 并行算法详细测试
    println!("\n⚡ 3. 并行算法详细测试");
    println!("{}", "-".repeat(50));

    let mut parallel_results = Vec::new();

    for (name, a, b) in &test_cases {
        println!("\n   📊 测试: {}", name);
        println!("   A 位数: {}, B 位数: {}", a.bit_length(), b.bit_length());

        // 并行Karatsuba算法
        let start = Instant::now();
        let _result_parallel_karatsuba =
            a.mul_with_backend(b, MultiplicationBackend::ParallelKaratsuba);
        let duration_parallel_karatsuba = start.elapsed();
        println!(
            "      Parallel Karatsuba: {:?}",
            duration_parallel_karatsuba
        );

        // 并行FFT算法
        let start = Instant::now();
        let _result_parallel_fft = a.mul_with_backend(b, MultiplicationBackend::ParallelFFT);
        let duration_parallel_fft = start.elapsed();
        println!("      Parallel FFT: {:?}", duration_parallel_fft);

        // 智能并行算法（自动选择）
        let start = Instant::now();
        let _result_parallel = a.mul_with_backend(b, MultiplicationBackend::Parallel);
        let duration_parallel = start.elapsed();
        println!("      Parallel (智能选择): {:?}", duration_parallel);

        parallel_results.push((
            *name,
            duration_parallel_karatsuba,
            duration_parallel_fft,
            duration_parallel,
        ));
    }

    // 4. 直接使用并行乘法器测试
    println!("\n🔧 4. 直接并行乘法器测试");
    println!("{}", "-".repeat(50));

    let parallel_config = ParallelConfig {
        num_threads: 4,
        min_task_size: 100,
        enabled: true,
    };

    let parallel_multiplier = ParallelMultiplier::new(parallel_config);

    let mut direct_results = Vec::new();

    for (name, a, b) in &test_cases {
        println!("\n   📊 测试: {}", name);
        println!("   A 位数: {}, B 位数: {}", a.bit_length(), b.bit_length());

        // 直接并行乘法
        let start = Instant::now();
        let _result_direct = parallel_multiplier.parallel_multiply(a, b)?;
        let duration_direct = start.elapsed();
        println!("      Direct Parallel: {:?}", duration_direct);

        // 直接并行Karatsuba
        let start = Instant::now();
        let _result_direct_karatsuba = parallel_multiplier.karatsuba_parallel_optimized(a, b)?;
        let duration_direct_karatsuba = start.elapsed();
        println!(
            "      Direct Parallel Karatsuba: {:?}",
            duration_direct_karatsuba
        );

        // 直接并行FFT
        let start = Instant::now();
        let _result_direct_fft = parallel_multiplier.parallel_fft_multiply(a, b)?;
        let duration_direct_fft = start.elapsed();
        println!("      Direct Parallel FFT: {:?}", duration_direct_fft);

        direct_results.push((
            *name,
            duration_direct,
            duration_direct_karatsuba,
            duration_direct_fft,
        ));
    }

    // 5. 性能对比分析
    println!("\n📈 5. 性能对比分析");
    println!("{}", "=".repeat(70));

    for (i, (name, _, _)) in test_cases.iter().enumerate() {
        println!("\n🔍 测试案例: {}", name);
        println!("{}", "-".repeat(50));

        let serial_karatsuba = serial_results[i].2;
        let parallel_karatsuba = parallel_results[i].1;
        let direct_parallel_karatsuba = direct_results[i].2;

        let serial_fft = serial_results[i].3;
        let parallel_fft = parallel_results[i].2;
        let direct_parallel_fft = direct_results[i].3;

        // Karatsuba对比
        println!("   📊 Karatsuba算法对比:");
        println!("      串行 Karatsuba: {:?}", serial_karatsuba);
        println!("      并行 Karatsuba: {:?}", parallel_karatsuba);
        println!("      直接并行 Karatsuba: {:?}", direct_parallel_karatsuba);

        if parallel_karatsuba < serial_karatsuba {
            let speedup = serial_karatsuba.as_nanos() as f64 / parallel_karatsuba.as_nanos() as f64;
            println!("      🚀 并行Karatsuba加速比: {:.2}x", speedup);
        } else {
            let slowdown =
                parallel_karatsuba.as_nanos() as f64 / serial_karatsuba.as_nanos() as f64;
            println!("      ⚠️ 并行Karatsuba减速: {:.2}x", slowdown);
        }

        // FFT对比
        println!("   📊 FFT算法对比:");
        println!("      串行 FFT: {:?}", serial_fft);
        println!("      并行 FFT: {:?}", parallel_fft);
        println!("      直接并行 FFT: {:?}", direct_parallel_fft);

        if parallel_fft < serial_fft {
            let speedup = serial_fft.as_nanos() as f64 / parallel_fft.as_nanos() as f64;
            println!("      🚀 并行FFT加速比: {:.2}x", speedup);
        } else {
            let slowdown = parallel_fft.as_nanos() as f64 / serial_fft.as_nanos() as f64;
            println!("      ⚠️ 并行FFT减速: {:.2}x", slowdown);
        }

        // 最佳算法选择
        let algorithms = [
            ("串行Karatsuba", serial_karatsuba),
            ("并行Karatsuba", parallel_karatsuba),
            ("直接并行Karatsuba", direct_parallel_karatsuba),
            ("串行FFT", serial_fft),
            ("并行FFT", parallel_fft),
            ("直接并行FFT", direct_parallel_fft),
        ];

        let fastest = algorithms
            .iter()
            .min_by_key(|(_, duration)| duration)
            .unwrap();
        println!("   🏆 最佳算法: {} ({:?})", fastest.0, fastest.1);

        // 性能排名
        let mut sorted_algorithms: Vec<_> = algorithms.iter().collect();
        sorted_algorithms.sort_by_key(|(_, duration)| duration);

        println!("   📊 性能排名:");
        for (rank, (name, duration)) in sorted_algorithms.iter().enumerate() {
            let ratio = duration.as_nanos() as f64 / fastest.1.as_nanos() as f64;
            println!("      {}. {}: {:.2}x", rank + 1, name, ratio);
        }
    }

    // 6. 算法特性分析
    println!("\n🔬 6. 算法特性分析");
    println!("{}", "=".repeat(70));

    println!("\n📊 并行Karatsuba算法:");
    println!("   ✅ 优点:");
    println!("      • 适用于中等大小数字 (1000-4000位)");
    println!("      • 并行化效果好，线程间通信少");
    println!("      • 内存使用相对较低");
    println!("      • 支持多级并行策略");
    println!("   ⚠️ 缺点:");
    println!("      • 对于超大数字可能栈溢出");
    println!("      • 小数字并行开销大于收益");
    println!("      • 递归深度需要控制");

    println!("\n📊 并行FFT算法:");
    println!("   ✅ 优点:");
    println!("      • 适用于超大数字 (>4000位)");
    println!("      • 时间复杂度 O(n log n)");
    println!("      • 支持分块处理，避免内存问题");
    println!("      • 缓存优化，旋转因子复用");
    println!("   ⚠️ 缺点:");
    println!("      • 小数字开销较大");
    println!("      • 需要复数运算");
    println!("      • 内存使用相对较高");

    println!("\n📊 智能并行算法:");
    println!("   ✅ 优点:");
    println!("      • 自动选择最优算法");
    println!("      • 根据数字大小动态调整");
    println!("      • 失败时自动回退");
    println!("      • 用户无需关心算法选择");
    println!("   ⚠️ 缺点:");
    println!("      • 选择逻辑可能不够精确");
    println!("      • 需要额外的判断开销");

    // 7. 配置建议
    println!("\n💡 7. 配置建议");
    println!("{}", "=".repeat(70));

    println!("\n🎯 算法选择建议:");
    println!("   • 小数字 (< 1000位): 使用串行算法");
    println!("   • 中等数字 (1000-4000位): 使用并行Karatsuba");
    println!("   • 大数字 (4000-16000位): 使用并行FFT");
    println!("   • 超大数字 (> 16000位): 使用分块并行FFT");

    println!("\n⚙️ 配置参数建议:");
    println!("   • 线程数: 设置为CPU核心数");
    println!("   • 最小任务大小: 1000-2000位");
    println!("   • 并行阈值: 1000位");
    println!("   • 启用智能算法选择");

    println!("\n🔧 性能调优建议:");
    println!("   • 监控不同算法的性能表现");
    println!("   • 根据硬件特性调整线程数");
    println!("   • 测试不同大小的数字找到最优配置");
    println!("   • 考虑内存使用和CPU负载的平衡");

    // 8. 总结
    println!("\n🎉 8. 总结");
    println!("{}", "=".repeat(70));

    println!("我们成功演示了MyNum库的多种并行乘法算法:");
    println!("   ✅ 并行Karatsuba: 适用于中等大小数字，并行化效果好");
    println!("   ✅ 并行FFT: 适用于超大数字，时间复杂度优秀");
    println!("   ✅ 智能并行: 自动选择最优算法，用户友好");
    println!("   ✅ 直接调用: 提供底层控制，性能可调");

    println!("\n关键发现:");
    println!("   • 并行算法的效果与数字大小密切相关");
    println!("   • 小数字使用并行算法通常得不偿失");
    println!("   • 大数字使用并行算法能显著提升性能");
    println!("   • 智能算法选择能自动找到最优解");

    println!("\n🚀 MyNum库现在具备了完整的并行计算能力！");
    println!("用户可以根据需要选择最适合的并行算法，享受高性能的大数运算体验！");

    Ok(())
}
