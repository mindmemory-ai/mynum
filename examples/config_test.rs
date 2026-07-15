//! 配置系统测试例子
//!
//! 验证重构后的配置系统是否正常工作

use mynum::{
    ComplexConfig, GlobalPrecisionConfig, MpfConfig, MpzMultiplicationConfig, MultiplicationBackend,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔧 MyNum 配置系统测试");
    println!("{}", "=".repeat(50));

    // 1. 测试全局精度配置
    println!("\n📊 1. 全局精度配置测试");
    println!("{}", "-".repeat(30));

    println!(
        "默认精度: {} 位",
        GlobalPrecisionConfig::get_default_precision()
    );
    println!(
        "最大精度: {} 位",
        GlobalPrecisionConfig::get_max_precision()
    );
    println!(
        "最小精度: {} 位",
        GlobalPrecisionConfig::get_min_precision()
    );

    // 修改精度
    GlobalPrecisionConfig::set_default_precision(512)?;
    println!(
        "修改后默认精度: {} 位",
        GlobalPrecisionConfig::get_default_precision()
    );

    // 重置
    GlobalPrecisionConfig::reset_to_default();
    println!(
        "重置后默认精度: {} 位",
        GlobalPrecisionConfig::get_default_precision()
    );

    // 2. 测试Mpz乘法配置
    println!("\n📊 2. Mpz乘法配置测试");
    println!("{}", "-".repeat(30));

    println!(
        "当前乘法后端: {:?}",
        MpzMultiplicationConfig::get_global_backend()
    );

    // 设置不同的后端
    MpzMultiplicationConfig::set_global_backend(MultiplicationBackend::Karatsuba)?;
    println!(
        "设置Karatsuba后端后: {:?}",
        MpzMultiplicationConfig::get_global_backend()
    );

    MpzMultiplicationConfig::set_global_backend(MultiplicationBackend::FFT)?;
    println!(
        "设置FFT后端后: {:?}",
        MpzMultiplicationConfig::get_global_backend()
    );

    // 测试并行配置
    MpzMultiplicationConfig::enable_parallel(None)?;
    println!(
        "并行计算状态: {}",
        MpzMultiplicationConfig::is_parallel_enabled()
    );

    MpzMultiplicationConfig::set_parallel_threshold(2000)?;
    println!(
        "并行阈值: {} 位",
        MpzMultiplicationConfig::get_parallel_threshold()
    );

    // 测试建议后端
    let suggested = MpzMultiplicationConfig::suggest_backend(1000, 1000);
    println!("1000位数字建议后端: {:?}", suggested);

    let suggested = MpzMultiplicationConfig::suggest_backend(5000, 5000);
    println!("5000位数字建议后端: {:?}", suggested);

    // 重置
    MpzMultiplicationConfig::reset_to_default();
    println!(
        "重置后乘法后端: {:?}",
        MpzMultiplicationConfig::get_global_backend()
    );

    // 3. 测试Mpf配置
    println!("\n📊 3. Mpf配置测试");
    println!("{}", "-".repeat(30));

    println!("默认精度: {} 位", MpfConfig::get_default_precision());
    println!("最大精度: {} 位", MpfConfig::get_max_precision());
    println!("最小精度: {} 位", MpfConfig::get_min_precision());

    // 修改精度
    MpfConfig::set_default_precision(128)?;
    println!("修改后默认精度: {} 位", MpfConfig::get_default_precision());

    // 测试舍入模式
    println!("当前舍入模式: {:?}", MpfConfig::get_rounding_mode());

    // 测试动态精度
    println!(
        "动态精度调整: {}",
        MpfConfig::is_dynamic_precision_enabled()
    );

    // 获取完整配置
    let config = MpfConfig::get_config();
    println!("完整配置: {:?}", config);

    // 重置
    MpfConfig::reset_to_default();
    println!("重置后默认精度: {} 位", MpfConfig::get_default_precision());

    // 4. 测试Complex配置
    println!("\n📊 4. Complex配置测试");
    println!("{}", "-".repeat(30));

    println!("默认精度: {} 位", ComplexConfig::get_default_precision());
    println!("最大精度: {} 位", ComplexConfig::get_max_precision());
    println!("最小精度: {} 位", ComplexConfig::get_min_precision());

    // 测试高精度模式
    println!(
        "高精度模式: {}",
        ComplexConfig::is_high_precision_mode_enabled()
    );
    ComplexConfig::enable_high_precision_mode();
    println!(
        "启用高精度模式后: {}",
        ComplexConfig::is_high_precision_mode_enabled()
    );

    // 测试缓存配置
    println!("缓存优化: {}", ComplexConfig::is_caching_enabled());
    ComplexConfig::disable_caching();
    println!("禁用缓存后: {}", ComplexConfig::is_caching_enabled());

    // 测试并行配置
    println!("并行计算: {}", ComplexConfig::is_parallel_enabled());
    ComplexConfig::enable_parallel();
    println!("启用并行后: {}", ComplexConfig::is_parallel_enabled());

    // 获取完整配置
    let precision_config = ComplexConfig::get_precision_config();
    let cache_config = ComplexConfig::get_cache_config();
    let parallel_config = ComplexConfig::get_parallel_config();

    println!("精度配置: {:?}", precision_config);
    println!("缓存配置: {:?}", cache_config);
    println!("并行配置: {:?}", parallel_config);

    // 测试预设配置
    println!("\n📊 5. 预设配置测试");
    println!("{}", "-".repeat(30));

    // 高精度配置
    ComplexConfig::high_precision();
    println!(
        "高精度配置 - 精度: {} 位, 高精度模式: {}",
        ComplexConfig::get_default_precision(),
        ComplexConfig::is_high_precision_mode_enabled()
    );

    // 性能优化配置
    ComplexConfig::performance();
    println!(
        "性能配置 - 精度: {} 位, 并行: {}",
        ComplexConfig::get_default_precision(),
        ComplexConfig::is_parallel_enabled()
    );

    // 内存优化配置
    ComplexConfig::memory_efficient();
    println!(
        "内存配置 - 精度: {} 位, 缓存: {}",
        ComplexConfig::get_default_precision(),
        ComplexConfig::is_caching_enabled()
    );

    // 重置
    ComplexConfig::reset_to_default();
    println!(
        "重置后默认精度: {} 位",
        ComplexConfig::get_default_precision()
    );

    println!("\n✅ 配置系统测试完成！");
    Ok(())
}
