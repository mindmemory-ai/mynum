//! Mpz模块配置
//!
//! 管理大整数运算的算法后端、并行计算等配置。

use crate::error::{Error, Result};
use core::sync::atomic::{AtomicUsize, Ordering};

/// 大整数乘法算法后端
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MultiplicationBackend {
    /// 基础乘法算法（适用于小数）
    Schoolbook,
    /// Karatsuba算法（适用于中等大小的数）
    Karatsuba,
    /// Toom-Cook 3-way算法（适用于大数）
    ToomCook3,
    /// Toom-Cook 4-way算法（适用于极大数）
    ToomCook4,
    /// 基于FFT的乘法（适用于超大数）
    FFT,
    /// 基于NTT（数论变换）的乘法
    NTT,
    /// 自适应算法（根据数的大小自动选择最优算法）
    #[default]
    Adaptive,
    /// 并行Karatsuba乘法算法（适用于大数）
    ParallelKaratsuba,
    /// 并行FFT乘法算法（适用于超大数）
    ParallelFFT,
    /// 并行乘法算法（智能选择并行Karatsuba或并行FFT）
    Parallel,
    /// 自定义后端（用户提供的实现）
    Custom,
}

/// 不同算法的切换阈值（以limb数量为单位）
#[derive(Debug, Clone)]
pub struct AlgorithmThresholds {
    /// 从基础算法切换到Karatsuba的阈值
    pub schoolbook_to_karatsuba: usize,
    /// 从Karatsuba切换到Toom-Cook 3的阈值
    pub karatsuba_to_toom3: usize,
    /// 从Toom-Cook 3切换到Toom-Cook 4的阈值
    pub toom3_to_toom4: usize,
    /// 从Toom-Cook切换到FFT的阈值
    pub toom_to_fft: usize,
    /// FFT算法中使用NTT的阈值
    pub fft_to_ntt: usize,
}

impl Default for AlgorithmThresholds {
    fn default() -> Self {
        Self {
            schoolbook_to_karatsuba: 32,
            karatsuba_to_toom3: 256,
            toom3_to_toom4: 1024,
            toom_to_fft: 4096,
            fft_to_ntt: 16384,
        }
    }
}

/// Mpz乘法配置
pub struct MpzMultiplicationConfig;

// 全局乘法后端设置
static GLOBAL_BACKEND: AtomicUsize = AtomicUsize::new(MultiplicationBackend::Adaptive as usize);
static PARALLEL_ENABLED: AtomicUsize = AtomicUsize::new(0); // 0 = false, 1 = true
static PARALLEL_THRESHOLD: AtomicUsize = AtomicUsize::new(1000); // 并行计算的最小位数阈值

// 全局算法阈值设置
static SCHOOLBOOK_TO_KARATSUBA: AtomicUsize = AtomicUsize::new(32);
static KARATSUBA_TO_TOOM3: AtomicUsize = AtomicUsize::new(256);
static TOOM3_TO_TOOM4: AtomicUsize = AtomicUsize::new(1024);
static TOOM_TO_FFT: AtomicUsize = AtomicUsize::new(4096);
static FFT_TO_NTT: AtomicUsize = AtomicUsize::new(16384);

impl MpzMultiplicationConfig {
    /// 创建默认配置（自适应算法）
    pub fn new() -> Self {
        Self
    }

    /// 设置全局乘法后端
    pub fn set_global_backend(backend: MultiplicationBackend) -> Result<()> {
        GLOBAL_BACKEND.store(backend as usize, Ordering::Relaxed);
        Ok(())
    }

    /// 获取当前全局乘法后端
    pub fn get_global_backend() -> MultiplicationBackend {
        let backend_val = GLOBAL_BACKEND.load(Ordering::Relaxed);
        match backend_val {
            0 => MultiplicationBackend::Schoolbook,
            1 => MultiplicationBackend::Karatsuba,
            2 => MultiplicationBackend::ToomCook3,
            3 => MultiplicationBackend::ToomCook4,
            4 => MultiplicationBackend::FFT,
            5 => MultiplicationBackend::NTT,
            6 => MultiplicationBackend::Adaptive,
            7 => MultiplicationBackend::ParallelKaratsuba,
            8 => MultiplicationBackend::ParallelFFT,
            9 => MultiplicationBackend::Parallel,
            _ => MultiplicationBackend::Custom,
        }
    }

    /// 设置算法切换阈值
    pub fn set_thresholds(thresholds: AlgorithmThresholds) -> Result<()> {
        if thresholds.schoolbook_to_karatsuba == 0
            || thresholds.karatsuba_to_toom3 <= thresholds.schoolbook_to_karatsuba
            || thresholds.toom3_to_toom4 <= thresholds.karatsuba_to_toom3
            || thresholds.toom_to_fft <= thresholds.toom3_to_toom4
            || thresholds.fft_to_ntt <= thresholds.toom_to_fft
        {
            return Err(Error::InvalidInput(
                "Thresholds must be positive and strictly increasing".into(),
            ));
        }
        SCHOOLBOOK_TO_KARATSUBA.store(thresholds.schoolbook_to_karatsuba, Ordering::Relaxed);
        KARATSUBA_TO_TOOM3.store(thresholds.karatsuba_to_toom3, Ordering::Relaxed);
        TOOM3_TO_TOOM4.store(thresholds.toom3_to_toom4, Ordering::Relaxed);
        TOOM_TO_FFT.store(thresholds.toom_to_fft, Ordering::Relaxed);
        FFT_TO_NTT.store(thresholds.fft_to_ntt, Ordering::Relaxed);
        Ok(())
    }

    /// 获取当前算法切换阈值
    pub fn get_thresholds() -> AlgorithmThresholds {
        AlgorithmThresholds {
            schoolbook_to_karatsuba: SCHOOLBOOK_TO_KARATSUBA.load(Ordering::Relaxed),
            karatsuba_to_toom3: KARATSUBA_TO_TOOM3.load(Ordering::Relaxed),
            toom3_to_toom4: TOOM3_TO_TOOM4.load(Ordering::Relaxed),
            toom_to_fft: TOOM_TO_FFT.load(Ordering::Relaxed),
            fft_to_ntt: FFT_TO_NTT.load(Ordering::Relaxed),
        }
    }

    /// 启用并行计算
    pub fn enable_parallel(_thread_count: Option<usize>) -> Result<()> {
        PARALLEL_ENABLED.store(1, Ordering::Relaxed);
        Ok(())
    }

    /// 禁用并行计算
    pub fn disable_parallel() {
        PARALLEL_ENABLED.store(0, Ordering::Relaxed);
    }

    /// 获取并行计算状态
    pub fn is_parallel_enabled() -> bool {
        PARALLEL_ENABLED.load(Ordering::Relaxed) != 0
    }

    /// 设置并行计算阈值（位数）
    pub fn set_parallel_threshold(threshold: usize) -> Result<()> {
        if threshold < 100 {
            return Err(Error::InvalidInput(
                "Parallel threshold must be at least 100 bits".into(),
            ));
        }
        PARALLEL_THRESHOLD.store(threshold, Ordering::Relaxed);
        Ok(())
    }

    /// 获取并行计算阈值
    pub fn get_parallel_threshold() -> usize {
        PARALLEL_THRESHOLD.load(Ordering::Relaxed)
    }

    /// 重置为默认配置
    pub fn reset_to_default() {
        GLOBAL_BACKEND.store(MultiplicationBackend::Adaptive as usize, Ordering::Relaxed);
        PARALLEL_ENABLED.store(0, Ordering::Relaxed);
        PARALLEL_THRESHOLD.store(1000, Ordering::Relaxed);
        SCHOOLBOOK_TO_KARATSUBA.store(32, Ordering::Relaxed);
        KARATSUBA_TO_TOOM3.store(256, Ordering::Relaxed);
        TOOM3_TO_TOOM4.store(1024, Ordering::Relaxed);
        TOOM_TO_FFT.store(4096, Ordering::Relaxed);
        FFT_TO_NTT.store(16384, Ordering::Relaxed);
    }

    /// 根据操作数大小建议最优后端
    pub fn suggest_backend(operand1_bits: usize, operand2_bits: usize) -> MultiplicationBackend {
        let max_bits = operand1_bits.max(operand2_bits);
        let limbs = max_bits.div_ceil(64); // 假设64位limb

        let thresholds = Self::get_thresholds();

        // 检查是否应该使用并行算法
        if Self::is_parallel_enabled() && max_bits >= Self::get_parallel_threshold() {
            return MultiplicationBackend::Parallel;
        }

        if limbs < thresholds.schoolbook_to_karatsuba {
            MultiplicationBackend::Schoolbook
        } else if limbs < thresholds.karatsuba_to_toom3 {
            MultiplicationBackend::Karatsuba
        } else if limbs < thresholds.toom3_to_toom4 {
            MultiplicationBackend::ToomCook3
        } else if limbs < thresholds.toom_to_fft {
            MultiplicationBackend::ToomCook4
        } else if limbs < thresholds.fft_to_ntt {
            MultiplicationBackend::FFT
        } else {
            MultiplicationBackend::NTT
        }
    }

    /// 性能基准测试（用于调优阈值）
    ///
    /// 运行内部时序逻辑来确定各乘法算法之间的实际交叉点。
    /// 返回校准后的阈值。校准耗时数秒；若需进度反馈，
    /// 可运行 `cargo run --example threshold_calibration`。
    ///
    /// 在 `#[cfg(not(feature = "std"))]` 环境下，返回当前存储的阈值。
    ///
    /// 注意：此函数返回校准阈值但不自动存储。
    /// 调用 `set_thresholds(benchmark_thresholds())` 来应用。
    pub fn benchmark_thresholds() -> AlgorithmThresholds {
        calibrate_thresholds()
    }
}

#[cfg(feature = "std")]
fn calibrate_thresholds() -> AlgorithmThresholds {
    let s_to_k = find_crossover(
        MultiplicationBackend::Schoolbook,
        MultiplicationBackend::Karatsuba,
        16,
    );
    let k_to_t3 = find_crossover(
        MultiplicationBackend::Karatsuba,
        MultiplicationBackend::ToomCook3,
        64,
    );
    let t3_to_t4 = find_crossover(
        MultiplicationBackend::ToomCook3,
        MultiplicationBackend::ToomCook4,
        256,
    );
    let t_to_f = find_crossover(
        MultiplicationBackend::ToomCook4,
        MultiplicationBackend::FFT,
        512,
    );

    let toom_to_fft_val = t_to_f.unwrap_or(4096);
    // NTT 目前委托给 mul_fft (见 src/mpz/multiplication.rs:20)，
    // 因此不进行计时。合成值为 FFT 交叉点的 4 倍，下限为 16384 limbs。
    let fft_to_ntt_val = std::cmp::max(toom_to_fft_val * 4, 16384);

    AlgorithmThresholds {
        schoolbook_to_karatsuba: s_to_k.unwrap_or(32),
        karatsuba_to_toom3: k_to_t3.unwrap_or(256),
        toom3_to_toom4: t3_to_t4.unwrap_or(1024),
        toom_to_fft: toom_to_fft_val,
        fft_to_ntt: fft_to_ntt_val,
    }
}

#[cfg(feature = "std")]
fn find_crossover(
    backend_a: MultiplicationBackend,
    backend_b: MultiplicationBackend,
    start_bits: usize,
) -> Option<usize> {
    use std::time::Instant;

    let mut bits = start_bits;
    while bits <= 32768 {
        let a = crate::Mpz::random_bits(bits).ok()?;
        let b = crate::Mpz::random_bits(bits).ok()?;

        // Warmup: 2 iterations
        for _ in 0..2 {
            let _ = a.mul_with_backend(&b, backend_a);
            let _ = a.mul_with_backend(&b, backend_b);
        }

        // Measured: 5 iterations, collect times
        let mut times_a = Vec::with_capacity(5);
        let mut times_b = Vec::with_capacity(5);

        for _ in 0..5 {
            let start = Instant::now();
            let _ = a.mul_with_backend(&b, backend_a);
            times_a.push(start.elapsed());

            let start = Instant::now();
            let _ = a.mul_with_backend(&b, backend_b);
            times_b.push(start.elapsed());
        }

        times_a.sort();
        times_b.sort();

        if times_b[2] < times_a[2] {
            return Some(bits.div_ceil(64).max(1));
        }

        bits *= 2;
    }
    None
}

#[cfg(not(feature = "std"))]
fn calibrate_thresholds() -> AlgorithmThresholds {
    MpzMultiplicationConfig::get_thresholds()
}

impl Default for MpzMultiplicationConfig {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thresholds_default() {
        MpzMultiplicationConfig::reset_to_default();
        let t = MpzMultiplicationConfig::get_thresholds();
        assert_eq!(t.schoolbook_to_karatsuba, 32);
        assert_eq!(t.karatsuba_to_toom3, 256);
    }

    #[test]
    fn test_thresholds_set_and_get() {
        MpzMultiplicationConfig::reset_to_default();
        let custom = AlgorithmThresholds {
            schoolbook_to_karatsuba: 64,
            karatsuba_to_toom3: 512,
            toom3_to_toom4: 2048,
            toom_to_fft: 8192,
            fft_to_ntt: 32768,
        };
        MpzMultiplicationConfig::set_thresholds(custom).unwrap();
        let t = MpzMultiplicationConfig::get_thresholds();
        assert_eq!(t.schoolbook_to_karatsuba, 64);
        assert_eq!(t.fft_to_ntt, 32768);
        MpzMultiplicationConfig::reset_to_default();
    }

    #[test]
    fn test_thresholds_validation_rejects_non_increasing() {
        let bad = AlgorithmThresholds {
            schoolbook_to_karatsuba: 100,
            karatsuba_to_toom3: 50, // smaller than previous
            toom3_to_toom4: 200,
            toom_to_fft: 300,
            fft_to_ntt: 400,
        };
        assert!(MpzMultiplicationConfig::set_thresholds(bad).is_err());
    }

    #[test]
    fn test_thresholds_validation_rejects_zero() {
        let bad = AlgorithmThresholds {
            schoolbook_to_karatsuba: 0,
            karatsuba_to_toom3: 256,
            toom3_to_toom4: 1024,
            toom_to_fft: 4096,
            fft_to_ntt: 16384,
        };
        assert!(MpzMultiplicationConfig::set_thresholds(bad).is_err());
    }
}
