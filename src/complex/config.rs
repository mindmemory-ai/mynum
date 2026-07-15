//! Complex模块配置
//!
//! 管理复数运算的精度、缓存、并行计算等配置。

use crate::error::{Error, Result};
use core::sync::atomic::{AtomicUsize, Ordering};

/// 复数运算精度配置
#[derive(Debug, Clone)]
pub struct ComplexPrecisionConfig {
    /// 默认精度（位数）
    pub default_precision: usize,
    /// 最大精度（位数）
    pub max_precision: usize,
    /// 最小精度（位数）
    pub min_precision: usize,
    /// 是否启用高精度模式
    pub high_precision_mode: bool,
    /// 精度调整步长
    pub precision_step: usize,
}

impl Default for ComplexPrecisionConfig {
    fn default() -> Self {
        Self {
            default_precision: 64,
            max_precision: 1_000_000,
            min_precision: 32,
            high_precision_mode: false,
            precision_step: 32,
        }
    }
}

/// 复数运算缓存配置
#[derive(Debug, Clone)]
pub struct ComplexCacheConfig {
    /// 是否启用缓存优化
    pub enable_caching: bool,
    /// 缓存大小限制
    pub cache_size_limit: usize,
    /// 缓存过期时间（毫秒）
    pub cache_expiry_ms: u64,
    /// 是否启用LRU缓存策略
    pub enable_lru: bool,
}

impl Default for ComplexCacheConfig {
    fn default() -> Self {
        Self {
            enable_caching: true,
            cache_size_limit: 1000,
            cache_expiry_ms: 300000, // 5分钟
            enable_lru: true,
        }
    }
}

/// 复数运算并行配置
#[derive(Debug, Clone)]
pub struct ComplexParallelConfig {
    /// 是否启用并行计算
    pub enable_parallel: bool,
    /// 并行计算阈值（位数）
    pub parallel_threshold: usize,
    /// 最大并行线程数
    pub max_threads: usize,
    /// 是否启用负载均衡
    pub enable_load_balancing: bool,
}

impl Default for ComplexParallelConfig {
    fn default() -> Self {
        Self {
            enable_parallel: false,
            parallel_threshold: 1000,
            max_threads: 4,
            enable_load_balancing: true,
        }
    }
}

/// 复数运算配置
#[derive(Debug, Clone)]
pub struct ComplexConfig;

// 全局复数配置设置
static COMPLEX_DEFAULT_PRECISION: AtomicUsize = AtomicUsize::new(64);
static COMPLEX_MAX_PRECISION: AtomicUsize = AtomicUsize::new(1_000_000);
static COMPLEX_MIN_PRECISION: AtomicUsize = AtomicUsize::new(32);
static COMPLEX_HIGH_PRECISION_MODE: AtomicUsize = AtomicUsize::new(0); // 0 = false, 1 = true
static COMPLEX_ENABLE_CACHING: AtomicUsize = AtomicUsize::new(1); // 1 = true, 0 = false
static COMPLEX_CACHE_SIZE_LIMIT: AtomicUsize = AtomicUsize::new(1000);
static COMPLEX_ENABLE_PARALLEL: AtomicUsize = AtomicUsize::new(0); // 0 = false, 1 = true
static COMPLEX_PARALLEL_THRESHOLD: AtomicUsize = AtomicUsize::new(1000);

impl ComplexConfig {
    /// 设置默认精度
    pub fn set_default_precision(precision: usize) -> Result<()> {
        let min_prec = COMPLEX_MIN_PRECISION.load(Ordering::Relaxed);
        let max_prec = COMPLEX_MAX_PRECISION.load(Ordering::Relaxed);

        if precision < min_prec || precision > max_prec {
            return Err(Error::InvalidInput(format!(
                "Precision must be between {} and {}",
                min_prec, max_prec
            )));
        }

        COMPLEX_DEFAULT_PRECISION.store(precision, Ordering::Relaxed);
        Ok(())
    }

    /// 获取当前默认精度
    pub fn get_default_precision() -> usize {
        COMPLEX_DEFAULT_PRECISION.load(Ordering::Relaxed)
    }

    /// 设置最大精度
    pub fn set_max_precision(max_precision: usize) -> Result<()> {
        let min_prec = COMPLEX_MIN_PRECISION.load(Ordering::Relaxed);
        let current_default = COMPLEX_DEFAULT_PRECISION.load(Ordering::Relaxed);

        if max_precision < min_prec || max_precision < current_default {
            return Err(Error::InvalidInput(format!(
                "Max precision must be at least {} and >= current default {}",
                min_prec, current_default
            )));
        }

        COMPLEX_MAX_PRECISION.store(max_precision, Ordering::Relaxed);
        Ok(())
    }

    /// 获取当前最大精度
    pub fn get_max_precision() -> usize {
        COMPLEX_MAX_PRECISION.load(Ordering::Relaxed)
    }

    /// 设置最小精度
    pub fn set_min_precision(min_precision: usize) -> Result<()> {
        let max_prec = COMPLEX_MAX_PRECISION.load(Ordering::Relaxed);
        let current_default = COMPLEX_DEFAULT_PRECISION.load(Ordering::Relaxed);

        if min_precision > max_prec || min_precision > current_default {
            return Err(Error::InvalidInput(format!(
                "Min precision must be <= {} and <= current default {}",
                max_prec, current_default
            )));
        }

        COMPLEX_MIN_PRECISION.store(min_precision, Ordering::Relaxed);
        Ok(())
    }

    /// 获取当前最小精度
    pub fn get_min_precision() -> usize {
        COMPLEX_MIN_PRECISION.load(Ordering::Relaxed)
    }

    /// 启用高精度模式
    pub fn enable_high_precision_mode() {
        COMPLEX_HIGH_PRECISION_MODE.store(1, Ordering::Relaxed);
    }

    /// 禁用高精度模式
    pub fn disable_high_precision_mode() {
        COMPLEX_HIGH_PRECISION_MODE.store(0, Ordering::Relaxed);
    }

    /// 检查是否启用高精度模式
    pub fn is_high_precision_mode_enabled() -> bool {
        COMPLEX_HIGH_PRECISION_MODE.load(Ordering::Relaxed) != 0
    }

    /// 启用缓存优化
    pub fn enable_caching() {
        COMPLEX_ENABLE_CACHING.store(1, Ordering::Relaxed);
    }

    /// 禁用缓存优化
    pub fn disable_caching() {
        COMPLEX_ENABLE_CACHING.store(0, Ordering::Relaxed);
    }

    /// 检查是否启用缓存优化
    pub fn is_caching_enabled() -> bool {
        COMPLEX_ENABLE_CACHING.load(Ordering::Relaxed) != 0
    }

    /// 设置缓存大小限制
    pub fn set_cache_size_limit(limit: usize) -> Result<()> {
        if limit < 100 {
            return Err(Error::InvalidInput(
                "Cache size limit must be at least 100".into(),
            ));
        }
        COMPLEX_CACHE_SIZE_LIMIT.store(limit, Ordering::Relaxed);
        Ok(())
    }

    /// 获取当前缓存大小限制
    pub fn get_cache_size_limit() -> usize {
        COMPLEX_CACHE_SIZE_LIMIT.load(Ordering::Relaxed)
    }

    /// 启用并行计算
    pub fn enable_parallel() {
        COMPLEX_ENABLE_PARALLEL.store(1, Ordering::Relaxed);
    }

    /// 禁用并行计算
    pub fn disable_parallel() {
        COMPLEX_ENABLE_PARALLEL.store(0, Ordering::Relaxed);
    }

    /// 检查是否启用并行计算
    pub fn is_parallel_enabled() -> bool {
        COMPLEX_ENABLE_PARALLEL.load(Ordering::Relaxed) != 0
    }

    /// 设置并行计算阈值
    pub fn set_parallel_threshold(threshold: usize) -> Result<()> {
        if threshold < 100 {
            return Err(Error::InvalidInput(
                "Parallel threshold must be at least 100 bits".into(),
            ));
        }
        COMPLEX_PARALLEL_THRESHOLD.store(threshold, Ordering::Relaxed);
        Ok(())
    }

    /// 获取当前并行计算阈值
    pub fn get_parallel_threshold() -> usize {
        COMPLEX_PARALLEL_THRESHOLD.load(Ordering::Relaxed)
    }

    /// 重置为默认配置
    pub fn reset_to_default() {
        COMPLEX_DEFAULT_PRECISION.store(64, Ordering::Relaxed);
        COMPLEX_MAX_PRECISION.store(1_000_000, Ordering::Relaxed);
        COMPLEX_MIN_PRECISION.store(32, Ordering::Relaxed);
        COMPLEX_HIGH_PRECISION_MODE.store(0, Ordering::Relaxed);
        COMPLEX_ENABLE_CACHING.store(1, Ordering::Relaxed);
        COMPLEX_CACHE_SIZE_LIMIT.store(1000, Ordering::Relaxed);
        COMPLEX_ENABLE_PARALLEL.store(0, Ordering::Relaxed);
        COMPLEX_PARALLEL_THRESHOLD.store(1000, Ordering::Relaxed);
    }

    /// 获取完整精度配置
    pub fn get_precision_config() -> ComplexPrecisionConfig {
        ComplexPrecisionConfig {
            default_precision: Self::get_default_precision(),
            max_precision: Self::get_max_precision(),
            min_precision: Self::get_min_precision(),
            high_precision_mode: Self::is_high_precision_mode_enabled(),
            precision_step: 32, // 固定步长
        }
    }

    /// 获取完整缓存配置
    pub fn get_cache_config() -> ComplexCacheConfig {
        ComplexCacheConfig {
            enable_caching: Self::is_caching_enabled(),
            cache_size_limit: Self::get_cache_size_limit(),
            cache_expiry_ms: 300000, // 固定值
            enable_lru: true,        // 固定值
        }
    }

    /// 获取完整并行配置
    pub fn get_parallel_config() -> ComplexParallelConfig {
        ComplexParallelConfig {
            enable_parallel: Self::is_parallel_enabled(),
            parallel_threshold: Self::get_parallel_threshold(),
            max_threads: 4,              // 固定值
            enable_load_balancing: true, // 固定值
        }
    }

    /// 创建高精度配置
    pub fn high_precision() -> Self {
        Self::enable_high_precision_mode();
        Self::set_default_precision(128).unwrap_or_default();
        Self::enable_caching();
        Self::set_cache_size_limit(2000).unwrap_or_default();
        Self::enable_parallel();
        Self
    }

    /// 创建性能优化配置
    pub fn performance() -> Self {
        Self::disable_high_precision_mode();
        Self::set_default_precision(64).unwrap_or_default();
        Self::enable_caching();
        Self::set_cache_size_limit(5000).unwrap_or_default();
        Self::enable_parallel();
        Self
    }

    /// 创建内存优化配置
    pub fn memory_efficient() -> Self {
        Self::disable_high_precision_mode();
        Self::set_default_precision(32).unwrap_or_default();
        Self::disable_caching();
        Self::set_cache_size_limit(100).unwrap_or_default();
        Self::disable_parallel();
        Self
    }
}

impl Default for ComplexConfig {
    fn default() -> Self {
        Self
    }
}
