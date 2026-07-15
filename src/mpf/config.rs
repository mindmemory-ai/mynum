//! Mpf模块配置
//!
//! 管理高精度浮点数运算的精度、舍入模式等配置。

use crate::error::{Error, Result};
use core::sync::atomic::{AtomicUsize, Ordering};

/// 浮点数舍入模式
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MpfRoundingMode {
    /// 向零舍入
    TowardZero,
    /// 向正无穷舍入
    TowardPositive,
    /// 向负无穷舍入
    TowardNegative,
    /// 向最近偶数舍入
    TowardNearest,
    /// 向最近舍入（四舍五入）
    #[default]
    TowardNearestTiesToEven,
    /// 向最近舍入（四舍五入，平局时向正无穷）
    TowardNearestTiesToPositive,
}

/// 浮点数精度配置
#[derive(Debug, Clone)]
pub struct MpfPrecisionConfig {
    /// 默认精度（位数）
    pub default_precision: usize,
    /// 最大精度（位数）
    pub max_precision: usize,
    /// 最小精度（位数）
    pub min_precision: usize,
    /// 是否启用动态精度调整
    pub enable_dynamic_precision: bool,
    /// 精度调整步长
    pub precision_step: usize,
}

impl Default for MpfPrecisionConfig {
    fn default() -> Self {
        Self {
            default_precision: 256,
            max_precision: 1_000_000,
            min_precision: 64,
            enable_dynamic_precision: true,
            precision_step: 64,
        }
    }
}

/// 浮点数运算配置
pub struct MpfConfig;

// 全局浮点数配置设置
static MPF_DEFAULT_PRECISION: AtomicUsize = AtomicUsize::new(256);
static MPF_MAX_PRECISION: AtomicUsize = AtomicUsize::new(1_000_000);
static MPF_MIN_PRECISION: AtomicUsize = AtomicUsize::new(64);
static MPF_ROUNDING_MODE: AtomicUsize =
    AtomicUsize::new(MpfRoundingMode::TowardNearestTiesToEven as usize);
static MPF_DYNAMIC_PRECISION: AtomicUsize = AtomicUsize::new(1); // 1 = true, 0 = false

impl MpfConfig {
    /// 设置默认精度
    pub fn set_default_precision(precision: usize) -> Result<()> {
        let min_prec = MPF_MIN_PRECISION.load(Ordering::Relaxed);
        let max_prec = MPF_MAX_PRECISION.load(Ordering::Relaxed);

        if precision < min_prec || precision > max_prec {
            return Err(Error::InvalidInput(format!(
                "Precision must be between {} and {}",
                min_prec, max_prec
            )));
        }

        MPF_DEFAULT_PRECISION.store(precision, Ordering::Relaxed);
        Ok(())
    }

    /// 获取当前默认精度
    pub fn get_default_precision() -> usize {
        MPF_DEFAULT_PRECISION.load(Ordering::Relaxed)
    }

    /// 设置最大精度
    pub fn set_max_precision(max_precision: usize) -> Result<()> {
        let min_prec = MPF_MIN_PRECISION.load(Ordering::Relaxed);
        let current_default = MPF_DEFAULT_PRECISION.load(Ordering::Relaxed);

        if max_precision < min_prec || max_precision < current_default {
            return Err(Error::InvalidInput(format!(
                "Max precision must be at least {} and >= current default {}",
                min_prec, current_default
            )));
        }

        MPF_MAX_PRECISION.store(max_precision, Ordering::Relaxed);
        Ok(())
    }

    /// 获取当前最大精度
    pub fn get_max_precision() -> usize {
        MPF_MAX_PRECISION.load(Ordering::Relaxed)
    }

    /// 设置最小精度
    pub fn set_min_precision(min_precision: usize) -> Result<()> {
        let max_prec = MPF_MAX_PRECISION.load(Ordering::Relaxed);
        let current_default = MPF_DEFAULT_PRECISION.load(Ordering::Relaxed);

        if min_precision > max_prec || min_precision > current_default {
            return Err(Error::InvalidInput(format!(
                "Min precision must be <= {} and <= current default {}",
                max_prec, current_default
            )));
        }

        MPF_MIN_PRECISION.store(min_precision, Ordering::Relaxed);
        Ok(())
    }

    /// 获取当前最小精度
    pub fn get_min_precision() -> usize {
        MPF_MIN_PRECISION.load(Ordering::Relaxed)
    }

    /// 设置舍入模式
    pub fn set_rounding_mode(mode: MpfRoundingMode) {
        MPF_ROUNDING_MODE.store(mode as usize, Ordering::Relaxed);
    }

    /// 获取当前舍入模式
    pub fn get_rounding_mode() -> MpfRoundingMode {
        let mode_val = MPF_ROUNDING_MODE.load(Ordering::Relaxed);
        match mode_val {
            0 => MpfRoundingMode::TowardZero,
            1 => MpfRoundingMode::TowardPositive,
            2 => MpfRoundingMode::TowardNegative,
            3 => MpfRoundingMode::TowardNearest,
            4 => MpfRoundingMode::TowardNearestTiesToEven,
            5 => MpfRoundingMode::TowardNearestTiesToPositive,
            _ => MpfRoundingMode::TowardNearestTiesToEven,
        }
    }

    /// 启用动态精度调整
    pub fn enable_dynamic_precision() {
        MPF_DYNAMIC_PRECISION.store(1, Ordering::Relaxed);
    }

    /// 禁用动态精度调整
    pub fn disable_dynamic_precision() {
        MPF_DYNAMIC_PRECISION.store(0, Ordering::Relaxed);
    }

    /// 检查是否启用动态精度调整
    pub fn is_dynamic_precision_enabled() -> bool {
        MPF_DYNAMIC_PRECISION.load(Ordering::Relaxed) != 0
    }

    /// 重置为默认配置
    pub fn reset_to_default() {
        MPF_DEFAULT_PRECISION.store(256, Ordering::Relaxed);
        MPF_MAX_PRECISION.store(1_000_000, Ordering::Relaxed);
        MPF_MIN_PRECISION.store(64, Ordering::Relaxed);
        MPF_ROUNDING_MODE.store(
            MpfRoundingMode::TowardNearestTiesToEven as usize,
            Ordering::Relaxed,
        );
        MPF_DYNAMIC_PRECISION.store(1, Ordering::Relaxed);
    }

    /// 获取完整配置
    pub fn get_config() -> MpfPrecisionConfig {
        MpfPrecisionConfig {
            default_precision: Self::get_default_precision(),
            max_precision: Self::get_max_precision(),
            min_precision: Self::get_min_precision(),
            enable_dynamic_precision: Self::is_dynamic_precision_enabled(),
            precision_step: 64, // 固定步长
        }
    }
}

impl Default for MpfConfig {
    fn default() -> Self {
        Self
    }
}
