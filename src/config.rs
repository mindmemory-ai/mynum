//! 全局配置模块
//!
//! 管理全局精度设置、舍入模式等通用配置。

#[cfg(not(feature = "std"))]
use alloc::{boxed::Box, string::String, vec::Vec};

#[cfg(feature = "std")]
use std::sync::{Mutex, OnceLock};

use crate::error::{Error, Result};
use core::marker::PhantomData;
use core::sync::atomic::{AtomicUsize, Ordering};

/// 可原子存储和加载的配置值
pub trait AtomicConfigValue: Copy + 'static {
    /// Convert this value to a `usize` for atomic storage.
    fn to_usize(self) -> usize;
    /// Create an instance of this value from its `usize` representation.
    fn from_usize(v: usize) -> Self;
}

impl AtomicConfigValue for usize {
    fn to_usize(self) -> usize {
        self
    }
    fn from_usize(v: usize) -> Self {
        v
    }
}

/// 原子配置存储 — 泛型包装 AtomicUsize
pub struct AtomicConfig<T: AtomicConfigValue> {
    inner: AtomicUsize,
    _phantom: PhantomData<T>,
}

impl<T: AtomicConfigValue> AtomicConfig<T> {
    /// Create a new `AtomicConfig` with the given default value.
    pub fn new(default: T) -> Self {
        Self {
            inner: AtomicUsize::new(default.to_usize()),
            _phantom: PhantomData,
        }
    }

    /// Load the current value.
    pub fn load(&self) -> T {
        T::from_usize(self.inner.load(Ordering::Relaxed))
    }

    /// Store a new value, notifying config change listeners.
    pub fn store(&self, value: T) {
        self.inner.store(value.to_usize(), Ordering::SeqCst);
        notify_config_change();
    }
}

/// Registry for config change callbacks.
#[cfg(feature = "std")]
type ConfigCallback = Box<dyn Fn() + Send + Sync>;

#[cfg(feature = "std")]
static CONFIG_CALLBACKS: OnceLock<Mutex<Vec<ConfigCallback>>> = OnceLock::new();

#[cfg(feature = "std")]
fn callbacks() -> &'static Mutex<Vec<ConfigCallback>> {
    CONFIG_CALLBACKS.get_or_init(|| Mutex::new(Vec::new()))
}

/// Register a callback to be invoked when any AtomicConfig value changes.
/// Returns a CallbackHandle; callbacks are unregistered when the handle is dropped.
#[cfg(feature = "std")]
pub struct CallbackHandle {
    index: usize,
}

#[cfg(feature = "std")]
impl Drop for CallbackHandle {
    fn drop(&mut self) {
        if let Ok(mut cbs) = callbacks().lock() {
            if self.index < cbs.len() {
                cbs[self.index] = Box::new(|| {});
            }
        }
    }
}

/// Register a callback to be invoked when any AtomicConfig value changes.
/// The returned handle unregisters the callback when dropped.
#[cfg(feature = "std")]
pub fn on_config_change<F: Fn() + Send + Sync + 'static>(f: F) -> CallbackHandle {
    let mut cbs = callbacks().lock().unwrap();
    let index = cbs.len();
    cbs.push(Box::new(f));
    CallbackHandle { index }
}

/// Called by AtomicConfig::store() after updating the value.
pub(crate) fn notify_config_change() {
    #[cfg(feature = "std")]
    if let Some(cbs) = CONFIG_CALLBACKS.get() {
        if let Ok(cbs) = cbs.lock() {
            for cb in cbs.iter() {
                cb();
            }
        }
    }
}

/// 舍入模式
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RoundingMode {
    /// 向零舍入
    TowardZero,
    /// 向正无穷舍入
    TowardPositive,
    /// 向负无穷舍入
    TowardNegative,
    /// 向最近偶数舍入
    #[default]
    TowardNearest,
}

/// 全局精度配置
pub struct GlobalPrecisionConfig;

// 全局精度设置（使用原子操作确保线程安全）
static DEFAULT_PRECISION: AtomicUsize = AtomicUsize::new(256);
static MAX_PRECISION: AtomicUsize = AtomicUsize::new(1_000_000);
static MIN_PRECISION: AtomicUsize = AtomicUsize::new(32);

impl AtomicConfigValue for RoundingMode {
    fn to_usize(self) -> usize {
        self as usize
    }
    fn from_usize(v: usize) -> Self {
        match v {
            0 => RoundingMode::TowardZero,
            1 => RoundingMode::TowardPositive,
            2 => RoundingMode::TowardNegative,
            _ => RoundingMode::TowardNearest,
        }
    }
}

// 将 RoundingMode 存储为 usize（变体判别值）
// 使用直接 AtomicUsize 因为 static 初始化需要 const
static ROUNDING_MODE: AtomicUsize = AtomicUsize::new(3); // TowardNearest = 3

impl GlobalPrecisionConfig {
    /// 设置全局默认精度
    pub fn set_default_precision(precision: usize) -> Result<()> {
        let min_prec = MIN_PRECISION.load(Ordering::Relaxed);
        let max_prec = MAX_PRECISION.load(Ordering::Relaxed);

        if precision < min_prec || precision > max_prec {
            return Err(Error::InvalidInput(format!(
                "Precision must be between {} and {}",
                min_prec, max_prec
            )));
        }

        DEFAULT_PRECISION.store(precision, Ordering::Relaxed);
        notify_config_change();
        Ok(())
    }

    /// 获取当前默认精度
    pub fn get_default_precision() -> usize {
        DEFAULT_PRECISION.load(Ordering::Relaxed)
    }

    /// 设置全局最大精度
    pub fn set_max_precision(max_precision: usize) -> Result<()> {
        let min_prec = MIN_PRECISION.load(Ordering::Relaxed);
        let current_default = DEFAULT_PRECISION.load(Ordering::Relaxed);

        if max_precision < min_prec || max_precision < current_default {
            return Err(Error::InvalidInput(format!(
                "Max precision must be at least {} and >= current default {}",
                min_prec, current_default
            )));
        }

        MAX_PRECISION.store(max_precision, Ordering::Relaxed);
        notify_config_change();
        Ok(())
    }

    /// 获取当前最大精度
    pub fn get_max_precision() -> usize {
        MAX_PRECISION.load(Ordering::Relaxed)
    }

    /// 设置全局最小精度
    pub fn set_min_precision(min_precision: usize) -> Result<()> {
        let max_prec = MAX_PRECISION.load(Ordering::Relaxed);
        let current_default = DEFAULT_PRECISION.load(Ordering::Relaxed);

        if min_precision > max_prec || min_precision > current_default {
            return Err(Error::InvalidInput(format!(
                "Min precision must be <= {} and <= current default {}",
                max_prec, current_default
            )));
        }

        MIN_PRECISION.store(min_precision, Ordering::Relaxed);
        notify_config_change();
        Ok(())
    }

    /// 获取当前最小精度
    pub fn get_min_precision() -> usize {
        MIN_PRECISION.load(Ordering::Relaxed)
    }

    /// 设置全局舍入模式
    pub fn set_rounding_mode(mode: RoundingMode) {
        ROUNDING_MODE.store(mode as usize, Ordering::Relaxed);
        notify_config_change();
    }

    /// 获取当前舍入模式
    pub fn get_rounding_mode() -> RoundingMode {
        match ROUNDING_MODE.load(Ordering::Relaxed) {
            0 => RoundingMode::TowardZero,
            1 => RoundingMode::TowardPositive,
            2 => RoundingMode::TowardNegative,
            _ => RoundingMode::TowardNearest,
        }
    }

    /// 重置为默认配置
    pub fn reset_to_default() {
        DEFAULT_PRECISION.store(256, Ordering::Relaxed);
        MAX_PRECISION.store(1_000_000, Ordering::Relaxed);
        MIN_PRECISION.store(32, Ordering::Relaxed);
        ROUNDING_MODE.store(RoundingMode::default() as usize, Ordering::Relaxed);
        notify_config_change();
    }
}

impl Default for GlobalPrecisionConfig {
    fn default() -> Self {
        Self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rounding_mode_default() {
        GlobalPrecisionConfig::reset_to_default();
        assert_eq!(
            GlobalPrecisionConfig::get_rounding_mode(),
            RoundingMode::TowardNearest
        );
    }

    #[test]
    fn test_rounding_mode_set_and_get() {
        GlobalPrecisionConfig::set_rounding_mode(RoundingMode::TowardZero);
        assert_eq!(
            GlobalPrecisionConfig::get_rounding_mode(),
            RoundingMode::TowardZero
        );

        GlobalPrecisionConfig::set_rounding_mode(RoundingMode::TowardPositive);
        assert_eq!(
            GlobalPrecisionConfig::get_rounding_mode(),
            RoundingMode::TowardPositive
        );

        GlobalPrecisionConfig::set_rounding_mode(RoundingMode::TowardNegative);
        assert_eq!(
            GlobalPrecisionConfig::get_rounding_mode(),
            RoundingMode::TowardNegative
        );

        // Reset after test to avoid polluting other tests
        GlobalPrecisionConfig::reset_to_default();
    }

    #[test]
    fn test_reset_restores_rounding_mode() {
        GlobalPrecisionConfig::set_rounding_mode(RoundingMode::TowardZero);
        GlobalPrecisionConfig::reset_to_default();
        assert_eq!(
            GlobalPrecisionConfig::get_rounding_mode(),
            RoundingMode::TowardNearest
        );
    }

    #[test]
    #[cfg(feature = "std")]
    fn test_config_change_callback() {
        use std::sync::atomic::{AtomicBool, Ordering};
        use std::sync::Arc;

        let was_called = Arc::new(AtomicBool::new(false));
        let was_called_clone = was_called.clone();

        let _handle = on_config_change(move || {
            was_called_clone.store(true, Ordering::SeqCst);
        });

        GlobalPrecisionConfig::set_default_precision(512).unwrap();
        assert!(was_called.load(Ordering::SeqCst));

        // Reset after test
        GlobalPrecisionConfig::reset_to_default();
    }
}
