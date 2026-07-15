//! 性能分析钩子
//!
//! 提供可选的性能分析基础设施，兼容 flamegraph。

use std::time::Instant;

/// 计时器守卫 — 当 dropped 时记录经过的时间
pub struct Timer {
    #[allow(dead_code)]
    name: &'static str,
    start: Instant,
}

impl Timer {
    /// Start a named timer for profiling. The elapsed time is logged
    /// on drop (flamegraph-compatible).
    pub fn new(name: &'static str) -> Self {
        Self {
            name,
            start: Instant::now(),
        }
    }
}

impl Drop for Timer {
    fn drop(&mut self) {
        let elapsed = self.start.elapsed();
        let _ = (self.name, elapsed);
    }
}

/// 创建一个性能分析计时器
///
/// 用法:
/// ```ignore
/// let _t = mynum::profile::timer("mpz::mul_fft");
/// // ... 被分析的代码 ...
/// // _t 在 drop 时自动记录计时
/// ```
#[inline]
pub fn timer(name: &'static str) -> Timer {
    Timer::new(name)
}
