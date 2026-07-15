//! Mpf 自适应精度模块
//!
//! 支持在迭代计算中动态调整精度以优化性能。

/// 自适应精度控制器
#[derive(Debug, Clone)]
pub struct AdaptivePrecision {
    /// 起始精度（位）
    pub start_precision: usize,
    /// 目标精度（位）
    pub target_precision: usize,
    /// 精度增长因子（每次迭代乘以此值）
    pub growth_factor: f64,
    /// 当前精度
    current_precision: usize,
    /// 当前迭代
    iteration: usize,
}

impl AdaptivePrecision {
    /// 创建新的自适应精度控制器
    ///
    /// 从 `start_precision` 开始，按 `growth_factor` 增长，
    /// 直到达到 `target_precision`。
    pub fn new(start_precision: usize, target_precision: usize, growth_factor: f64) -> Self {
        assert!(start_precision >= 32);
        assert!(target_precision >= start_precision);
        Self {
            start_precision,
            target_precision,
            growth_factor,
            current_precision: start_precision,
            iteration: 0,
        }
    }

    /// 获取当前迭代的推荐精度
    pub fn current_precision(&self) -> usize {
        self.current_precision.min(self.target_precision)
    }

    /// 当前迭代次数
    pub fn iteration(&self) -> usize {
        self.iteration
    }

    /// 进入下一次迭代（可能增加精度）
    pub fn advance(&mut self) {
        self.iteration += 1;
        let next =
            (self.start_precision as f64 * self.growth_factor.powi(self.iteration as i32)) as usize;
        self.current_precision = next.min(self.target_precision);
    }

    /// 是否已达到目标精度
    pub fn is_done(&self) -> bool {
        self.current_precision >= self.target_precision
    }

    /// 重置到起始状态
    pub fn reset(&mut self) {
        self.current_precision = self.start_precision;
        self.iteration = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_adaptive_precision_growth() {
        let mut ap = AdaptivePrecision::new(64, 1024, 2.0);
        assert_eq!(ap.current_precision(), 64);
        ap.advance();
        assert_eq!(ap.current_precision(), 128);
        ap.advance();
        assert_eq!(ap.current_precision(), 256);
        ap.advance();
        assert_eq!(ap.current_precision(), 512);
        ap.advance();
        assert_eq!(ap.current_precision(), 1024);
        assert!(ap.is_done());
    }

    #[test]
    fn test_adaptive_precision_capped() {
        let mut ap = AdaptivePrecision::new(64, 100, 2.0);
        ap.advance(); // 128 -> capped to 100
        assert_eq!(ap.current_precision(), 100);
        assert!(ap.is_done());
    }

    #[test]
    fn test_adaptive_precision_reset() {
        let mut ap = AdaptivePrecision::new(64, 256, 2.0);
        ap.advance(); // 128
        ap.advance(); // 256
        assert!(ap.is_done());
        ap.reset();
        assert_eq!(ap.current_precision(), 64);
        assert!(!ap.is_done());
    }
}
