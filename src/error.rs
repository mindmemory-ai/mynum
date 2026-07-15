//! 错误处理模块
//!
//! 定义了库中使用的所有错误类型和结果类型。

#[cfg(not(feature = "std"))]
use alloc::{format, string::String};

/// 库中使用的主要错误类型
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// 除零错误
    DivisionByZero,

    /// 无效的输入
    InvalidInput(String),

    /// 精度溢出
    PrecisionOverflow,

    /// 数值溢出
    Overflow,

    /// 数值下溢
    Underflow,

    /// 域错误（如负数开方）
    DomainError(String),

    /// 解析错误
    ParseError(String),

    /// 舍入错误
    RoundingError,

    /// 收敛错误（迭代算法未能在最大迭代次数内收敛）
    ConvergenceError(String),

    /// 功能未实现
    NotImplemented(String),

    /// 内存分配错误
    MemoryError,

    /// 其他错误
    Other(String),
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Error::DivisionByZero => write!(f, "Division by zero"),
            Error::InvalidInput(msg) => write!(f, "Invalid input: {}", msg),
            Error::PrecisionOverflow => write!(f, "Precision overflow"),
            Error::Overflow => write!(f, "Numerical overflow"),
            Error::Underflow => write!(f, "Numerical underflow"),
            Error::DomainError(msg) => write!(f, "Domain error: {}", msg),
            Error::ParseError(msg) => write!(f, "Parse error: {}", msg),
            Error::RoundingError => write!(f, "Rounding error"),
            Error::ConvergenceError(msg) => write!(f, "Convergence error: {}", msg),
            Error::NotImplemented(msg) => write!(f, "Not implemented: {}", msg),
            Error::MemoryError => write!(f, "Memory allocation error"),
            Error::Other(msg) => write!(f, "Error: {}", msg),
        }
    }
}

/// 错误严重级别，用于指导恢复决策。
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ErrorSeverity {
    /// 警告：结果可能不精确但仍可用（如精度损失）。
    Warning,
    /// 错误：操作失败但系统状态一致（可重试）。
    Error,
    /// 致命：系统状态可能不一致（不可恢复）。
    Fatal,
}

impl Error {
    /// 评估错误的严重级别以指导恢复决策。
    pub fn severity(&self) -> ErrorSeverity {
        match self {
            Error::Underflow => ErrorSeverity::Warning,
            Error::RoundingError => ErrorSeverity::Warning,
            Error::ConvergenceError(_) => ErrorSeverity::Error,
            Error::PrecisionOverflow => ErrorSeverity::Error,
            Error::Overflow => ErrorSeverity::Error,
            Error::DomainError(_) => ErrorSeverity::Error,
            Error::InvalidInput(_) => ErrorSeverity::Error,
            Error::ParseError(_) => ErrorSeverity::Error,
            Error::NotImplemented(_) => ErrorSeverity::Error,
            Error::DivisionByZero => ErrorSeverity::Fatal,
            Error::MemoryError => ErrorSeverity::Fatal,
            Error::Other(_) => ErrorSeverity::Error,
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None // 当前没有嵌套错误，为将来预留
    }
}

impl Error {
    /// 创建域错误
    pub fn domain(msg: impl Into<String>) -> Self {
        Error::DomainError(msg.into())
    }

    /// 创建无效输入错误
    pub fn invalid_input(msg: impl Into<String>) -> Self {
        Error::InvalidInput(msg.into())
    }

    /// 创建收敛错误
    pub fn convergence(msg: impl Into<String>) -> Self {
        Error::ConvergenceError(msg.into())
    }

    /// 创建解析错误
    pub fn parse(msg: impl Into<String>) -> Self {
        Error::ParseError(msg.into())
    }

    /// 为错误附加上下文
    pub fn with_context(self, context: impl Into<String>) -> Self {
        let ctx = context.into();
        match self {
            Error::DomainError(msg) => Error::DomainError(format!("{}: {}", ctx, msg)),
            Error::InvalidInput(msg) => Error::InvalidInput(format!("{}: {}", ctx, msg)),
            Error::ParseError(msg) => Error::ParseError(format!("{}: {}", ctx, msg)),
            Error::ConvergenceError(msg) => Error::ConvergenceError(format!("{}: {}", ctx, msg)),
            Error::Other(msg) => Error::Other(format!("{}: {}", ctx, msg)),
            other => Error::Other(format!("{}: {}", ctx, other)),
        }
    }
}

/// 可配置的错误恢复策略。
///
/// 定义在遇到可重试错误时如何重试操作。
/// `Warning` 和 `Error` 级别触发重试，`Fatal` 立即终止不重试。
#[derive(Debug, Clone)]
pub struct RecoveryStrategy {
    /// 最大尝试次数（含首次尝试）。
    pub max_retries: usize,
    /// 每次重试的精度增长因子（用于自适应精度循环）。
    pub precision_growth: f64,
}

impl Default for RecoveryStrategy {
    fn default() -> Self {
        Self {
            max_retries: 3,
            precision_growth: 2.0,
        }
    }
}

impl RecoveryStrategy {
    /// 创建新的恢复策略。
    ///
    /// `max_retries` 必须 >= 1。
    pub fn new(max_retries: usize, precision_growth: f64) -> Self {
        assert!(max_retries >= 1, "max_retries must be >= 1");
        Self {
            max_retries,
            precision_growth,
        }
    }

    /// 禁用恢复（仅尝试一次，失败立即返回）。
    pub fn no_recovery() -> Self {
        Self {
            max_retries: 1,
            precision_growth: 1.0,
        }
    }

    /// 从此策略创建 `AdaptivePrecision` 控制器。
    ///
    /// 使用 `self.precision_growth` 作为精度增长因子。
    pub fn to_adaptive(
        &self,
        start: usize,
        target: usize,
    ) -> crate::mpf::adaptive::AdaptivePrecision {
        crate::mpf::adaptive::AdaptivePrecision::new(start, target, self.precision_growth)
    }

    /// 重试一个可能失败的操作至多 `max_retries` 次。
    ///
    /// 若操作返回 `Ok`，立即返回。遇到 `Fatal` 错误立即返回；
    /// 其他错误重试直至尝试次数耗尽，然后返回最后一次错误。
    pub fn retry<F, T>(&self, mut f: F) -> Result<T>
    where
        F: FnMut() -> Result<T>,
    {
        let mut last_err = None;
        for _ in 0..self.max_retries {
            match f() {
                Ok(val) => return Ok(val),
                Err(err) => {
                    if err.severity() == ErrorSeverity::Fatal {
                        return Err(err);
                    }
                    last_err = Some(err);
                }
            }
        }
        Err(last_err.unwrap_or_else(|| Error::Other("retry exhausted".into())))
    }

    /// 每次重试时递增精度后重试。
    ///
    /// 每次失败后调用 `adaptive.advance()` 增加精度再重试。
    /// 适用于迭代算法（积分、级数求和、求根）。
    pub fn retry_with_precision<F, T>(
        &self,
        mut f: F,
        adaptive: &mut crate::mpf::adaptive::AdaptivePrecision,
    ) -> Result<T>
    where
        F: FnMut(usize) -> Result<T>,
    {
        let mut last_err = None;
        for attempt in 0..self.max_retries {
            let precision = adaptive.current_precision();
            match f(precision) {
                Ok(val) => return Ok(val),
                Err(err) => {
                    if err.severity() == ErrorSeverity::Fatal {
                        return Err(err);
                    }
                    last_err = Some(err);
                    if attempt + 1 < self.max_retries {
                        adaptive.advance();
                    }
                }
            }
        }
        Err(last_err.unwrap_or_else(|| Error::Other("precision retry exhausted".into())))
    }
}

/// 库中使用的结果类型
pub type Result<T> = core::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convenience_constructors() {
        let err = Error::domain("negative sqrt");
        assert_eq!(format!("{}", err), "Domain error: negative sqrt");

        let err = Error::convergence("Newton method");
        assert_eq!(format!("{}", err), "Convergence error: Newton method");
    }

    #[test]
    fn test_with_context() {
        let err = Error::domain("sqrt(-1)").with_context("gamma function");
        assert_eq!(format!("{}", err), "Domain error: gamma function: sqrt(-1)");
    }

    #[test]
    fn test_new_variants() {
        let err = Error::NotImplemented("matrix sqrt".into());
        assert_eq!(format!("{}", err), "Not implemented: matrix sqrt");

        let err = Error::ConvergenceError("max iterations exceeded".into());
        assert_eq!(
            format!("{}", err),
            "Convergence error: max iterations exceeded"
        );
    }

    #[test]
    fn test_severity_classification() {
        assert_eq!(Error::DivisionByZero.severity(), ErrorSeverity::Fatal);
        assert_eq!(Error::MemoryError.severity(), ErrorSeverity::Fatal);
        assert_eq!(Error::Underflow.severity(), ErrorSeverity::Warning);
        assert_eq!(Error::RoundingError.severity(), ErrorSeverity::Warning);
        assert_eq!(
            Error::ConvergenceError("test".into()).severity(),
            ErrorSeverity::Error
        );
        assert_eq!(
            Error::DomainError("test".into()).severity(),
            ErrorSeverity::Error
        );
    }

    #[test]
    fn test_recovery_strategy_success_first_try() {
        let strategy = RecoveryStrategy::default();
        let mut calls = 0;
        let result = strategy.retry(|| {
            calls += 1;
            Ok::<_, Error>(42)
        });
        assert_eq!(result, Ok(42));
        assert_eq!(calls, 1);
    }

    #[test]
    fn test_recovery_strategy_retries_on_error() {
        let strategy = RecoveryStrategy::new(3, 1.0);
        let mut calls = 0;
        let result: Result<i32> = strategy.retry(|| {
            calls += 1;
            if calls < 3 {
                Err(Error::ConvergenceError("not yet".into()))
            } else {
                Ok(99)
            }
        });
        assert_eq!(result, Ok(99));
        assert_eq!(calls, 3);
    }

    #[test]
    fn test_recovery_strategy_fatal_stops_immediately() {
        let strategy = RecoveryStrategy::new(5, 1.0);
        let mut calls = 0;
        let result: Result<i32> = strategy.retry(|| {
            calls += 1;
            Err(Error::DivisionByZero)
        });
        assert!(result.is_err());
        assert_eq!(calls, 1);
    }

    #[test]
    fn test_recovery_strategy_exhausted_returns_last_error() {
        let strategy = RecoveryStrategy::new(2, 1.0);
        let result: Result<i32> =
            strategy.retry(|| Err(Error::ConvergenceError("always fails".into())));
        assert!(result.is_err());
        match result {
            Err(Error::ConvergenceError(msg)) => assert!(msg.contains("always fails")),
            _ => panic!("expected ConvergenceError"),
        }
    }

    #[test]
    fn test_recovery_strategy_no_recovery_fails_immediately() {
        let strategy = RecoveryStrategy::no_recovery();
        let mut calls = 0;
        let result: Result<i32> = strategy.retry(|| {
            calls += 1;
            Err(Error::ConvergenceError("fail".into()))
        });
        assert!(result.is_err());
        assert_eq!(calls, 1); // no retry
    }

    #[test]
    fn test_recovery_strategy_warning_triggers_retry() {
        // Warning errors should trigger retry, not short-circuit (only Fatal short-circuits)
        let strategy = RecoveryStrategy::new(3, 1.0);
        let mut calls = 0;
        let result: Result<i32> = strategy.retry(|| {
            calls += 1;
            if calls < 2 {
                Err(Error::Underflow)
            } else {
                Ok(7)
            }
        });
        assert_eq!(result, Ok(7));
        assert_eq!(calls, 2); // retried after Warning
    }

    #[test]
    fn test_retry_with_precision_advances_and_succeeds() {
        use crate::mpf::adaptive::AdaptivePrecision;
        let strategy = RecoveryStrategy::new(4, 2.0);
        let mut adaptive = AdaptivePrecision::new(64, 512, 2.0);
        let mut calls = 0;
        let result: Result<i32> = strategy.retry_with_precision(
            |precision| {
                calls += 1;
                if precision < 256 {
                    Err(Error::ConvergenceError("need more precision".into()))
                } else {
                    Ok(precision as i32)
                }
            },
            &mut adaptive,
        );
        assert!(result.is_ok());
        assert!(
            calls >= 2,
            "should have retried at least once after precision advance"
        );
        assert!(
            adaptive.current_precision() > 64,
            "precision should have advanced"
        );
    }
}
