//! MyNum - 高精度数值计算库
//!
//! 这是一个工业级别的高精度数值计算库，支持：
//! - [`mpz`] - 大整数运算（包含并行乘法支持）
//! - [`mpf`] - 高精度浮点数运算
//! - [`complex`] - 高精度复数运算
//! - [`linalg`] - 线性代数（矩阵运算、LU 分解）
//! - [`signal`] - 信号处理（FFT/iFFT、卷积、自相关）
//! - 全局配置管理

#![warn(unused_imports)]
#![warn(unused_variables)]

pub mod algorithm;
pub mod complex;
pub mod config;
pub mod error;
pub mod linalg;
pub mod mpf;
pub mod mpz;
pub mod profile;
pub mod signal;

pub use complex::{Complex, ComplexConfig, ComplexPrecisionConfig};
pub use config::GlobalPrecisionConfig;
pub use error::{Error, ErrorSeverity, RecoveryStrategy, Result};
pub use mpf::{Mpf, MpfConfig, MpfRoundingMode};
pub use mpz::{AlgorithmThresholds, Mpz, MpzMultiplicationConfig, MultiplicationBackend};

// 基础编译检查
#[test]
fn it_works() {
    assert!(true); // 基础编译检查
}
