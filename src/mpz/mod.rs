//! 大整数运算模块
//!
//! 提供任意精度的整数运算功能。

pub mod arithmetic;
pub mod bitwise;
pub mod comparison;
pub mod config;
pub mod conversion;
pub mod core;
pub mod multiplication;
pub mod number_theory;
pub mod parallel;
pub mod pool;
pub mod random;
pub mod traits;

pub use config::{AlgorithmThresholds, MpzMultiplicationConfig, MultiplicationBackend};
pub use core::Mpz;
pub use parallel::{
    get_global_parallel, set_global_parallel_config, ParallelConfig, ParallelMultiplier,
};

// 重新导出主要功能（暂时注释掉未使用的导入）
// pub use arithmetic::*;
// pub use multiplication::*;
// pub use comparison::*;
// pub use number_theory::*;
// pub use random::*;
