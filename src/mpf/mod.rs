//! Mpf 模块 - 高精度浮点数运算
//!
//! 基于 Mpz 实现的高精度浮点数运算，支持任意精度的浮点数计算。

pub mod adaptive;
pub mod arithmetic;
pub mod comparison;
pub mod config;
pub mod constants;
pub mod conversion;
pub mod cordic;
pub mod core;
pub mod differentiation;
pub mod elementary;
pub mod integration;
pub mod random;
pub mod solve;
pub mod special;
pub mod traits;

pub use config::{MpfConfig, MpfPrecisionConfig, MpfRoundingMode};
pub use cordic::{Cordic, CordicConfig};
pub use core::Mpf;
