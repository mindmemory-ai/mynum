//! 高精度复数运算模块
//!
//! 这是一个工业级别的复数运算库，支持：
//! - 高精度复数算术运算
//! - 复数三角函数和双曲函数
//! - 复数指数和对数函数
//! - 复数通用数学函数
//! - 复数配置管理

pub mod arithmetic;
pub mod config;
pub mod core;
pub mod exponential;
pub mod functions;
pub mod hyperbolic;
pub mod integration;
pub mod polynomial;
pub mod trigonometry;

pub use config::{
    ComplexCacheConfig, ComplexConfig, ComplexParallelConfig, ComplexPrecisionConfig,
};
pub use core::Complex;

/// 复数运算的常量
pub mod constants {
    use super::core::Complex;
    use crate::mpf::Mpf;

    /// 复数单位 i
    pub fn i() -> Complex {
        Complex::from_imag(Mpf::from_i64(1, 64))
    }

    /// 复数单位 -i
    pub fn neg_i() -> Complex {
        Complex::from_imag(Mpf::from_i64(-1, 64))
    }

    /// 复数 1
    pub fn one() -> Complex {
        Complex::from_real(Mpf::from_i64(1, 64))
    }

    /// 复数 -1
    pub fn neg_one() -> Complex {
        Complex::from_real(Mpf::from_i64(-1, 64))
    }

    /// 复数 0
    pub fn zero() -> Complex {
        Complex::new()
    }

    /// 复数 π
    pub fn pi(precision: usize) -> Complex {
        Complex::from_real(Mpf::pi(precision))
    }

    /// 复数 e
    pub fn e(precision: usize) -> Complex {
        Complex::from_real(Mpf::e(precision))
    }

    /// 复数 2π
    pub fn two_pi(precision: usize) -> Complex {
        let pi_val = Mpf::pi(precision);
        let two = Mpf::from_i64(2, precision);
        Complex::from_real(pi_val.mul(&two))
    }

    /// 复数 π/2
    pub fn pi_over_2(precision: usize) -> Complex {
        let pi_val = Mpf::pi(precision);
        let two = Mpf::from_i64(2, precision);
        Complex::from_real(pi_val.div(&two).unwrap())
    }

    /// 复数 π/4
    pub fn pi_over_4(precision: usize) -> Complex {
        let pi_val = Mpf::pi(precision);
        let four = Mpf::from_i64(4, precision);
        Complex::from_real(pi_val.div(&four).unwrap())
    }
}
