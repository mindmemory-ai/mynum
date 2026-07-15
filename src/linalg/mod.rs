//! 线性代数模块
//!
//! 提供基于 Mpf 的矩阵运算和分解。

pub mod decomposition;
pub mod matrix;

pub use decomposition::{
    cholesky, condition_number, eigenvalues, expm, logm, lu_decomposition, qr_decomposition, sqrtm,
    svd, EigenDecomposition, LUDecomposition, QRDecomposition, SVD,
};
