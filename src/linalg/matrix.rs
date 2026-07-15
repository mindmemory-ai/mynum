//! 密集矩阵类型和基本运算

use crate::error::{Error, Result};
use crate::mpf::Mpf;

/// 密集矩阵（以 Mpf 为元素）
#[derive(Debug, Clone)]
pub struct Matrix {
    pub(crate) rows: usize,
    pub(crate) cols: usize,
    pub(crate) data: Vec<Mpf>,
}

impl Matrix {
    /// 创建 rows × cols 的零矩阵
    pub fn zero(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![Mpf::new(); rows * cols],
        }
    }

    /// 从 Vec<Vec<f64>> 创建矩阵
    pub fn from_f64_vec(data: Vec<Vec<f64>>, precision: usize) -> Self {
        let rows = data.len();
        let cols = data.first().map_or(0, |r| r.len());
        let mut flat = Vec::with_capacity(rows * cols);
        for row in &data {
            for &val in row {
                flat.push(Mpf::from_f64(val, precision));
            }
        }
        Self {
            rows,
            cols,
            data: flat,
        }
    }

    /// 获取元素 (i, j)
    pub fn get(&self, i: usize, j: usize) -> Option<&Mpf> {
        if i < self.rows && j < self.cols {
            Some(&self.data[i * self.cols + j])
        } else {
            None
        }
    }

    /// 设置元素 (i, j)
    pub fn set(&mut self, i: usize, j: usize, val: Mpf) -> Result<()> {
        if i >= self.rows || j >= self.cols {
            return Err(Error::invalid_input("index out of bounds"));
        }
        self.data[i * self.cols + j] = val;
        Ok(())
    }

    /// 矩阵加法
    pub fn add(&self, other: &Matrix) -> Result<Matrix> {
        if self.rows != other.rows || self.cols != other.cols {
            return Err(Error::invalid_input("matrix dimensions must match"));
        }
        let mut result = self.clone();
        for i in 0..result.data.len() {
            result.data[i] = result.data[i].add(&other.data[i]);
        }
        Ok(result)
    }

    /// 矩阵减法
    pub fn sub(&self, other: &Matrix) -> Result<Matrix> {
        if self.rows != other.rows || self.cols != other.cols {
            return Err(Error::invalid_input("matrix dimensions must match"));
        }
        let mut result = self.clone();
        for i in 0..result.data.len() {
            result.data[i] = result.data[i].sub(&other.data[i]);
        }
        Ok(result)
    }

    /// 矩阵乘法
    pub fn mul(&self, other: &Matrix) -> Result<Matrix> {
        if self.cols != other.rows {
            return Err(Error::invalid_input("matrix dimensions incompatible"));
        }
        let mut result = Matrix::zero(self.rows, other.cols);
        for i in 0..self.rows {
            for k in 0..self.cols {
                let aik = &self.data[i * self.cols + k];
                if aik.is_zero() {
                    continue;
                }
                for j in 0..other.cols {
                    let bkj = &other.data[k * other.cols + j];
                    if bkj.is_zero() {
                        continue;
                    }
                    let prod = aik.mul(bkj);
                    let idx = i * result.cols + j;
                    result.data[idx] = result.data[idx].add(&prod);
                }
            }
        }
        Ok(result)
    }

    /// 标量乘法
    pub fn scale(&self, scalar: &Mpf) -> Matrix {
        let mut result = self.clone();
        for val in &mut result.data {
            *val = val.mul(scalar);
        }
        result
    }

    /// 转置
    pub fn transpose(&self) -> Matrix {
        let mut result = Matrix::zero(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                result.data[j * self.rows + i] = self.data[i * self.cols + j].clone();
            }
        }
        result
    }

    /// 单位矩阵
    pub fn identity(n: usize, precision: usize) -> Self {
        let mut result = Matrix::zero(n, n);
        let one = Mpf::from_i64(1, precision);
        for i in 0..n {
            result.data[i * n + i] = one.clone();
        }
        result
    }

    /// Return the number of rows.
    pub fn rows(&self) -> usize {
        self.rows
    }
    /// Return the number of columns.
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Frobenius norm: sqrt(sum of squares of all elements)
    pub fn norm_frobenius(&self) -> Result<Mpf> {
        let mut sum_sq = Mpf::new();
        for val in &self.data {
            sum_sq = sum_sq.add(&val.mul(val));
        }
        sum_sq.sqrt()
    }

    /// Check if A == A^T (square and symmetric)
    pub fn is_symmetric(&self) -> bool {
        if self.rows != self.cols {
            return false;
        }
        for i in 0..self.rows {
            for j in 0..i {
                if self.get(i, j).unwrap().cmp(self.get(j, i).unwrap())
                    != core::cmp::Ordering::Equal
                {
                    return false;
                }
            }
        }
        true
    }

    /// Extract diagonal elements
    pub fn diag(&self) -> Vec<Mpf> {
        let n = self.rows.min(self.cols);
        (0..n)
            .map(|i| self.data[i * self.cols + i].clone())
            .collect()
    }

    /// Create diagonal matrix from diagonal elements
    pub fn from_diag(diag: &[Mpf], _precision: usize) -> Self {
        let n = diag.len();
        let mut m = Matrix::zero(n, n);
        for (i, d) in diag.iter().enumerate() {
            m.data[i * n + i] = d.clone();
        }
        m
    }

    /// Sum of diagonal elements
    pub fn trace(&self) -> Mpf {
        let n = self.rows.min(self.cols);
        let mut sum = Mpf::new();
        for i in 0..n {
            sum = sum.add(&self.data[i * self.cols + i]);
        }
        sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_matrix_add() {
        let a = Matrix::from_f64_vec(vec![vec![1.0, 2.0], vec![3.0, 4.0]], 64);
        let b = Matrix::from_f64_vec(vec![vec![5.0, 6.0], vec![7.0, 8.0]], 64);
        let c = a.add(&b).unwrap();
        assert!((c.get(0, 0).unwrap().to_f64().unwrap() - 6.0).abs() < 1e-10);
        assert!((c.get(1, 1).unwrap().to_f64().unwrap() - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_matrix_mul() {
        let a = Matrix::from_f64_vec(vec![vec![1.0, 2.0], vec![3.0, 4.0]], 64);
        let b = Matrix::from_f64_vec(vec![vec![2.0, 0.0], vec![1.0, 2.0]], 64);
        let c = a.mul(&b).unwrap();
        assert!((c.get(0, 0).unwrap().to_f64().unwrap() - 4.0).abs() < 1e-10);
        assert!((c.get(0, 1).unwrap().to_f64().unwrap() - 4.0).abs() < 1e-10);
        assert!((c.get(1, 0).unwrap().to_f64().unwrap() - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_transpose() {
        let a = Matrix::from_f64_vec(vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]], 64);
        let t = a.transpose();
        assert_eq!(t.rows(), 3);
        assert_eq!(t.cols(), 2);
        assert!((t.get(0, 1).unwrap().to_f64().unwrap() - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_norm_frobenius() {
        let a = Matrix::from_f64_vec(vec![vec![3.0, 4.0], vec![0.0, 0.0]], 64);
        let norm = a.norm_frobenius().unwrap();
        assert!((norm.to_f64().unwrap() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_is_symmetric() {
        let a = Matrix::from_f64_vec(vec![vec![1.0, 2.0], vec![2.0, 3.0]], 64);
        assert!(a.is_symmetric());
        let b = Matrix::from_f64_vec(vec![vec![1.0, 2.0], vec![3.0, 4.0]], 64);
        assert!(!b.is_symmetric());
    }

    #[test]
    fn test_diag() {
        let a = Matrix::from_f64_vec(vec![vec![1.0, 2.0], vec![3.0, 4.0]], 64);
        let d = a.diag();
        assert_eq!(d.len(), 2);
        assert!((d[0].to_f64().unwrap() - 1.0).abs() < 1e-10);
        assert!((d[1].to_f64().unwrap() - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_trace() {
        let a = Matrix::from_f64_vec(vec![vec![1.0, 2.0], vec![3.0, 4.0]], 64);
        let t = a.trace();
        assert!((t.to_f64().unwrap() - 5.0).abs() < 1e-10);
    }
}
