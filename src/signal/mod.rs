//! 信号处理模块
//!
//! FFT/iFFT、卷积、自相关。

use crate::error::{Error, Result};
use crate::mpf::Mpf;
use std::f64::consts::PI;

/// 计算实值序列的 FFT（返回复数结果）
///
/// 输入长度必须是 2 的幂。
pub fn fft_real(data: &[Mpf]) -> Result<Vec<(f64, f64)>> {
    let n = data.len();
    if n == 0 || !n.is_power_of_two() {
        return Err(Error::invalid_input("FFT length must be a power of two"));
    }

    let mut complex: Vec<(f64, f64)> = data
        .iter()
        .map(|x| (x.to_f64().unwrap_or(0.0), 0.0f64))
        .collect();

    // 位逆序重排
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            complex.swap(i, j);
        }
    }

    // Cooley-Tukey 迭代 FFT
    let mut len = 2;
    while len <= n {
        let angle = -2.0 * PI / len as f64;
        let wlen = (angle.cos(), angle.sin());
        for i in (0..n).step_by(len) {
            let mut w = (1.0f64, 0.0f64);
            for j in 0..len / 2 {
                let u = complex[i + j];
                let v = (
                    complex[i + j + len / 2].0 * w.0 - complex[i + j + len / 2].1 * w.1,
                    complex[i + j + len / 2].0 * w.1 + complex[i + j + len / 2].1 * w.0,
                );
                complex[i + j] = (u.0 + v.0, u.1 + v.1);
                complex[i + j + len / 2] = (u.0 - v.0, u.1 - v.1);
                let w_next = (w.0 * wlen.0 - w.1 * wlen.1, w.0 * wlen.1 + w.1 * wlen.0);
                w = w_next;
            }
        }
        len <<= 1;
    }

    Ok(complex)
}

/// 逆 FFT
pub fn ifft_real(spectrum: &[(f64, f64)]) -> Result<Vec<Mpf>> {
    let n = spectrum.len();
    if n == 0 || !n.is_power_of_two() {
        return Err(Error::invalid_input("IFFT length must be a power of two"));
    }

    // 共轭
    let mut complex: Vec<(f64, f64)> = spectrum.iter().map(|(re, im)| (*re, -im)).collect();

    // FFT（与 fft_real 相同的 Cooley-Tukey，角度符号相反）
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            complex.swap(i, j);
        }
    }

    let mut len = 2;
    while len <= n {
        let angle = -2.0 * PI / len as f64;
        let wlen = (angle.cos(), angle.sin());
        for i in (0..n).step_by(len) {
            let mut w = (1.0f64, 0.0f64);
            for j in 0..len / 2 {
                let u = complex[i + j];
                let v = (
                    complex[i + j + len / 2].0 * w.0 - complex[i + j + len / 2].1 * w.1,
                    complex[i + j + len / 2].0 * w.1 + complex[i + j + len / 2].1 * w.0,
                );
                complex[i + j] = (u.0 + v.0, u.1 + v.1);
                complex[i + j + len / 2] = (u.0 - v.0, u.1 - v.1);
                let w_next = (w.0 * wlen.0 - w.1 * wlen.1, w.0 * wlen.1 + w.1 * wlen.0);
                w = w_next;
            }
        }
        len <<= 1;
    }

    // 缩放 + 取实部
    let scale = 1.0 / n as f64;
    Ok(complex
        .iter()
        .map(|(re, _)| Mpf::from_f64(re * scale, 64))
        .collect())
}

/// 直接卷积: result[k] = Σ a[i] * b[k-i]
pub fn convolution_direct(a: &[Mpf], b: &[Mpf]) -> Vec<Mpf> {
    let n = a.len() + b.len() - 1;
    let mut result = vec![Mpf::new(); n];
    for (i, ai) in a.iter().enumerate() {
        for (j, bj) in b.iter().enumerate() {
            result[i + j] = result[i + j].add(&ai.mul(bj));
        }
    }
    result
}

/// 自相关: R[k] = Σ x[i] * x[i+k]
pub fn autocorrelation(data: &[Mpf]) -> Vec<Mpf> {
    let n = data.len();
    let mut result = vec![Mpf::new(); n];
    for k in 0..n {
        let mut sum = Mpf::new();
        for i in 0..n - k {
            sum = sum.add(&data[i].mul(&data[i + k]));
        }
        result[k] = sum;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fft_ifft_roundtrip() {
        let data: Vec<Mpf> = (0..8).map(|i| Mpf::from_f64(i as f64, 64)).collect();
        let spectrum = fft_real(&data).unwrap();
        let recovered = ifft_real(&spectrum).unwrap();
        for (orig, rec) in data.iter().zip(recovered.iter()) {
            let diff = (orig.to_f64().unwrap() - rec.to_f64().unwrap()).abs();
            assert!(
                diff < 1e-10,
                "FFT roundtrip error: {} vs {}",
                orig.to_f64().unwrap(),
                rec.to_f64().unwrap()
            );
        }
    }

    #[test]
    fn test_convolution() {
        let a: Vec<Mpf> = vec![1.0, 2.0, 3.0]
            .iter()
            .map(|&x| Mpf::from_f64(x, 64))
            .collect();
        let b: Vec<Mpf> = vec![4.0, 5.0]
            .iter()
            .map(|&x| Mpf::from_f64(x, 64))
            .collect();
        let result = convolution_direct(&a, &b);
        // [4, 13, 22, 15]
        assert!((result[0].to_f64().unwrap() - 4.0).abs() < 1e-10);
        assert!((result[1].to_f64().unwrap() - 13.0).abs() < 1e-10);
        assert!((result[2].to_f64().unwrap() - 22.0).abs() < 1e-10);
        assert!((result[3].to_f64().unwrap() - 15.0).abs() < 1e-10);
    }
}
