//! Mpf 实轴数值积分与极限模块
//!
//! 提供实值函数在实数轴上的数值积分算法和极限计算。

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

use crate::error::{Error, Result};
use crate::mpf::Mpf;
use crate::mpz::Mpz;

/// 梯形法则数值积分
///
/// ∫_a^b f(x) dx ≈ (h/2) * [f(a) + 2*Σf(x_i) + f(b)]
pub fn trapezoidal<F>(f: &F, a: &Mpf, b: &Mpf, n: usize) -> Result<Mpf>
where
    F: Fn(&Mpf) -> Result<Mpf>,
{
    if n == 0 {
        return Err(Error::invalid_input("n must be positive"));
    }
    let n_mpf = Mpf::from_mpz(Mpz::from_i64(n as i64), a.precision());
    let h = b.sub(a).div(&n_mpf)?;

    let mut sum = f(a)?.add(&f(b)?);
    let two = Mpf::from_mpz(Mpz::from_i64(2), a.precision());

    for i in 1..n {
        let i_mpf = Mpf::from_mpz(Mpz::from_i64(i as i64), a.precision());
        let x = a.add(&h.mul(&i_mpf));
        sum = sum.add(&f(&x)?.mul(&two));
    }

    h.mul(&sum).div(&two)
}

/// 辛普森法则数值积分
///
/// ∫_a^b f(x) dx ≈ (h/3) * [f(a) + 4f(x_1) + 2f(x_2) + ... + f(b)]
/// n 必须是偶数。
pub fn simpson<F>(f: &F, a: &Mpf, b: &Mpf, n: usize) -> Result<Mpf>
where
    F: Fn(&Mpf) -> Result<Mpf>,
{
    if n == 0 || !n.is_multiple_of(2) {
        return Err(Error::invalid_input("n must be a positive even number"));
    }
    let n_mpf = Mpf::from_mpz(Mpz::from_i64(n as i64), a.precision());
    let h = b.sub(a).div(&n_mpf)?;

    let mut sum = f(a)?.add(&f(b)?);
    let two = Mpf::from_mpz(Mpz::from_i64(2), a.precision());
    let four = Mpf::from_mpz(Mpz::from_i64(4), a.precision());

    for i in 1..n {
        let i_mpf = Mpf::from_mpz(Mpz::from_i64(i as i64), a.precision());
        let x = a.add(&h.mul(&i_mpf));
        let coeff = if i % 2 == 0 { &two } else { &four };
        sum = sum.add(&f(&x)?.mul(coeff));
    }

    let three = Mpf::from_mpz(Mpz::from_i64(3), a.precision());
    h.mul(&sum).div(&three)
}

/// 五阶高斯-勒让德求积
///
/// 固定的 5 点规则，通过变量变换适应 [a, b]。
pub fn gauss_legendre_5<F>(f: &F, a: &Mpf, b: &Mpf) -> Result<Mpf>
where
    F: Fn(&Mpf) -> Result<Mpf>,
{
    let nodes: [f64; 5] = [
        -0.906_179_845_938_664,
        -0.5384693101056831,
        0.0,
        0.5384693101056831,
        0.906_179_845_938_664,
    ];
    let weights: [f64; 5] = [
        0.2369268850561891,
        0.4786286704993665,
        0.5688888888888889,
        0.4786286704993665,
        0.2369268850561891,
    ];

    let prec = a.precision();
    let half = Mpf::from_f64(0.5, prec);
    let b_minus_a = b.sub(a);
    let b_plus_a = b.add(a);

    let mut sum = Mpf::new();
    for i in 0..5 {
        let node = Mpf::from_f64(nodes[i], prec);
        let weight = Mpf::from_f64(weights[i], prec);
        let x = half.mul(&b_minus_a).mul(&node).add(&half.mul(&b_plus_a));
        sum = sum.add(&f(&x)?.mul(&weight));
    }

    Ok(half.mul(&b_minus_a).mul(&sum))
}

/// Adaptive Gauss-Kronrod (G7-K15) quadrature.
///
/// Recursively subdivides [a,b] until the difference between the 7-point
/// Gauss and 15-point Kronrod estimates is below `tol`.
///
/// Default tol: 1e-14 (at current precision). Default max_depth: 20.
///
/// # Example
///
/// ```
/// use mynum::Mpf;
/// use mynum::mpf::integration::quad;
/// let a = Mpf::new();
/// let b = Mpf::from_i64(1, 64);
/// let f = |x: &Mpf| Ok(x.mul(x));
/// let result = quad(&f, &a, &b, None, None).unwrap();
/// ```
pub fn quad<F>(f: &F, a: &Mpf, b: &Mpf, tol: Option<&Mpf>, max_depth: Option<usize>) -> Result<Mpf>
where
    F: Fn(&Mpf) -> Result<Mpf>,
{
    let precision = a.precision();
    let default_tol = Mpf::from_f64(1e-14, precision);
    let tol_ref = tol.unwrap_or(&default_tol);
    let depth = max_depth.unwrap_or(20);

    quad_recursive(f, a, b, tol_ref, depth)
}

/// Recursive subdivision for adaptive Gauss-Kronrod quadrature.
///
/// Evaluates the G7-K15 pair on [a,b]. If the error (|K15 - G7|) is below
/// `tol` or `remaining_depth` has been exhausted, returns the K15 estimate.
/// Otherwise splits the interval at the midpoint and recurses on each half.
fn quad_recursive<F>(f: &F, a: &Mpf, b: &Mpf, tol: &Mpf, remaining_depth: usize) -> Result<Mpf>
where
    F: Fn(&Mpf) -> Result<Mpf>,
{
    let precision = a.precision();

    // G7-K15 nodes (positive half, symmetric about 0).
    // G7 nodes are a subset at indices 0, 2, 4, 6.
    // K15 adds nodes at indices 1, 3, 5, 7.
    let nodes: [f64; 8] = [
        0.0,
        0.20778495500789848,
        0.4058451513773972,
        0.5860872354676911,
        0.7415311855993945,
        0.8648644233597691,
        0.9491079123427585,
        0.9914553711208126,
    ];
    // K15 weights for interval [-1, 1]
    let weights_k15: [f64; 8] = [
        0.20948214108472782,
        0.20443294007529889,
        0.19035057806478542,
        0.169_004_726_639_267_9,
        0.14065325971552592,
        0.10479001032225019,
        0.06309209262997855,
        0.022935322010529224,
    ];
    // G7 weights for interval [-1, 1]
    // Order matches G7 nodes at indices 0, 2, 4, 6 of nodes[]
    let weights_g7: [f64; 4] = [
        0.4179591836734694,  // node 0.0
        0.3818300505051189,  // node 0.4058451513773972
        0.27970539148927664, // node 0.7415311855993945
        0.1294849661688697,  // node 0.9491079123427585
    ];

    // G7 node index -> G7 weight index mapping.
    // Indices 1, 3, 5, 7 are K15-only (not part of G7).
    const G7_IDX: [usize; 8] = [
        0,          // index 0: G7 weight 0
        usize::MAX, // index 1: K15-only
        1,          // index 2: G7 weight 1
        usize::MAX, // index 3: K15-only
        2,          // index 4: G7 weight 2
        usize::MAX, // index 5: K15-only
        3,          // index 6: G7 weight 3
        usize::MAX, // index 7: K15-only
    ];

    let half = Mpf::from_f64(0.5, precision);
    let b_minus_a_half = b.sub(a).mul(&half);
    let b_plus_a_half = b.add(a).mul(&half);

    let mut sum_k15 = Mpf::new();
    let mut sum_g7 = Mpf::new();

    // Center node (index 0) — no symmetric counterpart.
    let x0 = b_plus_a_half.clone();
    let f0 = f(&x0)?;
    sum_k15 = sum_k15.add(&f0.mul(&Mpf::from_f64(weights_k15[0], precision)));
    sum_g7 = sum_g7.add(&f0.mul(&Mpf::from_f64(weights_g7[0], precision)));

    // Remaining positive nodes (indices 1..8) — symmetric ±x.
    for i in 1..8 {
        let node = Mpf::from_f64(nodes[i], precision);
        let x_plus = b_plus_a_half.add(&b_minus_a_half.mul(&node));
        let x_minus = b_plus_a_half.sub(&b_minus_a_half.mul(&node));

        let f_sum = f(&x_plus)?.add(&f(&x_minus)?);

        // K15 contribution
        let wk15 = Mpf::from_f64(weights_k15[i], precision);
        sum_k15 = sum_k15.add(&f_sum.mul(&wk15));

        // G7 contribution (only for G7 nodes)
        let g7_idx = G7_IDX[i];
        if g7_idx != usize::MAX {
            let wg7 = Mpf::from_f64(weights_g7[g7_idx], precision);
            sum_g7 = sum_g7.add(&f_sum.mul(&wg7));
        }
    }

    let result_k15 = b_minus_a_half.mul(&sum_k15);
    let result_g7 = b_minus_a_half.mul(&sum_g7);
    let error = result_k15.sub(&result_g7).abs();

    if error.cmp(tol) == core::cmp::Ordering::Less || remaining_depth == 0 {
        Ok(result_k15)
    } else {
        let mid = a.add(b).div(&Mpf::from_mpz(Mpz::from_i64(2), precision))?;
        let left = quad_recursive(f, a, &mid, tol, remaining_depth - 1)?;
        let right = quad_recursive(f, &mid, b, tol, remaining_depth - 1)?;
        Ok(left.add(&right))
    }
}

/// Sum a converging series with convergence acceleration (Wynn's epsilon algorithm).
///
/// `f(n)` returns the n-th term (n starts at `start`). Terms are summed until
/// `|term| < tol` or `max_terms` terms are evaluated, then Wynn's epsilon
/// transformation is applied to the partial sum sequence to accelerate convergence.
///
/// # Example
///
/// ```
/// use mynum::Mpf;
/// use mynum::mpf::integration::nsum;
/// use mynum::error::Result;
///
/// // Sum Σ 1/2^n from n=0 to n=10 (approaches 2)
/// let f = |n: u64| -> Result<Mpf> {
///     let prec = 64;
///     Ok(Mpf::from_i64(1, prec).div(&Mpf::from_i64(1i64 << (n as u32), prec))?)
/// };
/// let result = nsum(&f, 0, 10, None).unwrap();
/// assert!((result.to_f64().unwrap() - 2.0).abs() < 0.01);
/// ```
pub fn nsum<F>(f: &F, start: u64, max_terms: usize, tol: Option<&Mpf>) -> Result<Mpf>
where
    F: Fn(u64) -> Result<Mpf>,
{
    if max_terms == 0 {
        return Ok(Mpf::new());
    }
    let precision = tol.map(|t| t.precision()).unwrap_or(64);
    let threshold = tol.cloned().unwrap_or(Mpf::from_f64(1e-14, precision));

    let mut partial_sums: Vec<Mpf> = Vec::new();
    let mut sum = Mpf::new();

    for n in start..start + max_terms as u64 {
        let term = f(n)?;
        sum = sum.add(&term);
        partial_sums.push(sum.clone());
        if term.abs().cmp(&threshold) == core::cmp::Ordering::Less && n > start + 1 {
            break;
        }
    }

    // Wynn's epsilon algorithm for convergence acceleration
    // ε^{(0)}_k = S_k, ε^{(1)}_k = 1/(S_{k+1} - S_k), ε^{(m+1)}_k = ε^{(m-1)}_{k+1} + 1/(ε^{(m)}_{k+1} - ε^{(m)}_k)
    let n = partial_sums.len();
    if n < 3 {
        return Ok(partial_sums.last().unwrap().clone());
    }

    let mut e: Vec<Vec<Option<Mpf>>> = vec![vec![None; n]; n];
    for k in 0..n {
        e[0][k] = Some(partial_sums[k].clone());
    }

    for m in 1..n {
        for k in 0..n - m {
            let d = e[m - 1][k + 1]
                .as_ref()
                .unwrap()
                .sub(e[m - 1][k].as_ref().unwrap());
            if d.is_zero() {
                e[m][k] = e[m - 1][k + 1].clone();
            } else {
                let inv = Mpf::from_i64(1, precision).div(&d)?;
                e[m][k] = Some(if m >= 2 {
                    e[m - 2][k + 1].as_ref().unwrap().add(&inv)
                } else {
                    inv
                });
            }
        }
    }

    // Return the best accelerated estimate (deepest even-order epsilon value)
    for m in (2..n).step_by(2).rev() {
        if let Some(ref val) = e[m][0] {
            return Ok(val.clone());
        }
    }
    Ok(partial_sums.last().unwrap().clone())
}

/// 使用 Richardson 外推法计算数值极限
///
/// 计算当 x → x0 时 f(x) 的极限。通过序列 f(x0 + h_n) 进行
/// Richardson 外推，其中 h_n = 2^{-n} * direction，逐步消除
/// O(h), O(h^2), O(h^3) 等误差项。
///
/// `direction` 默认为 -1.0（从左侧逼近），也可设为 +1.0（从右侧逼近）。
///
/// # Example
///
/// ```
/// use mynum::Mpf;
/// use mynum::mpf::integration::limit;
/// let x0 = Mpf::new();
/// let f = |x: &Mpf| x.sin()?.div(x);
/// let result = limit(&f, &x0, None).unwrap();
/// assert!((result.to_f64().unwrap() - 1.0).abs() < 1e-8);
/// ```
pub fn limit<F>(f: &F, x0: &Mpf, direction: Option<f64>) -> Result<Mpf>
where
    F: Fn(&Mpf) -> Result<Mpf>,
{
    let precision = x0.precision();
    let dir = direction.unwrap_or(-1.0);
    let initial_h = Mpf::from_f64(dir * 0.5, precision);

    // Richardson extrapolation table
    // D_{n,0} = f(x0 + h/2^n)
    // D_{n,k} = D_{n,k-1} + (D_{n,k-1} - D_{n-1,k-1}) / (2^k - 1)
    let max_rows = 10;
    let mut table: Vec<Vec<Mpf>> = Vec::new();
    let two = Mpf::from_mpz(Mpz::from_i64(2), precision);
    let one = Mpf::from_mpz(Mpz::from_i64(1), precision);

    for n in 0..max_rows {
        let h = initial_h.div(&two.pow(n as u32)?)?;
        let x = x0.add(&h);
        let val = f(&x)?;

        let mut row = vec![val];
        if !table.is_empty() {
            for k in 1..=n {
                let two_pow_k = Mpf::from_mpz(Mpz::from_i64(1i64 << k), precision);
                let denom = two_pow_k.sub(&one);
                let d_n_k_minus_1 = &row[k - 1];
                let d_n_minus_1_k_minus_1 = &table[n - 1][k - 1];
                let delta = d_n_k_minus_1.sub(d_n_minus_1_k_minus_1);
                let d_n_k = d_n_k_minus_1.add(&delta.div(&denom)?);
                row.push(d_n_k);
            }
        }

        table.push(row);

        // Check convergence: compare the last two rows' best estimates
        if n >= 3 {
            let curr = &table[n][n.min(table[n].len() - 1)];
            let prev = &table[n - 1][(n - 1).min(table[n - 1].len() - 1)];
            let diff = curr.sub(prev).abs();
            if diff.cmp(&Mpf::from_f64(1e-14, precision)) == core::cmp::Ordering::Less {
                return Ok(curr.clone());
            }
        }
    }

    let last_row = table.last().unwrap();
    Ok(last_row[last_row.len() - 1].clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mpz::Mpz;

    #[test]
    fn test_trapezoidal_linear() {
        let a = Mpf::new();
        let b = Mpf::from_i64(1, 64);
        let f = |x: &Mpf| Ok(x.clone());
        let result = trapezoidal(&f, &a, &b, 100).unwrap();
        assert!((result.to_f64().unwrap() - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_simpson_quadratic() {
        let a = Mpf::new();
        let b = Mpf::from_i64(1, 64);
        let f = |x: &Mpf| Ok(x.mul(x));
        let result = simpson(&f, &a, &b, 100).unwrap();
        assert!((result.to_f64().unwrap() - 1.0 / 3.0).abs() < 0.001);
    }

    #[test]
    fn test_gauss_legendre_polynomial() {
        let a = Mpf::from_i64(-1, 64);
        let b = Mpf::from_i64(1, 64);
        let f = |x: &Mpf| {
            let x2 = x.mul(x);
            Ok(x2.mul(&x2))
        };
        let result = gauss_legendre_5(&f, &a, &b).unwrap();
        assert!((result.to_f64().unwrap() - 0.4).abs() < 1e-10);
    }

    #[test]
    fn test_quad_x_squared() {
        let a = Mpf::new();
        let b = Mpf::from_i64(1, 64);
        let f = |x: &Mpf| Ok(x.mul(x));
        let result = quad(&f, &a, &b, None, None).unwrap();
        assert!((result.to_f64().unwrap() - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_quad_sin() {
        let a = Mpf::new();
        let b = Mpf::pi(64);
        let f = |x: &Mpf| x.sin();
        let result = quad(&f, &a, &b, None, None).unwrap();
        assert!((result.to_f64().unwrap() - 2.0).abs() < 1e-8);
    }

    #[test]
    fn test_quad_invalid_interval() {
        let a = Mpf::from_i64(1, 64);
        let b = Mpf::new();
        let f = |x: &Mpf| Ok(x.clone());
        let result = quad(&f, &a, &b, None, None).unwrap();
        assert!((result.to_f64().unwrap() + 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_nsum_basel() {
        // Σ 1/n^2 = π^2/6 ≈ 1.6449340668
        // Use larger max_terms since epsilon is less effective for monotone series
        let f = |n: u64| -> Result<Mpf> {
            let n_mpf = Mpf::from_mpz(Mpz::from_i64(n as i64), 64);
            Ok(Mpf::from_mpz(Mpz::from_i64(1), 64).div(&n_mpf.mul(&n_mpf))?)
        };
        let result = nsum(&f, 1, 2000, None).unwrap();
        let val = result.to_f64().unwrap();
        let expected = 1.6449340668482264;
        eprintln!(
            "test_nsum_basel: val={:.12} expected={:.12} diff={:.12e}",
            val,
            expected,
            (val - expected).abs()
        );
        assert!((val - expected).abs() < 1e-4);
    }

    #[test]
    fn test_nsum_alternating() {
        // Σ (-1)^{n+1}/n = ln(2) ≈ 0.6931471806
        let f = |n: u64| -> Result<Mpf> {
            let n_mpf = Mpf::from_mpz(Mpz::from_i64(n as i64), 64);
            let term = Mpf::from_mpz(Mpz::from_i64(1), 64).div(&n_mpf)?;
            if n % 2 == 0 {
                Ok(term.neg())
            } else {
                Ok(term)
            }
        };
        let result = nsum(&f, 1, 1000, None).unwrap();
        let val = result.to_f64().unwrap();
        let expected = 0.6931471805599453;
        eprintln!(
            "test_nsum_alternating: val={:.12} expected={:.12} diff={:.12e}",
            val,
            expected,
            (val - expected).abs()
        );
        assert!((val - expected).abs() < 1e-6);
    }

    #[test]
    fn test_limit_sin_over_x() {
        // lim_{x->0} sin(x)/x = 1
        let f = |x: &Mpf| x.sin()?.div(x);
        let x0 = Mpf::new();
        let result = limit(&f, &x0, None).unwrap();
        assert!((result.to_f64().unwrap() - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_limit_exp_minus_1_over_x() {
        // lim_{x->0} (e^x - 1)/x = 1
        let one = Mpf::from_mpz(Mpz::from_i64(1), 128);
        let f = |x: &Mpf| Ok(x.exp()?.sub(&one).div(x)?);
        let result = limit(&f, &Mpf::with_precision(128), None).unwrap();
        assert!((result.to_f64().unwrap() - 1.0).abs() < 1e-8);
    }
}
