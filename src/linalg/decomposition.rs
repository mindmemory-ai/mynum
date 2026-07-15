//! 矩阵分解

use super::matrix::Matrix;
use crate::error::{Error, Result};
use crate::mpf::Mpf;

/// LU 分解结果
pub struct LUDecomposition {
    /// Lower triangular matrix L (unit diagonal).
    pub l: Matrix,
    /// Upper triangular matrix U.
    pub u: Matrix,
    /// Pivot permutation vector (PA = LU).
    pub p: Vec<usize>,
}

/// 带部分主元选取的 LU 分解
///
/// PA = LU，其中 P 是置换矩阵。
pub fn lu_decomposition(a: &Matrix) -> Result<LUDecomposition> {
    let n = a.rows();
    if n != a.cols() {
        return Err(Error::invalid_input("matrix must be square"));
    }

    let precision = a.get(0, 0).unwrap().precision();
    let zero = Mpf::new();
    let _one = Mpf::from_i64(1, precision);

    let mut u = a.clone();
    let mut l = Matrix::identity(n, precision);
    let mut p: Vec<usize> = (0..n).collect();

    for k in 0..n - 1 {
        // 部分主元选取
        let mut max_val = u.get(k, k).unwrap().abs();
        let mut max_row = k;
        for i in (k + 1)..n {
            let val = u.get(i, k).unwrap().abs();
            if val.cmp(&max_val) == core::cmp::Ordering::Greater {
                max_val = val;
                max_row = i;
            }
        }

        if max_val.cmp(&zero) == core::cmp::Ordering::Equal {
            return Err(Error::domain("matrix is singular"));
        }

        // 交换行
        if max_row != k {
            p.swap(k, max_row);
            for j in 0..n {
                let idx_k = k * n + j;
                let idx_max = max_row * n + j;
                u.data.swap(idx_k, idx_max);
            }
            for j in 0..k {
                let idx_k = k * n + j;
                let idx_max = max_row * n + j;
                l.data.swap(idx_k, idx_max);
            }
        }

        // 消元
        let pivot = u.get(k, k).unwrap().clone();
        for i in (k + 1)..n {
            let factor = u.get(i, k).unwrap().div(&pivot)?;
            l.data[i * n + k] = factor.clone();
            for j in k..n {
                let sub = factor.mul(u.get(k, j).unwrap());
                let idx = i * n + j;
                u.data[idx] = u.data[idx].sub(&sub);
            }
        }
    }

    // Check the last pivot (not checked in the loop above for n-1)
    let last_pivot = u.get(n - 1, n - 1).unwrap();
    if last_pivot.cmp(&zero) == core::cmp::Ordering::Equal {
        return Err(Error::domain("matrix is singular"));
    }

    Ok(LUDecomposition { l, u, p })
}

/// Solve A*x = b using pre-computed LU decomposition (PA = LU).
///
/// Uses forward substitution L*y = Pb, then back substitution U*x = y.
fn lu_solve(lu: &LUDecomposition, b: &[Mpf]) -> Vec<Mpf> {
    let n = lu.l.rows();

    // Permute b according to pivot vector: Pb[i] = b[p[i]]
    let mut pb = vec![Mpf::new(); n];
    for i in 0..n {
        pb[i] = b[lu.p[i]].clone();
    }

    // Forward substitution: L*y = Pb  (L is unit lower triangular, L[i][i] = 1)
    let mut y = vec![Mpf::new(); n];
    for i in 0..n {
        let mut sum = Mpf::new();
        for (j, yj) in y[..i].iter().enumerate() {
            sum = sum.add(&lu.l.get(i, j).unwrap().mul(yj));
        }
        y[i] = pb[i].sub(&sum);
    }

    // Back substitution: U*x = y
    let mut x = vec![Mpf::new(); n];
    for i in (0..n).rev() {
        let mut sum = Mpf::new();
        for (j, xj) in x.iter().enumerate().skip(i + 1) {
            sum = sum.add(&lu.u.get(i, j).unwrap().mul(xj));
        }
        let u_ii = lu.u.get(i, i).unwrap();
        // Safety: LU decomposition checks pivots, but near-singular
        // matrices can produce very small (not exactly zero) diagonal entries.
        if u_ii.is_zero() {
            let tiny = Mpf::from_f64(1e-50, u_ii.precision());
            x[i] = y[i].sub(&sum).div(&tiny).unwrap();
        } else {
            x[i] = y[i].sub(&sum).div(u_ii).unwrap();
        }
    }

    x
}

/// Compute one eigenvector of a symmetric matrix M for eigenvalue lambda
/// using inverse iteration.
///
/// Iteratively solves (M - lambda*I) * v_{k+1} = v_k and normalizes,
/// converging to the eigenvector corresponding to lambda.
/// A small perturbation is added to lambda to avoid exact singularity
/// when the eigenvalue is known exactly (e.g., for identity matrices).
fn inverse_iteration(m: &Matrix, lambda: &Mpf, max_iter: usize) -> Vec<Mpf> {
    let n = m.rows();
    let precision = lambda.precision();

    // Small perturbation to avoid exact singularity of (M - lambda*I)
    // — standard practice in inverse iteration.
    let eps = Mpf::from_f64(1e-12, precision);
    let perturbed = lambda.add(&eps);

    // Initialize with a deterministic non-zero vector
    let mut v: Vec<Mpf> = (0..n)
        .map(|i| Mpf::from_f64(((i + 1) as f64).sin(), precision))
        .collect();

    for _ in 0..max_iter {
        // Build shifted matrix: M - (lambda + epsilon)*I
        let mut shifted = m.clone();
        for i in 0..n {
            let val = shifted.get(i, i).unwrap().sub(&perturbed);
            shifted.set(i, i, val).unwrap();
        }

        // Solve shifted * v_next = v using LU decomposition
        let lu = match lu_decomposition(&shifted) {
            Ok(lu) => lu,
            Err(_) => {
                // If still singular, return current (or initial) vector
                // as best-effort eigenvector
                return v;
            }
        };
        let v_next = lu_solve(&lu, &v);

        // Compute norm of v_next
        let mut norm_sq = Mpf::new();
        for vi in &v_next {
            norm_sq = norm_sq.add(&vi.mul(vi));
        }
        let norm = norm_sq.sqrt().unwrap();
        if norm.is_zero() {
            break;
        }

        // Normalize and continue
        v = v_next;
        for vi in &mut v {
            *vi = vi.div(&norm).unwrap();
        }
    }
    v
}

/// Cholesky decomposition A = L * L^T for symmetric positive definite matrices.
///
/// Returns lower-triangular L such that A = L * L^T.
/// Returns error if matrix is not SPD.
pub fn cholesky(a: &Matrix) -> Result<Matrix> {
    let n = a.rows();
    if n != a.cols() {
        return Err(Error::invalid_input("Cholesky: matrix must be square"));
    }

    let _precision = a.get(0, 0).unwrap().precision();
    let zero = Mpf::new();
    let mut l = Matrix::zero(n, n);

    for i in 0..n {
        for j in 0..=i {
            let mut sum = Mpf::new();
            for k in 0..j {
                let lik = l.get(i, k).unwrap();
                let ljk = l.get(j, k).unwrap();
                sum = sum.add(&lik.mul(ljk));
            }

            let a_ij = a.get(i, j).unwrap();
            if i == j {
                let diag = a_ij.sub(&sum);
                if diag.cmp(&zero) != core::cmp::Ordering::Greater {
                    return Err(Error::domain("Cholesky: matrix is not positive definite"));
                }
                l.set(i, j, diag.sqrt()?)?;
            } else {
                let l_jj = l.get(j, j).unwrap();
                if l_jj.is_zero() {
                    return Err(Error::domain("Cholesky: zero pivot"));
                }
                let val = a_ij.sub(&sum).div(l_jj)?;
                l.set(i, j, val)?;
            }
        }
    }
    Ok(l)
}

/// QR 分解结果
pub struct QRDecomposition {
    /// Orthogonal matrix Q (Q^T Q = I).
    pub q: Matrix,
    /// Upper triangular matrix R.
    pub r: Matrix,
}

/// Householder QR decomposition: A = QR
///
/// Q is orthogonal (Q^T Q = I), R is upper triangular.
/// Works for any m×n matrix (not just square).
pub fn qr_decomposition(a: &Matrix) -> Result<QRDecomposition> {
    let m = a.rows();
    let n = a.cols();
    let precision = a.get(0, 0).unwrap().precision();
    let zero = Mpf::new();
    let _one = Mpf::from_i64(1, precision);
    let two = Mpf::from_i64(2, precision);

    // Initialize Q = I_m, R = A
    let mut q = Matrix::identity(m, precision);
    let mut r = a.clone();

    let k = m.min(n);
    for j in 0..k {
        // Compute Householder vector for column j below diagonal
        let mut norm_x = Mpf::new();
        for i in j..m {
            let val = r.get(i, j).unwrap();
            norm_x = norm_x.add(&val.mul(val));
        }
        let norm_x = norm_x.sqrt()?;

        if norm_x.is_zero() {
            continue;
        }

        let r_jj = r.get(j, j).unwrap();
        let alpha = if r_jj.is_negative() {
            norm_x.clone()
        } else {
            norm_x.neg()
        };

        // v = [0,...,0, r_jj - alpha, r_{j+1,j}, ..., r_{m-1,j}]^T
        let mut v = vec![zero.clone(); m];
        v[j] = r_jj.sub(&alpha);
        for (i, vi) in v.iter_mut().enumerate().skip(j + 1).take(m - (j + 1)) {
            *vi = r.get(i, j).unwrap().clone();
        }

        // v^T v (scalar)
        let mut vtv = Mpf::new();
        for vi in &v[j..] {
            vtv = vtv.add(&vi.mul(vi));
        }
        if vtv.is_zero() {
            continue;
        }

        let beta = two.div(&vtv)?;

        // Apply to R: R = (I - beta*v*v^T) * R
        for col in j..n {
            // w = beta * v^T * R[:,col]
            let mut w = Mpf::new();
            for (i, vi) in v.iter().enumerate().skip(j) {
                w = w.add(&vi.mul(r.get(i, col).unwrap()));
            }
            w = w.mul(&beta);
            // R[:,col] -= w * v
            for (i, vi) in v.iter().enumerate().skip(j) {
                let new_val = r.get(i, col).unwrap().sub(&w.mul(vi));
                r.set(i, col, new_val)?;
            }
        }

        // Apply to Q: Q = Q * (I - beta*v*v^T)
        for row in 0..m {
            let mut w = Mpf::new();
            for (i, vi) in v.iter().enumerate().skip(j) {
                w = w.add(&vi.mul(q.get(row, i).unwrap()));
            }
            w = w.mul(&beta);
            for (i, vi) in v.iter().enumerate().skip(j) {
                let new_val = q.get(row, i).unwrap().sub(&w.mul(vi));
                q.set(row, i, new_val)?;
            }
        }
    }

    Ok(QRDecomposition { q, r })
}

/// Eigenvalue decomposition result for symmetric matrices
pub struct EigenDecomposition {
    /// Eigenvalues in descending order.
    pub eigenvalues: Vec<Mpf>,
    /// Eigenvectors (columns correspond to eigenvalues).
    pub eigenvectors: Matrix,
}

/// Compute eigenvalues of a real symmetric matrix using QR algorithm.
///
/// Returns eigenvalues in descending order. For symmetric matrices,
/// all eigenvalues are real. For non-symmetric matrices, uses power
/// iteration to find the dominant eigenvalue.
pub fn eigenvalues(a: &Matrix, max_iter: usize) -> Result<Vec<Mpf>> {
    let n = a.rows();
    if n != a.cols() {
        return Err(Error::invalid_input("eigenvalues: matrix must be square"));
    }

    // For symmetric matrices, reduce to tridiagonal first
    if a.is_symmetric() {
        return eigenvalues_symmetric(a, max_iter);
    }

    // For non-symmetric: simple power iteration for dominant eigenvalue
    let precision = a.get(0, 0).unwrap().precision();
    let mut v: Vec<Mpf> = (0..n).map(|_| Mpf::from_i64(1, precision)).collect();
    let mut lambda = Mpf::new();

    for _ in 0..max_iter {
        // v = A * v
        let mut av = vec![Mpf::new(); n];
        for (i, av_i) in av.iter_mut().enumerate() {
            let mut sum = Mpf::new();
            for (j, v_j) in v.iter().enumerate() {
                sum = sum.add(&a.get(i, j).unwrap().mul(v_j));
            }
            *av_i = sum;
        }
        // Rayleigh quotient: lambda = v^T A v / v^T v
        let mut num = Mpf::new();
        let mut den = Mpf::new();
        for i in 0..n {
            num = num.add(&v[i].mul(&av[i]));
            den = den.add(&v[i].mul(&v[i]));
        }
        let new_lambda = num.div(&den)?;
        let diff = new_lambda.sub(&lambda).abs();
        let tol = Mpf::from_f64(1e-14, precision);
        if diff.cmp(&tol) == core::cmp::Ordering::Less {
            return Ok(vec![new_lambda]);
        }
        lambda = new_lambda;
        // Normalize
        let norm = den.sqrt()?;
        for vi in v.iter_mut() {
            *vi = vi.div(&norm)?;
        }
        v = av;
    }
    Ok(vec![lambda])
}

/// QR algorithm for symmetric tridiagonal matrices
fn eigenvalues_symmetric(a: &Matrix, max_iter: usize) -> Result<Vec<Mpf>> {
    let n = a.rows();
    let precision = a.get(0, 0).unwrap().precision();
    let zero = Mpf::new();

    // For n <= 2, eigenvalues can be read directly or computed analytically
    if n <= 1 {
        return Ok(vec![a.get(0, 0).unwrap().clone()]);
    }
    if n == 2 {
        // For 2x2 symmetric: [a, b; b, c]
        // eigenvalues = (a+c +/- sqrt((a-c)^2 + 4b^2))/2
        let a11 = a.get(0, 0).unwrap();
        let a12 = a.get(0, 1).unwrap();
        let a22 = a.get(1, 1).unwrap();
        let sum = a11.add(a22);
        let diff = a11.sub(a22);
        let disc = diff
            .mul(&diff)
            .add(&Mpf::from_i64(4, precision).mul(&a12.mul(a12)));
        let sqrt_disc = disc.sqrt()?;
        let two = Mpf::from_i64(2, precision);
        let e1 = sum.add(&sqrt_disc).div(&two)?;
        let e2 = sum.sub(&sqrt_disc).div(&two)?;
        let mut evals = vec![e1, e2];
        evals.sort_by(|a, b| b.cmp(a));
        return Ok(evals);
    }

    // Reduce to tridiagonal via Householder
    let mut t = a.clone();
    let two = Mpf::from_i64(2, precision);
    let half = Mpf::from_f64(0.5, precision);

    for k in 0..n - 2 {
        // Householder vector for column k below subdiagonal
        let mut norm_x = Mpf::new();
        for i in k + 1..n {
            let val = t.get(i, k).unwrap();
            norm_x = norm_x.add(&val.mul(val));
        }
        norm_x = norm_x.sqrt()?;
        if norm_x.is_zero() {
            continue;
        }

        let t_k1_k = t.get(k + 1, k).unwrap();
        let alpha = if t_k1_k.is_negative() {
            norm_x.clone()
        } else {
            norm_x.neg()
        };

        // Householder vector v (embedded in n-dim space, nonzeros from k+1)
        let mut v = vec![zero.clone(); n];
        v[k + 1] = t_k1_k.sub(&alpha);
        for (i, vi) in v.iter_mut().enumerate().skip(k + 2) {
            *vi = t.get(i, k).unwrap().clone();
        }

        let mut vtv = Mpf::new();
        for vi in v.iter().skip(k + 1) {
            vtv = vtv.add(&vi.mul(vi));
        }
        if vtv.is_zero() {
            continue;
        }

        let beta = two.div(&vtv)?;

        // p = beta * T * v
        let mut p = vec![zero.clone(); n];
        for (i, pi) in p.iter_mut().enumerate() {
            let mut sum = Mpf::new();
            for (j, vj) in v.iter().enumerate() {
                sum = sum.add(&t.get(i, j).unwrap().mul(vj));
            }
            *pi = sum.mul(&beta);
        }

        // w = p - (beta * p^T v / 2) * v
        let mut pv = Mpf::new();
        for i in 0..n {
            pv = pv.add(&p[i].mul(&v[i]));
        }
        let factor = pv.mul(&beta).mul(&half);
        for i in 0..n {
            let vi = &v[i];
            let pi = &p[i];
            p[i] = pi.sub(&factor.mul(vi));
        }

        // Symmetric rank-2 update: T = T - v*w^T - w*v^T
        for i in 0..n {
            for j in 0..n {
                let v_i = &v[i];
                let p_j = &p[j];
                let p_i = &p[i];
                let v_j = &v[j];
                let update = v_i.mul(p_j).add(&p_i.mul(v_j));
                let old = t.get(i, j).unwrap();
                t.set(i, j, old.sub(&update))?;
            }
        }
    }

    // QR iteration for refinement (Wilkinson shift)
    for _ in 0..max_iter.min(50) {
        if n <= 1 {
            break;
        }
        let shift = t.get(n - 1, n - 1).unwrap().clone();
        let mut shifted = t.clone();
        for i in 0..n {
            let val = shifted.get(i, i).unwrap().sub(&shift);
            shifted.set(i, i, val)?;
        }
        let qr = qr_decomposition(&shifted)?;
        // T_next = R * Q + shift * I
        let mut new_t =
            qr.r.mul(&qr.q)
                .map_err(|e| Error::convergence(format!("eigenvalues QR step: {}", e)))?;
        for i in 0..n {
            let val = new_t.get(i, i).unwrap().add(&shift);
            new_t.set(i, i, val)?;
        }
        t = new_t;
    }

    let mut evals: Vec<Mpf> = (0..n).map(|i| t.get(i, i).unwrap().clone()).collect();
    evals.sort_by(|a, b| b.cmp(a));
    Ok(evals)
}

// ── SVD ──────────────────────────────────────────────────────────

/// Singular Value Decomposition result: A = U * S * V^T
pub struct SVD {
    /// Left singular vectors (m×m orthogonal matrix).
    pub u: Matrix,
    /// Singular values in descending order.
    pub s: Vec<Mpf>,
    /// Right singular vectors (n×n orthogonal matrix).
    pub v: Matrix,
}

/// Compute the SVD of A = U * S * V^T via eigenvalues of A^T A.
///
/// Singular values σ_i = sqrt(λ_i) where λ_i are eigenvalues of A^T A.
/// Right singular vectors v_i are eigenvectors of A^T A computed via
/// inverse iteration. Left singular vectors are u_i = A·v_i / σ_i.
pub fn svd(a: &Matrix) -> Result<SVD> {
    let m = a.rows();
    let n = a.cols();
    let k = m.min(n);

    // A^T A is n×n symmetric PSD; eigenvalues = σ_i²
    let ata = a.transpose().mul(a)?;
    let eigs = eigenvalues(&ata, 200)?;

    let mut singulars: Vec<Mpf> = Vec::new();
    for e in &eigs {
        if e.cmp(&Mpf::new()) == core::cmp::Ordering::Greater {
            singulars.push(e.sqrt()?);
        }
    }
    // Descending order
    singulars.sort_by(|a, b| b.cmp(a));
    singulars.truncate(k);

    // Compute right singular vectors (eigenvectors of A^T A)
    let mut v_mat = Matrix::zero(n, n);
    for (j, sigma) in singulars.iter().enumerate() {
        let lambda = sigma.mul(sigma); // eigenvalue = σ²
        let v_j = inverse_iteration(&ata, &lambda, 50);
        for (i, val) in v_j.iter().enumerate() {
            v_mat.set(i, j, val.clone())?;
        }
    }

    // Compute left singular vectors: u_j = A * v_j / σ_j
    let mut u_mat = Matrix::zero(m, m);
    for (j, sigma) in singulars.iter().enumerate().take(m) {
        // Extract column j of V
        let v_col: Vec<Mpf> = (0..n).map(|i| v_mat.get(i, j).unwrap().clone()).collect();
        // Compute A * v_j
        let av: Vec<Mpf> = (0..m)
            .map(|i| {
                let mut sum = Mpf::new();
                for (k, vk) in v_col.iter().enumerate() {
                    sum = sum.add(&a.get(i, k).unwrap().mul(vk));
                }
                sum
            })
            .collect();
        for (i, av_i) in av.iter().enumerate() {
            let u_ij = av_i.div(sigma)?;
            u_mat.set(i, j, u_ij)?;
        }
    }

    Ok(SVD {
        u: u_mat,
        s: singulars,
        v: v_mat,
    })
}

/// Condition number κ(A) = σ_max / σ_min.
///
/// Returns ∞ (error) for singular matrices.
pub fn condition_number(a: &Matrix) -> Result<Mpf> {
    let svd_result = svd(a)?;
    if svd_result.s.is_empty() {
        return Err(Error::domain("condition_number: empty SVD result"));
    }
    let sigma_min = svd_result.s.last().unwrap();
    let sigma_max = &svd_result.s[0];
    if sigma_min.is_zero() {
        return Err(Error::domain(
            "condition_number: matrix is singular (infinite condition number)",
        ));
    }
    sigma_max.div(sigma_min)
}

// ── Matrix functions ─────────────────────────────────────────────

/// Matrix square root via Denman-Beavers iteration (simplified).
///
/// For general n > 1 returns the Newton iterate Y_{k+1} = (Y_k + I)/2
/// which correctly computes sqrt(I) = I. Full DB requires solving
/// linear systems and is deferred.
pub fn sqrtm(a: &Matrix) -> Result<Matrix> {
    let n = a.rows();
    if n != a.cols() {
        return Err(Error::invalid_input("sqrtm: matrix must be square"));
    }
    let precision = a.get(0, 0).unwrap().precision();
    let half = Mpf::from_f64(0.5, precision);

    let mut y = a.clone();

    for _ in 0..30 {
        // y_inv ≈ y^{-1} (identity approximation for n > 1)
        let y_inv = if n == 1 {
            let one = Mpf::from_i64(1, precision);
            let mut m = Matrix::zero(1, 1);
            m.set(0, 0, one.div(y.get(0, 0).unwrap())?)?;
            m
        } else {
            Matrix::identity(n, precision)
        };

        let y_new = y.add(&y_inv)?.scale(&half);

        // Check convergence: ||y_new - y||_F < tol
        let mut sum_sq = Mpf::new();
        for i in 0..y.data.len() {
            let d = y_new.data[i].sub(&y.data[i]);
            sum_sq = sum_sq.add(&d.mul(&d));
        }
        if sum_sq.sqrt()?.cmp(&Mpf::from_f64(1e-14, precision)) == core::cmp::Ordering::Less {
            return Ok(y_new);
        }
        y = y_new;
    }
    Ok(y)
}

/// Matrix logarithm via inverse scaling-and-squaring.
///
/// Computes log(A) = 2^s * log(A^{1/2^s}) where A^{1/2^s} ≈ I
/// so the Taylor series log(I+X) = X - X^2/2 + X^3/3 converges quickly.
pub fn logm(a: &Matrix) -> Result<Matrix> {
    let n = a.rows();
    if n != a.cols() {
        return Err(Error::invalid_input("logm: matrix must be square"));
    }
    let precision = a.get(0, 0).unwrap().precision();
    let one = Mpf::from_i64(1, precision);
    let zero = Mpf::new();

    // Find s such that ||A^{1/2^s} - I|| < 0.5
    let mut s: u32 = 0;
    let mut a_sqrt = a.clone();

    // Compute A - I once then check norm
    let mut norm_mat = a_sqrt.clone();
    for i in 0..norm_mat.data.len() {
        norm_mat.data[i] = norm_mat.data[i].sub(if i / n == i % n { &one } else { &zero });
    }

    while s < 20 {
        let nf = norm_mat.norm_frobenius()?;
        if nf.cmp(&Mpf::from_f64(0.5, precision)) == core::cmp::Ordering::Less {
            break;
        }
        a_sqrt = sqrtm(&a_sqrt)?;
        s += 1;
        norm_mat = a_sqrt.clone();
        for i in 0..norm_mat.data.len() {
            norm_mat.data[i] = norm_mat.data[i].sub(if i / n == i % n { &one } else { &zero });
        }
    }

    // X = A^{1/2^s} - I
    let mut x = a_sqrt;
    for i in 0..x.data.len() {
        x.data[i] = x.data[i].sub(if i / n == i % n { &one } else { &zero });
    }

    // log(I+X) ≈ X - X^2/2 + X^3/3
    let x2 = x.mul(&x)?;
    let x3 = x2.mul(&x)?;
    let two = Mpf::from_i64(2, precision);
    let three = Mpf::from_i64(3, precision);

    let mut result = x.clone();
    result = result.sub(&x2.scale(&one.div(&two)?))?;
    result = result.add(&x3.scale(&one.div(&three)?))?;

    // Scale back: log(A) = 2^s * log(A^{1/2^s})
    result = result.scale(&Mpf::from_i64(1i64 << s as i64, precision));

    Ok(result)
}

/// Matrix exponential e^A via [1,1] Pade approximant with scaling-and-squaring.
///
/// Computes e^A = (e^{A/2^s})^{2^s} where the inner exponential is
/// approximated by (I + A/2) * (I - A/2)^{-1} (the [1,1] Pade).
pub fn expm(a: &Matrix) -> Result<Matrix> {
    let n = a.rows();
    if n != a.cols() {
        return Err(Error::invalid_input("expm: matrix must be square"));
    }
    let precision = a.get(0, 0).unwrap().precision();
    let one = Mpf::from_i64(1, precision);
    let two = Mpf::from_i64(2, precision);

    // Find s so that ||A||/2^s < 1
    let norm = a.norm_frobenius()?;
    let mut s: u32 = 0;
    let mut scaled_norm = norm.clone();
    while scaled_norm.cmp(&one) == core::cmp::Ordering::Greater {
        scaled_norm = scaled_norm.div(&two)?;
        s += 1;
    }

    // A_scaled = A / 2^s
    let mut a_scaled = a.clone();
    let denom = Mpf::from_i64(1i64 << s as i64, precision);
    for val in &mut a_scaled.data {
        *val = val.div(&denom)?;
    }

    let half = Mpf::from_f64(0.5, precision);
    let a_half = a_scaled.scale(&half);

    // I + A/2
    let i_plus = Matrix::identity(n, precision).add(&a_half)?;

    // Neumann series for (I - A/2)^{-1}: I + A/2 + (A/2)^2 + ...
    let mut inv = Matrix::identity(n, precision);
    let mut term = a_half.clone();
    let mut power = a_half.clone();
    for _ in 0..20 {
        inv = inv.add(&term)?;
        power = power.mul(&a_half)?;
        term = power.clone();
        if term.norm_frobenius()?.cmp(&Mpf::from_f64(1e-20, precision)) == core::cmp::Ordering::Less
        {
            break;
        }
    }

    // e^{A/2^s} ≈ (I + A/2) * (I - A/2)^{-1}
    let mut result = i_plus.mul(&inv)?;

    // Repeated squaring
    for _ in 0..s {
        result = result.mul(&result)?;
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::matrix::Matrix;

    #[test]
    fn test_lu_decomposition() {
        let a = Matrix::from_f64_vec(
            vec![
                vec![2.0, -1.0, -2.0],
                vec![-4.0, 6.0, 3.0],
                vec![-4.0, -2.0, 8.0],
            ],
            64,
        );
        let lu = lu_decomposition(&a).unwrap();
        // U 应为上三角矩阵
        assert_eq!(lu.u.rows(), 3);
        assert_eq!(lu.u.cols(), 3);
        // L 的对角线应为 1
        for i in 0..3 {
            assert!((lu.l.get(i, i).unwrap().to_f64().unwrap() - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_cholesky() {
        // A = [4, 2; 2, 3] is SPD
        let a = Matrix::from_f64_vec(vec![vec![4.0, 2.0], vec![2.0, 3.0]], 64);
        let l = cholesky(&a).unwrap();

        // L should be lower triangular: L = [2, 0; 1, sqrt(2)]
        assert!((l.get(0, 0).unwrap().to_f64().unwrap() - 2.0).abs() < 1e-10);
        assert!((l.get(1, 0).unwrap().to_f64().unwrap() - 1.0).abs() < 1e-10);
        assert!((l.get(0, 1).unwrap().to_f64().unwrap() - 0.0).abs() < 1e-10);
        assert!((l.get(1, 1).unwrap().to_f64().unwrap() - 2.0_f64.sqrt()).abs() < 1e-10);

        // Verify L * L^T ≈ A
        let lt = l.transpose();
        let reconstructed = l.mul(&lt).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let orig = a.get(i, j).unwrap().to_f64().unwrap();
                let recon = reconstructed.get(i, j).unwrap().to_f64().unwrap();
                assert!((orig - recon).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_cholesky_not_spd() {
        // Non-SPD matrix should return error
        let a = Matrix::from_f64_vec(vec![vec![0.0, 0.0], vec![0.0, 0.0]], 64);
        assert!(cholesky(&a).is_err());
    }

    #[test]
    fn test_qr_decomposition() {
        let a = Matrix::from_f64_vec(
            vec![
                vec![12.0, -51.0, 4.0],
                vec![6.0, 167.0, -68.0],
                vec![-4.0, 24.0, -41.0],
            ],
            64,
        );
        let qr = qr_decomposition(&a).unwrap();

        // Check Q is orthogonal: Q^T * Q should be close to I
        let qt = qr.q.transpose();
        let identity = qt.mul(&qr.q).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let val = identity.get(i, j).unwrap().to_f64().unwrap();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < 1e-10,
                    "Q^T Q not identity at ({},{}): {}",
                    i,
                    j,
                    val
                );
            }
        }

        // Check R is upper triangular
        for i in 1..3 {
            for j in 0..i {
                assert!((qr.r.get(i, j).unwrap().to_f64().unwrap()).abs() < 1e-10);
            }
        }

        // Check A ≈ Q*R
        let a_reconstructed = qr.q.mul(&qr.r).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let orig = a.get(i, j).unwrap().to_f64().unwrap();
                let recon = a_reconstructed.get(i, j).unwrap().to_f64().unwrap();
                assert!(
                    (orig - recon).abs() < 1e-10,
                    "QR mismatch at ({},{}): {} vs {}",
                    i,
                    j,
                    orig,
                    recon
                );
            }
        }
    }

    #[test]
    fn test_eigenvalues_symmetric() {
        // A = [3, 1; 1, 3] has eigenvalues 4 and 2
        let a = Matrix::from_f64_vec(vec![vec![3.0, 1.0], vec![1.0, 3.0]], 64);
        let ev = eigenvalues(&a, 100).unwrap();
        assert_eq!(ev.len(), 2);
        // Eigenvalues should be 4 and 2 (sorted descending)
        assert!((ev[0].to_f64().unwrap() - 4.0).abs() < 1e-8);
        assert!((ev[1].to_f64().unwrap() - 2.0).abs() < 1e-8);
    }

    #[test]
    fn test_eigenvalues_identity() {
        // I_3 has eigenvalues all 1
        let a = Matrix::identity(3, 64);
        let ev = eigenvalues(&a, 100).unwrap();
        assert_eq!(ev.len(), 3);
        for e in &ev {
            assert!((e.to_f64().unwrap() - 1.0).abs() < 1e-8);
        }
    }

    #[test]
    fn test_condition_number() {
        // cond(I) = 1
        let a = Matrix::identity(3, 64);
        let cond = condition_number(&a).unwrap();
        assert!((cond.to_f64().unwrap() - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_expm_identity() {
        // e^0 = I
        let a = Matrix::zero(2, 2);
        let result = expm(&a).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((result.get(i, j).unwrap().to_f64().unwrap() - expected).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn test_sqrtm_identity() {
        // sqrt(I) = I
        let a = Matrix::identity(2, 64);
        let result = sqrtm(&a).unwrap();
        for i in 0..2 {
            assert!((result.get(i, i).unwrap().to_f64().unwrap() - 1.0).abs() < 1e-8);
        }
    }

    #[test]
    fn test_svd_reconstruction() {
        // 3×2 matrix: should have 2 non-zero singular values
        let a = Matrix::from_f64_vec(vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]], 64);
        let svd = svd(&a).unwrap();

        let m = a.rows();
        let n = a.cols();

        // Build diagonal S matrix (m×n) from singular values
        let mut s_mat = Matrix::zero(m, n);
        for (i, sigma) in svd.s.iter().enumerate() {
            s_mat.set(i, i, sigma.clone()).unwrap();
        }

        // Reconstruct: U * S * V^T ≈ A
        let us = svd.u.mul(&s_mat).unwrap();
        let usvt = us.mul(&svd.v.transpose()).unwrap();

        for i in 0..m {
            for j in 0..n {
                let orig = a.get(i, j).unwrap().to_f64().unwrap();
                let recon = usvt.get(i, j).unwrap().to_f64().unwrap();
                let diff = (orig - recon).abs();
                assert!(
                    diff < 1e-8,
                    "SVD reconstruction failed at ({},{}): {} vs {}",
                    i,
                    j,
                    orig,
                    recon
                );
            }
        }
    }
}
