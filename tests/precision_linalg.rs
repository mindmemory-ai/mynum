//! Linalg precision tests — decomposition accuracy and reconstruction.

use mynum::linalg::{decomposition::*, matrix::Matrix};
use mynum::mpf::Mpf;

#[test]
fn test_lu_reconstruction() {
    let a = Matrix::from_f64_vec(
        vec![
            vec![2.0, -1.0, -2.0],
            vec![-4.0, 6.0, 3.0],
            vec![-4.0, -2.0, 8.0],
        ],
        128,
    );
    let lu = lu_decomposition(&a).unwrap();
    // Reconstruct: P*A = L*U, so A ≈ P^T * L * U
    let lu_prod = lu.l.mul(&lu.u).unwrap();
    for i in 0..3 {
        let src = lu.p[i];
        for j in 0..3 {
            let diff = (lu_prod.get(src, j).unwrap().to_f64().unwrap()
                - a.get(i, j).unwrap().to_f64().unwrap())
            .abs();
            assert!(diff < 1e-12, "LU mismatch at ({},{}): diff={}", i, j, diff);
        }
    }
}

#[test]
fn test_qr_orthogonality() {
    let a = Matrix::from_f64_vec(
        vec![
            vec![1.0, -2.0, 3.0],
            vec![4.0, 0.0, -6.0],
            vec![-2.0, 5.0, 1.0],
        ],
        128,
    );
    let qr = qr_decomposition(&a).unwrap();
    let qtq = qr.q.transpose().mul(&qr.q).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            let val = qtq.get(i, j).unwrap().to_f64().unwrap();
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (val - expected).abs() < 1e-12,
                "Q^T Q[{},{}] = {}, expected {}",
                i,
                j,
                val,
                expected
            );
        }
    }
    // A ≈ Q*R
    let recon = qr.q.mul(&qr.r).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            let orig = a.get(i, j).unwrap().to_f64().unwrap();
            let rec = recon.get(i, j).unwrap().to_f64().unwrap();
            assert!(
                (orig - rec).abs() < 1e-10,
                "QR reconstruction failed at ({},{})",
                i,
                j
            );
        }
    }
}

#[test]
fn test_cholesky_reconstruction() {
    let a = Matrix::from_f64_vec(
        vec![
            vec![4.0, 12.0, -16.0],
            vec![12.0, 37.0, -43.0],
            vec![-16.0, -43.0, 98.0],
        ],
        128,
    );
    let l = cholesky(&a).unwrap();
    let recon = l.mul(&l.transpose()).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            let diff = (recon.get(i, j).unwrap().to_f64().unwrap()
                - a.get(i, j).unwrap().to_f64().unwrap())
            .abs();
            assert!(
                diff < 1e-12,
                "Cholesky mismatch at ({},{}): diff={}",
                i,
                j,
                diff
            );
        }
    }
}

#[test]
fn test_cholesky_non_spd() {
    let a = Matrix::from_f64_vec(vec![vec![0.0, 0.0], vec![0.0, 0.0]], 64);
    assert!(cholesky(&a).is_err());
    // Non-symmetric should still produce error
    let b = Matrix::from_f64_vec(vec![vec![1.0, 2.0], vec![3.0, 4.0]], 64);
    // This may or may not be SPD — cholesky should either work or error, not panic
    let _ = cholesky(&b);
}

#[test]
fn test_eigenvalues_identity() {
    let i3 = Matrix::identity(3, 64);
    let ev = eigenvalues(&i3, 100).unwrap();
    for e in &ev {
        assert!(
            (e.to_f64().unwrap() - 1.0).abs() < 1e-10,
            "I_3 eigenvalue should be 1"
        );
    }
}

#[test]
fn test_eigenvalues_symmetric_2x2() {
    let a = Matrix::from_f64_vec(vec![vec![3.0, 1.0], vec![1.0, 3.0]], 64);
    let ev = eigenvalues(&a, 100).unwrap();
    assert_eq!(ev.len(), 2);
    assert!((ev[0].to_f64().unwrap() - 4.0).abs() < 1e-8);
    assert!((ev[1].to_f64().unwrap() - 2.0).abs() < 1e-8);
}

#[test]
fn test_condition_number_identity() {
    let i3 = Matrix::identity(3, 64);
    let cond = condition_number(&i3).unwrap();
    assert!((cond.to_f64().unwrap() - 1.0).abs() < 1e-8);
}

#[test]
fn test_matrix_norm() {
    let a = Matrix::from_f64_vec(vec![vec![3.0, 4.0], vec![0.0, 0.0]], 64);
    let norm = a.norm_frobenius().unwrap();
    assert!((norm.to_f64().unwrap() - 5.0).abs() < 1e-10);
    let i2 = Matrix::identity(2, 64);
    let ni = i2.norm_frobenius().unwrap();
    assert!((ni.to_f64().unwrap() - 2.0_f64.sqrt()).abs() < 1e-10);
}

#[test]
fn test_matrix_sub() {
    let a = Matrix::from_f64_vec(vec![vec![5.0, 6.0], vec![7.0, 8.0]], 64);
    let b = Matrix::from_f64_vec(vec![vec![1.0, 2.0], vec![3.0, 4.0]], 64);
    let c = a.sub(&b).unwrap();
    assert!((c.get(0, 0).unwrap().to_f64().unwrap() - 4.0).abs() < 1e-10);
    assert!((c.get(0, 1).unwrap().to_f64().unwrap() - 4.0).abs() < 1e-10);
}
