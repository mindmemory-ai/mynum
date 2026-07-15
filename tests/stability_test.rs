//! 长时间运行稳定性测试
//!
//! 运行方式:
//!   快速烟雾测试: cargo test --test stability_test -- --ignored   (Ctrl+C after 30s)
//!   完整 72h:     cargo test --test stability_test -- --ignored --test-threads 1

use mynum::{Mpf, Mpz};

#[test]
#[ignore]
fn stability_multiplication_loop() {
    let iterations = 10_000_000;
    let mut a = Mpz::from_i64(2);
    let b = Mpz::from_i64(3);
    let two = Mpz::from_i64(2);

    for i in 0..iterations {
        a = a.mul(&b);
        if i % 100 == 0 {
            a = a.div(&two).unwrap();
        }
        if i % 100000 == 0 {
            eprintln!("mpz iteration {}/{}", i, iterations);
        }
    }
    assert!(!a.is_zero());
    eprintln!("mpz stability test complete: {}", a.to_string(10));
}

#[test]
#[ignore]
fn stability_float_operations() {
    let iterations = 5_000_000;
    let mut x = Mpf::from_f64(1.0, 128);
    let two = Mpf::from_i64(2, 128);

    for i in 0..iterations {
        x = x.mul(&two);
        if x.to_f64().unwrap_or(0.0) > 1e100 {
            x = Mpf::from_f64(1.0, 128);
        }
        if i % 50000 == 0 {
            eprintln!("float iteration {}/{}", i, iterations);
        }
    }
    eprintln!("mpf stability test complete: {}", x.to_string(10));
}

#[test]
#[ignore]
fn stability_matrix_operations() {
    use mynum::linalg::matrix::Matrix;

    let iterations = 100_000;
    let a = Matrix::from_f64_vec(vec![vec![1.0, 2.0], vec![3.0, 4.0]], 64);

    for i in 0..iterations {
        let _b = a.add(&a).unwrap();
        let _c = a.mul(&a).unwrap();
        if i % 10000 == 0 {
            eprintln!("matrix iteration {}/{}", i, iterations);
        }
    }
    eprintln!("matrix stability test complete");
}
