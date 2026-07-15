use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mynum::Mpz;

fn bench_mynum_vs_numbigint_multiplication(c: &mut Criterion) {
    let mut group = c.benchmark_group("Multiplication Comparison");

    let sizes = [64, 256, 1024, 4096];

    for &size in &sizes {
        let a_mynum = Mpz::random_bits(size).unwrap();
        let b_mynum = Mpz::random_bits(size).unwrap();

        let a_str = a_mynum.to_string(10);
        let b_str = b_mynum.to_string(10);
        let a_nb = num_bigint::BigInt::parse_bytes(a_str.as_bytes(), 10).unwrap();
        let b_nb = num_bigint::BigInt::parse_bytes(b_str.as_bytes(), 10).unwrap();

        group.bench_function(&format!("mynum_{}_bits", size), |bencher| {
            bencher.iter(|| black_box(a_mynum.mul(&b_mynum)));
        });

        group.bench_function(&format!("num_bigint_{}_bits", size), |bencher| {
            bencher.iter(|| black_box(&a_nb * &b_nb));
        });
    }

    group.finish();
}

fn bench_mynum_vs_numbigint_addition(c: &mut Criterion) {
    let mut group = c.benchmark_group("Addition Comparison");

    let sizes = [64, 256, 1024, 4096];

    for &size in &sizes {
        let a = Mpz::random_bits(size).unwrap();
        let b = Mpz::random_bits(size).unwrap();

        let a_str = a.to_string(10);
        let b_str = b.to_string(10);
        let a_nb = num_bigint::BigInt::parse_bytes(a_str.as_bytes(), 10).unwrap();
        let b_nb = num_bigint::BigInt::parse_bytes(b_str.as_bytes(), 10).unwrap();

        group.bench_function(&format!("mynum_add_{}_bits", size), |bencher| {
            bencher.iter(|| black_box(a.add(&b)));
        });

        group.bench_function(&format!("num_bigint_add_{}_bits", size), |bencher| {
            bencher.iter(|| black_box(&a_nb + &b_nb));
        });
    }

    group.finish();
}

criterion_group!(
    comparative_benches,
    bench_mynum_vs_numbigint_multiplication,
    bench_mynum_vs_numbigint_addition
);
criterion_main!(comparative_benches);
