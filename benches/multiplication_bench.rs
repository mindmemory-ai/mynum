use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mynum::{Mpz, MultiplicationBackend};

fn bench_schoolbook_multiplication(c: &mut Criterion) {
    let mut group = c.benchmark_group("Schoolbook Multiplication");

    // 测试不同大小的数
    let sizes = [32, 64, 128, 256];

    for &size in &sizes {
        let a = Mpz::random_bits(size).unwrap();
        let b = Mpz::random_bits(size).unwrap();

        group.bench_function(&format!("{}_bits", size), |bencher| {
            bencher.iter(|| black_box(a.mul_schoolbook(&b)));
        });
    }

    group.finish();
}

fn bench_karatsuba_multiplication(c: &mut Criterion) {
    let mut group = c.benchmark_group("Karatsuba Multiplication");

    // 测试不同大小的数
    let sizes = [64, 128, 256, 512];

    for &size in &sizes {
        let a = Mpz::random_bits(size).unwrap();
        let b = Mpz::random_bits(size).unwrap();

        group.bench_function(&format!("{}_bits", size), |bencher| {
            bencher.iter(|| black_box(a.mul_karatsuba(&b)));
        });
    }

    group.finish();
}

fn bench_adaptive_multiplication(c: &mut Criterion) {
    let mut group = c.benchmark_group("Adaptive Multiplication");

    // 测试不同大小的数
    let sizes = [32, 64, 128, 256, 512, 1024];

    for &size in &sizes {
        let a = Mpz::random_bits(size).unwrap();
        let b = Mpz::random_bits(size).unwrap();

        group.bench_function(&format!("{}_bits", size), |bencher| {
            bencher.iter(|| black_box(a.mul_with_backend(&b, MultiplicationBackend::Adaptive)));
        });
    }

    group.finish();
}

fn bench_multiplication_backends(c: &mut Criterion) {
    let mut group = c.benchmark_group("Multiplication Backends Comparison");

    let a = Mpz::random_bits(256).unwrap();
    let b = Mpz::random_bits(256).unwrap();

    group.bench_function("Schoolbook", |bencher| {
        bencher.iter(|| black_box(a.mul_with_backend(&b, MultiplicationBackend::Schoolbook)));
    });

    group.bench_function("Karatsuba", |bencher| {
        bencher.iter(|| black_box(a.mul_with_backend(&b, MultiplicationBackend::Karatsuba)));
    });

    group.bench_function("Adaptive", |bencher| {
        bencher.iter(|| black_box(a.mul_with_backend(&b, MultiplicationBackend::Adaptive)));
    });

    group.finish();
}

fn bench_random_number_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("Random Number Generation");

    group.bench_function("64_bits", |bencher| {
        bencher.iter(|| black_box(Mpz::random_bits(64).unwrap()));
    });

    group.bench_function("128_bits", |bencher| {
        bencher.iter(|| black_box(Mpz::random_bits(128).unwrap()));
    });

    group.bench_function("256_bits", |bencher| {
        bencher.iter(|| black_box(Mpz::random_bits(256).unwrap()));
    });

    group.finish();
}

fn bench_threshold_calibration(c: &mut Criterion) {
    let mut group = c.benchmark_group("Threshold Calibration");

    // Test each backend at sizes spanning the default thresholds
    let test_sizes: &[usize] = &[16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192];
    let backends: &[(MultiplicationBackend, &str)] = &[
        (MultiplicationBackend::Schoolbook, "schoolbook"),
        (MultiplicationBackend::Karatsuba, "karatsuba"),
        (MultiplicationBackend::ToomCook3, "toom3"),
        (MultiplicationBackend::ToomCook4, "toom4"),
        (MultiplicationBackend::FFT, "fft"),
    ];

    for &size in test_sizes {
        let a = Mpz::random_bits(size).unwrap();
        let b = Mpz::random_bits(size).unwrap();

        for &(backend, name) in backends {
            group.bench_function(&format!("{}_{}_bits", name, size), |bencher| {
                bencher.iter(|| black_box(a.mul_with_backend(&b, backend)));
            });
        }
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_schoolbook_multiplication,
    bench_karatsuba_multiplication,
    bench_adaptive_multiplication,
    bench_multiplication_backends,
    bench_random_number_generation,
    bench_threshold_calibration
);
criterion_main!(benches);
