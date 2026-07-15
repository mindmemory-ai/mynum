# Changelog

All notable changes to mynum are documented in this file.

## [0.5.3] — 2026-06-21

### Added
- `ErrorSeverity` enum (Warning/Error/Fatal) for classifying error impact
- `RecoveryStrategy` struct with `retry()`, `retry_with_precision()`, and
  `to_adaptive()` methods for systematic error recovery
- `benchmark_thresholds()` calibration using `std::time::Instant` to measure
  actual crossover points between multiplication backends
- `examples/threshold_calibration.rs` with `--apply` flag
- `#[cfg_attr(tarpaulin, ignore)]` on 5 large tests to prevent timeout under
  tarpaulin instrumentation

### Changed
- CI coverage threshold: `cargo tarpaulin --fail-under 70` → `80`
- PLAN.md archived to `docs/archive/PLAN-2026-06-21-completed.md` (all
  actionable items 100% complete)

## [0.5.2] — 2026-06-20

### Added
- 18 new special functions: `expint_ei`, `erfi`, `sinint_si`, `cosint_ci`,
  `fresnel_s`, `fresnel_c`, `struve_h`, `lommel_s`, `whittaker_m`,
  `sinhint_shi`, `coshint_chi`, `bessel_j_real`, `bessel_y_real`,
  `bessel_i_real`, `bessel_k_real`, `elliptic_pi`, `jacobi_sn`, `jacobi_cn`,
  `jacobi_dn`
- 48 precision_special tests

### Fixed
- `gamma()`: `is_integer()` → `is_exact_integer()` f64 integer truncation
- `expint_ei`: optimal truncation for asymptotic expansion
- `fresnel_aux_fg`: g coefficient corrected per DLMF 7.12.2

## [0.5.1] — 2026-06-20

### Added
- `Mpf::from_str` arbitrary-precision decimal parsing (>15 digits)
- `ParallelExecutor` unified abstraction over rayon + std::thread
- Meijer G for |z| ≥ 1 via Slater inversion formula
- Full SVD U/V matrices via inverse iteration eigenvectors
- `spherical_harmonic_real`, `lerchphi`
- `AtomicConfig::on_change()` callback registration

### Changed
- `pow_big` → `pow(&Mpz)`, `pow(u32)` → `pow_u32`, added `Mpf::pow_big`

## [0.5.0] — 2026-06-20

### Added
- Bessel/Airy: 15 functions (complete J/Y/I/K/Hankel/half-integer families)
- Numerical: `quad` (adaptive Gauss-Kronrod), `nsum` (Wynn ε), `limit` (Richardson)
- Gamma+: `loggamma`, `digamma`, `polygamma`, `polylog`, `hurwitz_zeta`, `lambert_w`
- Linalg: QR, Cholesky, eigenvalues, SVD, condition_number, expm/logm/sqrtm
- Hypergeometric: `hyp0f1`, `hyp1f1`, `hyp2f0`, `hyper_pfq`, `meijerg`
- Matrix utilities: `norm_frobenius`, `is_symmetric`, `diag`, `trace`, `sub`

### Fixed
- CORDIC sin/cos angle reduction, sin() small-angle approximation
- gamma() reflection formula truncation

## [0.4.1] — 2026-06-19

### Added
- `gamma_lower`, `gamma_upper` (incomplete Gamma)
- `elliptic_k_parameter(k)` = K(k²)
- `WorkStealingQueue` parallel load balancing

### Fixed
- Bessel J₀/J₁ Taylor convergence loop
- `Mpf::cmp()` clone overhead via magnitude fast path

## [0.4.0] — 2026-06

### Added
- Miri: 0 undefined behavior
- 72h stability test harness
- `AtomicConfig<T>` generic configuration
- `Mpz::pow_big` large-exponent power
- Python bindings via PyO3 (PyMpz, PyMpf) + maturin

## [0.3.0] — 2026-06

### Added
- FFT twiddle factor cache (OnceLock)
- Thread-local limb object pool
- `AdaptivePrecision` controller
- Linear algebra: `Matrix<Mpf>`, LU decomposition
- Signal processing: FFT/iFFT, convolution, autocorrelation
- Flamegraph-compatible profiling hooks

## [0.2.0] — 2026-06

### Added
- Error convenience constructors + `with_context()` chain
- CORDIC OnceLock arctan table + adaptive iterations
- Beta function, numerical integration/differentiation/solve
- Fuzz targets in CI, comparative benchmarks vs num-bigint

## [0.1.0] — 2026-06

### Added
- Initial release: `Mpz`, `Mpf`, `Complex` types
- 6+ multiplication backends (Schoolbook → NTT)
- CORDIC-based trigonometry
- Linear algebra (Matrix, LU, QR, Cholesky, SVD)
- Signal processing (FFT/iFFT)
- Configuration system, CLI benchmarks, fuzz targets
- CI: test, clippy, fmt, fuzz
- `no_std` support via alloc
