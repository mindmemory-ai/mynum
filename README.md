# MyNum — High-Precision Arbitrary-Precision Arithmetic Library

[![CI](https://github.com/gA4ss/mynum/actions/workflows/ci.yml/badge.svg)](https://github.com/gA4ss/mynum/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/mynum.svg)](https://crates.io/crates/mynum)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

MyNum is a high-performance, arbitrary-precision numeric library written in Rust.
It provides three core types — `Mpz` (arbitrary-precision integer), `Mpf`
(configurable-precision float), and `Complex` — plus 68+ special functions,
linear algebra, signal processing, and adaptive numerical integration.

## Features

### Core Types
- **`Mpz`** — arbitrary-precision integer with adaptive multiplication
  (Schoolbook → Karatsuba → Toom-Cook 3/4 → FFT → NTT)
- **`Mpf`** — configurable-precision float (`mantissa * 2^exponent`, precision
  in bits)
- **`Complex`** — complex numbers with full elementary + special function support

### Mathematical Functions
- **Elementary:** `exp`, `ln`, `sqrt`, `cbrt`, `sin`, `cos`, `tan` and inverses,
  `sinh`, `cosh`, `tanh` and inverses
- **Gamma family:** `gamma`, `loggamma`, `digamma`, `polygamma`, `gamma_lower`,
  `gamma_upper`, `beta`, `log_beta`
- **Bessel (complete):** `J₀/J₁/Jₙ`, `Y₀/Y₁/Yₙ`, `I₀/I₁/Iₙ`, `K₀/K₁/Kₙ`,
  Hankel `H¹/H²`, half-integer, real-order `Jᵥ/Yᵥ/Iᵥ/Kᵥ`
- **Airy:** `Ai`, `Bi`
- **Error:** `erf`, `erfc`, `erfi`
- **Elliptic:** `K(k)`, `E(k)`, `Π(n,m)`, `jacobi_sn/cn/dn`
- **Zeta:** `zeta`, `hurwitz_zeta`
- **Polylog / Lambert W:** `polylog`, `lambert_w`
- **Hypergeometric:** `₂F₁`, `₀F₁`, `₁F₁`, `₂F₀`, `pFq`, Meijer G
- **Integrals:** `Ei`, `Si`, `Ci`, `Shi`, `Chi`, Fresnel `S(x)`/`C(x)`
- **Advanced:** `Struve`, `Lommel`, `Whittaker`, `spherical_harmonic`, `lerchphi`

### Numerical Methods
- **Integration:** Trapezoidal, Simpson, Gauss-Legendre 5-point, adaptive
  Gauss-Kronrod G7-K15 (`quad`)
- **Series:** `nsum` with Wynn epsilon acceleration
- **Limits:** Richardson extrapolation
- **Differentiation:** forward/central difference + Richardson
- **Root finding:** bisection, Newton, secant

### Linear Algebra (`linalg`)
- `Matrix<Mpf>` with LU, QR (Householder), Cholesky, SVD, eigenvalues
- Matrix functions: `expm` (Padé + scaling/squaring), `logm`, `sqrtm`
- `norm_frobenius`, `condition_number`, `is_symmetric`, `trace`, `diag`

### Signal Processing (`signal`)
- FFT / iFFT, convolution, autocorrelation

### Other
- **Number theory:** GCD, LCM, primality testing (Miller-Rabin), factorization
  (trial division, Pollard's rho), Tonelli-Shanks
- **`no_std` compatible** — uses `alloc` as fallback when `std` is disabled
- **Configurable precision** — `GlobalPrecisionConfig` with thread-safe
  `AtomicUsize` backing
- **Adaptive precision** — `AdaptivePrecision` controller for iterative
  algorithms
- **Error recovery** — `ErrorSeverity` (Warning/Error/Fatal) +
  `RecoveryStrategy` with retry and adaptive-precision retry

## Installation

```toml
[dependencies]
mynum = "0.5"
```

MSRV: Rust 1.70+.

## Quick Start

### Basic arithmetic

```rust
use mynum::{Mpz, Mpf, Complex};

// Arbitrary-precision integers
let a = Mpz::from_str("12345678901234567890", 10).unwrap();
let b = Mpz::from_str("98765432109876543210", 10).unwrap();
let sum = a.add(&b);
println!("{} + {} = {}", a, b, sum);

// Configurable-precision floats (256 bits)
let x = Mpf::from_f64(std::f64::consts::PI, 256);
let y = Mpf::from_f64(std::f64::consts::E, 256);
println!("π × e = {}", x.mul(&y));

// Complex numbers
let z1 = Complex::from_f64(3.0, 4.0, 64);
let z2 = Complex::from_f64(1.0, 2.0, 64);
println!("(3+4i) × (1+2i) = {}", z1.mul(&z2));
```

### Special functions

```rust
use mynum::mpf::special;
use mynum::Mpf;

let x = Mpf::from_f64(1.5, 256);

// Gamma function
let g = special::gamma(&x).unwrap();
println!("Γ(1.5) = {}", g);

// Bessel J₀
let j0 = special::bessel_j0(&x).unwrap();
println!("J₀(1.5) = {}", j0);

// Riemann zeta
let z = special::zeta(&Mpf::from_f64(2.0, 256)).unwrap();
println!("ζ(2) = {}", z); // π²/6

// Error function
let e = special::erf(&Mpf::from_f64(1.0, 256)).unwrap();
println!("erf(1) = {}", e);
```

### Adaptive Gauss-Kronrod integration

```rust
use mynum::mpf::integration::quad;
use mynum::Mpf;

let result = quad(
    &|x| Ok(x.sin()?),
    &Mpf::from_f64(0.0, 128),
    &Mpf::from_f64(std::f64::consts::PI, 128),
    None,  // tolerance
    None,  // max iterations
).unwrap();
println!("∫₀^π sin(x) dx = {}", result); // ~2.0
```

### Matrix linear algebra

```rust
use mynum::linalg::Matrix;
use mynum::Mpf;

// Create a 2×2 matrix
let mut m = Matrix::new(2, 2);
m.set(0, 0, Mpf::from_f64(4.0, 128));
m.set(0, 1, Mpf::from_f64(1.0, 128));
m.set(1, 0, Mpf::from_f64(1.0, 128));
m.set(1, 1, Mpf::from_f64(3.0, 128));

// LU decomposition
let (l, u) = m.lu_decomposition().unwrap();
println!("L = {}", l);
println!("U = {}", u);
```

### Parallel multiplication

```rust
use mynum::{Mpz, MpzMultiplicationConfig};

// Enable parallel multiplication for large operands
MpzMultiplicationConfig::enable_parallel(None).unwrap();

let a = Mpz::from_str(&("1".repeat(10000)), 10).unwrap();
let b = Mpz::from_str(&("2".repeat(10000)), 10).unwrap();
let result = a.mul(&b); // uses parallel Karatsuba or FFT automatically
```

### Error recovery

```rust
use mynum::{Error, ErrorSeverity, RecoveryStrategy};

// Check severity before retrying
let err = Error::ConvergenceError("max iterations".into());
assert_eq!(err.severity(), ErrorSeverity::Error); // retryable

let err = Error::DivisionByZero;
assert_eq!(err.severity(), ErrorSeverity::Fatal); // no retry

// Retry with RecoveryStrategy
let strategy = RecoveryStrategy::new(3, 2.0);
let result = strategy.retry(|| {
    // fallible operation that may need retries
    Ok::<_, Error>(42)
});
assert_eq!(result, Ok(42));
```

## Configuration

```rust
use mynum::{MpzMultiplicationConfig, AlgorithmThresholds};

// Set custom multiplication backend thresholds
let thresholds = AlgorithmThresholds {
    schoolbook_to_karatsuba: 32,
    karatsuba_to_toom3: 256,
    toom3_to_toom4: 1024,
    toom_to_fft: 4096,
    fft_to_ntt: 16384,
};
MpzMultiplicationConfig::set_thresholds(thresholds).unwrap();

// Benchmark to find optimal thresholds on current hardware
let calibrated = MpzMultiplicationConfig::benchmark_thresholds();
MpzMultiplicationConfig::set_thresholds(calibrated).unwrap();
```

## Examples

```bash
cargo run --example basic_demo
cargo run --example advanced_math
cargo run --example threshold_calibration -- --apply
```

## Crate Features

| Feature | Default | Description |
|---------|---------|-------------|
| `std` | ✅ | Standard library support |
| `parallel` | — | Rayon-based parallelism |
| `rand-std` | ✅ | Std `OsRng` for cryptographic random |

For `no_std`, disable defaults:
```toml
[dependencies]
mynum = { version = "0.5", default-features = false }
```

## License

MIT — see [LICENSE](LICENSE).
