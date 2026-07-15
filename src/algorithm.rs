//! Shared numerical algorithm helpers.
//!
//! Constants and utilities shared across arithmetic modules
//! (mpz, mpf, etc.) to avoid duplication.
//!
//! # Newton iteration
//!
//! Both [`crate::mpz`] and [`crate::mpf`] contain Newton-method implementations
//! for sqrt and reciprocal. The constants below provide canonical iteration
//! bounds and convergence thresholds. See the individual arithmetic modules
//! for the type-specific Newton loops.

/// Maximum Newton iterations for sqrt (shared by mpz and mpf).
///
/// Used as a safety bound to prevent infinite loops when convergence is slow.
// Unused until existing Newton loops migrate to these shared constants.
#[allow(dead_code)]
pub(crate) const NEWTON_SQRT_MAX_ITER: usize = 50;

/// Convergence threshold for Newton (relative change < 1ulp).
///
/// When the relative change between successive Newton iterates drops below
/// this threshold, the iteration is considered converged.
// Unused until existing Newton loops migrate to these shared constants.
#[allow(dead_code)]
pub(crate) const NEWTON_CONVERGENCE_THRESHOLD: f64 = 1e-30;
