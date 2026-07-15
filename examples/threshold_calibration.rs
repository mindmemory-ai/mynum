//! Multiplication backend threshold calibration tool.
//!
//! Runs internal timing benchmarks to find optimal crossover points
//! between adjacent multiplication backends for the current hardware.
//!
//! Usage:
//!   cargo run --example threshold_calibration          # display calibrated thresholds
//!   cargo run --example threshold_calibration --apply  # apply to global config

use mynum::MpzMultiplicationConfig;

fn main() {
    println!("=== MyNum Threshold Calibration ===\n");

    let thresholds = MpzMultiplicationConfig::benchmark_thresholds();

    println!("Calibrated thresholds (limbs):");
    println!(
        "  Schoolbook -> Karatsuba:  {}",
        thresholds.schoolbook_to_karatsuba
    );
    println!(
        "  Karatsuba  -> Toom-Cook 3: {}",
        thresholds.karatsuba_to_toom3
    );
    println!(
        "  Toom-Cook 3 -> Toom-Cook 4: {}",
        thresholds.toom3_to_toom4
    );
    println!("  Toom-Cook 4 -> FFT:         {}", thresholds.toom_to_fft);
    println!("  FFT        -> NTT:          {}", thresholds.fft_to_ntt);
    println!();

    let args: Vec<String> = std::env::args().collect();
    if args.len() > 1 && args[1] == "--apply" {
        match MpzMultiplicationConfig::set_thresholds(thresholds) {
            Ok(_) => println!("Thresholds applied successfully."),
            Err(e) => eprintln!("Error: Failed to apply thresholds: {}", e),
        }
    } else {
        println!("Pass --apply to store these thresholds in global config.");
    }
}
