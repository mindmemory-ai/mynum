#![no_main]

use libfuzzer_sys::fuzz_target;
use mynum::Mpz;

fuzz_target!(|data: &[u8]| {
    // 将输入解释为两个任意精度整数，检查算术不变式
    if data.len() < 2 {
        return;
    }

    // 使用不同的字节切片作为两个 Mpz 的种子
    let mid = data.len() / 2;
    let a_bytes = &data[..mid];
    let b_bytes = &data[mid..];

    // 将字节解释为无符号整数（小端序）
    let a = Mpz::from_limbs(
        a_bytes
            .chunks(8)
            .map(|chunk| {
                let mut arr = [0u8; 8];
                let len = chunk.len().min(8);
                arr[..len].copy_from_slice(&chunk[..len]);
                u64::from_le_bytes(arr)
            })
            .collect(),
        false,
    );
    let b = Mpz::from_limbs(
        b_bytes
            .chunks(8)
            .map(|chunk| {
                let mut arr = [0u8; 8];
                let len = chunk.len().min(8);
                arr[..len].copy_from_slice(&chunk[..len]);
                u64::from_le_bytes(arr)
            })
            .collect(),
        false,
    );

    // 不变式 1: a + b - b == a
    let sum = a.add(&b);
    let recovered = sum.sub(&b);
    assert_eq!(recovered, a, "addition not invertible");

    // 不变式 2: (a * b) / b == a（当 b != 0）
    if !b.is_zero() {
        let product = a.mul(&b);
        let quotient = product.div(&b).unwrap();
        assert_eq!(quotient, a, "multiplication not invertible");
    }

    // 不变式 3: a + b == b + a
    assert_eq!(a.add(&b), b.add(&a), "addition not commutative");
});
