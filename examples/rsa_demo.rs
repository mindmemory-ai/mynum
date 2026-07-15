//! RSA 密钥生成和加密演示
//!
//! 使用 mynum 的高精度整数生成 RSA 密钥对并执行加解密。

use mynum::mpz::Mpz;

fn main() {
    println!("=== RSA 2048-bit 密钥生成演示 ===\n");

    // 生成两个 1024 位素数 p, q
    println!("生成 1024 位素数 p...");
    let p = Mpz::random_prime(1024).unwrap();
    println!("p = {}... ({} 位)", &p.to_string(10)[..40], p.bit_length());

    println!("生成 1024 位素数 q...");
    let q = Mpz::random_prime(1024).unwrap();
    println!("q = {}... ({} 位)", &q.to_string(10)[..40], q.bit_length());

    // n = p * q
    let n = p.mul(&q);
    println!("\nn = p × q ({} 位)", n.bit_length());

    // φ(n) = (p-1)(q-1)
    let one = Mpz::from_i64(1);
    let phi = p.sub(&one).mul(&q.sub(&one));

    // 公钥指数 e = 65537
    let e = Mpz::from_i64(65537);
    println!("公钥 e = 65537");

    // 私钥指数 d = e^(-1) mod φ(n)
    let d = e.mod_inverse(&phi).unwrap();
    println!("私钥 d 已计算 ({} 位)", d.bit_length());

    // 加密测试
    let message = Mpz::from_str("12345678901234567890", 10).unwrap();
    println!("\n原始消息: {}", message);

    let encrypted = Mpz::mod_pow(&message, &e, &n);
    println!("加密后: {}...", &encrypted.to_string(10)[..50]);

    let decrypted = Mpz::mod_pow(&encrypted, &d, &n);
    println!("解密后: {}", decrypted);

    assert_eq!(message, decrypted);
    println!("\n✅ RSA 加解密成功！");
}
