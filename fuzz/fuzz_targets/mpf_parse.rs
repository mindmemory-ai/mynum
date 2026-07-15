#![no_main]

use libfuzzer_sys::fuzz_target;
use std::str;

fuzz_target!(|data: &[u8]| {
    // 尝试将任意字节解析为十进制字符串，不应 panic
    if let Ok(s) = str::from_utf8(data) {
        let s = s.trim();
        if s.is_empty() || s.len() > 1000 {
            return;
        }
        // from_str 绝不应 panic — 只应返回 Ok 或 Err
        let _ = mynum::Mpf::from_str(s, 10);
    }
});
