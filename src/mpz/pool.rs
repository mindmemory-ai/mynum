//! Mpz limb 对象池
//!
//! 复用 limb 向量以减少重复分配。

use std::cell::RefCell;
use std::collections::VecDeque;

const MAX_POOL_SIZE: usize = 64;
const MAX_LIMB_LEN: usize = 4096;

thread_local! {
    static LIMB_POOL: RefCell<VecDeque<Vec<u64>>> = const { RefCell::new(VecDeque::new()) };
}

/// 从池中获取一个 limb 向量（或分配新的）
pub fn acquire_limbs(capacity: usize) -> Vec<u64> {
    if capacity > MAX_LIMB_LEN {
        return Vec::with_capacity(capacity);
    }
    LIMB_POOL.with(|pool| {
        let mut pool = pool.borrow_mut();
        pool.pop_front()
            .unwrap_or_else(|| Vec::with_capacity(capacity))
    })
}

/// 将 limb 向量归还池中
pub fn release_limbs(mut limbs: Vec<u64>) {
    if limbs.capacity() > MAX_LIMB_LEN {
        return;
    }
    limbs.clear();
    LIMB_POOL.with(|pool| {
        let mut pool = pool.borrow_mut();
        if pool.len() < MAX_POOL_SIZE {
            pool.push_back(limbs);
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pool_acquire_release() {
        let limbs = acquire_limbs(64);
        let cap = limbs.capacity();
        assert!(cap >= 64);
        release_limbs(limbs);
        let reused = acquire_limbs(64);
        assert!(reused.capacity() >= 64);
    }

    #[test]
    fn test_pool_oversized_not_cached() {
        let limbs = acquire_limbs(10000);
        assert!(limbs.capacity() >= 10000);
        // 过大的向量不应被缓存（capacity > MAX_LIMB_LEN）
        release_limbs(limbs);
        // 下一个获取应该是新分配
        let fresh = acquire_limbs(64);
        assert!(fresh.capacity() >= 64);
    }
}
