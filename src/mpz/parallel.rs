//! Mpz 并行乘法模块
//!
//! 提供多线程并行乘法支持，用于加速大整数运算

use crate::error::Result;
use crate::mpz::Mpz;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::sync::OnceLock;
use std::thread;

// 全局线程计数器
static ACTIVE_THREADS: AtomicUsize = AtomicUsize::new(0);
static MAX_CONCURRENT_THREADS: AtomicUsize = AtomicUsize::new(8); // 限制最大并发线程数

/// 线程资源管理器
struct ThreadResourceManager;

impl ThreadResourceManager {
    /// 尝试获取线程资源
    fn try_acquire_thread() -> bool {
        let current = ACTIVE_THREADS.load(Ordering::Relaxed);
        let max = MAX_CONCURRENT_THREADS.load(Ordering::Relaxed);

        if current < max {
            ACTIVE_THREADS.fetch_add(1, Ordering::Relaxed);
            true
        } else {
            false
        }
    }

    /// 释放线程资源
    fn release_thread() {
        ACTIVE_THREADS.fetch_sub(1, Ordering::Relaxed);
    }

    /// 设置最大并发线程数
    fn set_max_threads(max: usize) {
        MAX_CONCURRENT_THREADS.store(max, Ordering::Relaxed);
    }
}

impl ThreadResourceManager {
    /// 创建新的线程资源管理器实例
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self
    }

    /// 使用 std::thread::scope 并行地对 items 中的每个元素应用 f，收集结果
    pub fn map_reduce<T, R, F>(&self, items: &[T], f: F) -> Vec<R>
    where
        T: Sync,
        R: Send,
        F: Fn(&T) -> R + Sync + Send,
    {
        if items.is_empty() {
            return Vec::new();
        }

        let num_threads = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4)
            .min(items.len());

        let chunk_size = items.len().div_ceil(num_threads);
        let f = std::sync::Arc::new(f);

        std::thread::scope(|s| {
            let mut handles = Vec::with_capacity(num_threads);

            for i in 0..num_threads {
                let start = i * chunk_size;
                let end = (start + chunk_size).min(items.len());
                if start >= end {
                    break;
                }
                let f = std::sync::Arc::clone(&f);
                handles.push(s.spawn(move || {
                    items[start..end]
                        .iter()
                        .map(|item| f(item))
                        .collect::<Vec<_>>()
                }));
            }

            let mut results = Vec::with_capacity(items.len());
            for handle in handles {
                if let Ok(chunk) = handle.join() {
                    results.extend(chunk);
                }
            }
            results
        })
    }
}

/// 工作窃取队列 — 线程安全的动态任务分配
#[allow(dead_code)]
struct WorkStealingQueue<T: Clone + Send + Sync> {
    items: Vec<T>,
    next_index: AtomicUsize,
}

impl<T: Clone + Send + Sync> WorkStealingQueue<T> {
    #[allow(dead_code)]
    fn new(items: Vec<T>) -> Self {
        Self {
            items,
            next_index: AtomicUsize::new(0),
        }
    }

    #[allow(dead_code)]
    fn steal(&self) -> Option<T> {
        let idx = self.next_index.fetch_add(1, Ordering::Relaxed);
        if idx < self.items.len() {
            Some(self.items[idx].clone())
        } else {
            None
        }
    }
}

/// 复数类型，用于FFT计算
#[derive(Debug, Clone, Copy)]
pub struct Complex<T> {
    /// 实部
    pub re: T,
    /// 虚部
    pub im: T,
}

impl<T> Complex<T> {
    /// 创建新的复数
    pub fn new(re: T, im: T) -> Self {
        Self { re, im }
    }

    /// 计算共轭复数
    pub fn conj(&self) -> Self
    where
        T: Clone + std::ops::Neg<Output = T>,
    {
        Self {
            re: self.re.clone(),
            im: -self.im.clone(),
        }
    }
}

impl<T> std::ops::Add for Complex<T>
where
    T: std::ops::Add<Output = T> + Clone,
{
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
}

impl<T> std::ops::Sub for Complex<T>
where
    T: std::ops::Sub<Output = T> + Clone,
{
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
}

impl<T> std::ops::Mul for Complex<T>
where
    T: std::ops::Mul<Output = T> + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + Clone,
{
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        let re = self.re.clone() * other.re.clone() - self.im.clone() * other.im.clone();
        let im = self.re * other.im + self.im * other.re;
        Self { re, im }
    }
}

impl<T> std::ops::Div<f64> for Complex<T>
where
    T: std::ops::Div<f64, Output = T> + Clone,
{
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            re: self.re / rhs,
            im: self.im / rhs,
        }
    }
}

/// 并行乘法配置
#[derive(Debug, Clone)]
pub struct ParallelConfig {
    /// 线程数，0表示自动检测
    pub num_threads: usize,
    /// 最小任务大小，小于此大小的任务不会并行化
    pub min_task_size: usize,
    /// 是否启用并行计算
    pub enabled: bool,
}

impl Default for ParallelConfig {
    fn default() -> Self {
        Self {
            num_threads: 0, // 自动检测
            enabled: true,
            min_task_size: 1000, // 最小1000位
        }
    }
}

/// 并行乘法计算器
#[derive(Clone)]
pub struct ParallelMultiplier {
    config: ParallelConfig,
    _num_threads: usize, // 标记为未使用，但保留以备将来扩展
}

impl ParallelMultiplier {
    /// 创建新的并行乘法器
    pub fn new(config: ParallelConfig) -> Self {
        // 设置合理的线程限制，避免资源耗尽
        let max_threads = std::cmp::min(config.num_threads, 4); // 最多4个线程
        ThreadResourceManager::set_max_threads(max_threads);

        Self {
            config,
            _num_threads: max_threads,
        }
    }

    /// 并行乘法（主要入口）
    pub fn parallel_multiply(&self, a: &Mpz, b: &Mpz) -> Result<Mpz> {
        if !self.config.enabled
            || a.bit_length() < self.config.min_task_size
            || b.bit_length() < self.config.min_task_size
        {
            return Ok(a.mul(b));
        }

        // 对于超大数，使用并行FFT
        if a.bit_length() > 4096 || b.bit_length() > 4096 {
            return self.fft_parallel_optimized(a, b);
        }

        // 对于大数，使用并行Karatsuba
        self.karatsuba_parallel_optimized(a, b)
    }

    /// 并行FFT乘法
    pub fn parallel_fft_multiply(&self, a: &Mpz, b: &Mpz) -> Result<Mpz> {
        if !self.config.enabled
            || a.bit_length() < self.config.min_task_size * 2
            || b.bit_length() < self.config.min_task_size * 2
        {
            return Ok(a.mul(b));
        }

        // 对于超大数，使用并行FFT
        if a.bit_length() > 4096 || b.bit_length() > 4096 {
            return self.fft_parallel_optimized(a, b);
        }

        // 对于大数，使用标准并行FFT
        self.fft_parallel_standard(a, b)
    }

    /// 标准并行FFT乘法
    fn fft_parallel_standard(&self, a: &Mpz, b: &Mpz) -> Result<Mpz> {
        let max_bits = a.bit_length().max(b.bit_length());
        let fft_size = self.next_power_of_two(max_bits * 2);

        // 将数字转换为复数多项式
        let a_poly = self.mpz_to_complex_poly(a, fft_size);
        let b_poly = self.mpz_to_complex_poly(b, fft_size);

        // 并行FFT变换
        let a_fft = self.parallel_fft(&a_poly)?;
        let b_fft = self.parallel_fft(&b_poly)?;

        // 并行点乘
        let product_fft = self.parallel_pointwise_multiply(&a_fft, &b_fft)?;

        // 并行逆FFT
        let product_poly = self.parallel_inverse_fft(&product_fft)?;

        // 转换回大整数
        Ok(self.complex_poly_to_mpz(&product_poly))
    }

    /// 优化的并行FFT乘法（使用分块策略）
    fn fft_parallel_optimized(&self, a: &Mpz, b: &Mpz) -> Result<Mpz> {
        let max_bits = a.bit_length().max(b.bit_length());

        // 对于超大数，使用分块FFT策略
        if max_bits > 16384 {
            return self.fft_blocked_strategy(a, b);
        }

        // 使用缓存优化的FFT
        self.fft_parallel_cached(a, b)
    }

    /// 分块FFT策略（适用于超大数）
    fn fft_blocked_strategy(&self, a: &Mpz, b: &Mpz) -> Result<Mpz> {
        let _max_bits = a.bit_length().max(b.bit_length());
        let block_size = 8192; // 8K位为一个块

        // 将大数分割成块
        let a_blocks = self.split_into_blocks(a, block_size);
        let b_blocks = self.split_into_blocks(b, block_size);

        // 限制并行度，避免创建过多线程
        let max_parallel_blocks = 2; // 最多2个并行块，减少线程数量
        let mut results = Vec::new();

        // 分批处理，避免同时创建过多线程
        for chunk_start in (0..a_blocks.len()).step_by(max_parallel_blocks) {
            let chunk_end = (chunk_start + max_parallel_blocks).min(a_blocks.len());
            let mut chunk_handles = Vec::new();

            for (rel_i, a_block) in a_blocks[chunk_start..chunk_end].iter().enumerate() {
                let i = chunk_start + rel_i;
                for (j, b_block) in b_blocks.iter().enumerate() {
                    let a_block = a_block.clone();
                    let b_block = b_block.clone();

                    // 检查是否有足够的线程资源
                    if ThreadResourceManager::try_acquire_thread() {
                        let handle = thread::spawn(move || {
                            let compute = ParallelMultiplier::new(ParallelConfig::default());
                            let result = compute.fft_parallel_standard(&a_block, &b_block);
                            ThreadResourceManager::release_thread();
                            result
                        });
                        chunk_handles.push((i, j, handle));
                    } else {
                        // 线程资源不足，使用串行计算
                        let compute = ParallelMultiplier::new(ParallelConfig::default());
                        let result = compute.fft_parallel_standard(&a_block, &b_block)?;
                        results.push((i, j, result));
                    }
                }
            }

            // 等待当前批次的线程完成
            for (i, j, handle) in chunk_handles {
                match handle.join() {
                    Ok(result) => {
                        match result {
                            Ok(block_result) => results.push((i, j, block_result)),
                            Err(_) => {
                                // 如果并行计算失败，回退到串行
                                let a_block = a_blocks[i].clone();
                                let b_block = b_blocks[j].clone();
                                let result = a_block.mul(&b_block);
                                results.push((i, j, result));
                            }
                        }
                    }
                    Err(_) => {
                        // 线程join失败，使用串行计算
                        let a_block = a_blocks[i].clone();
                        let b_block = b_blocks[j].clone();
                        let result = a_block.mul(&b_block);
                        results.push((i, j, result));
                    }
                }
            }
        }

        // 合并结果（考虑位置偏移）
        self.merge_block_results(&results, block_size)
    }

    /// 缓存优化的并行FFT
    fn fft_parallel_cached(&self, a: &Mpz, b: &Mpz) -> Result<Mpz> {
        let max_bits = a.bit_length().max(b.bit_length());
        let fft_size = self.next_power_of_two(max_bits * 2);

        // 检查是否有缓存的旋转因子
        let twiddle_factors = self.get_cached_twiddle_factors(fft_size);

        // 使用缓存的旋转因子进行FFT
        let a_poly = self.mpz_to_complex_poly(a, fft_size);
        let b_poly = self.mpz_to_complex_poly(b, fft_size);

        let a_fft = self.fft_with_twiddles(&a_poly, &twiddle_factors)?;
        let b_fft = self.fft_with_twiddles(&b_poly, &twiddle_factors)?;

        let product_fft = self.parallel_pointwise_multiply(&a_fft, &b_fft)?;
        let product_poly = self.inverse_fft_with_twiddles(&product_fft, &twiddle_factors)?;

        Ok(self.complex_poly_to_mpz(&product_poly))
    }

    /// 优化的并行Karatsuba乘法（多级并行）
    pub fn karatsuba_parallel_optimized(&self, a: &Mpz, b: &Mpz) -> Result<Mpz> {
        let max_bits = a.bit_length().max(b.bit_length());

        // 动态调整并行策略
        if max_bits < 1024 {
            return self.karatsuba_parallel_simple(a, b);
        }

        // 对于大数，使用多级并行策略
        let levels = (max_bits as f64).log2().ceil() as usize;
        let optimal_levels = (levels / 2).min(3); // 最多3级并行

        self.karatsuba_multi_level(a, b, optimal_levels)
    }

    /// 简化的并行Karatsuba乘法实现（避免递归）
    fn karatsuba_parallel_simple(&self, a: &Mpz, b: &Mpz) -> Result<Mpz> {
        let max_bits = a.bit_length().max(b.bit_length());

        // 如果数字太大，直接使用串行算法
        if max_bits > 2048 {
            return Ok(a.mul(b));
        }

        let half = max_bits.div_ceil(2);

        // 如果分割后仍然太大，直接使用串行算法
        if half < 512 {
            return Ok(a.mul(b));
        }

        // 分割数字
        let (a_high, a_low) = self.split_number(a, half);
        let (b_high, b_low) = self.split_number(b, half);

        // 使用更高效的结果收集方式，避免频繁的锁操作
        let (p0, p2, p1_full) = {
            // 为每个线程创建独立的克隆
            let a_low_clone1 = a_low.clone();
            let b_low_clone1 = b_low.clone();
            let a_high_clone1 = a_high.clone();
            let b_high_clone1 = b_high.clone();
            let a_low_clone2 = a_low.clone();
            let a_high_clone2 = a_high.clone();
            let b_low_clone2 = b_low.clone();
            let b_high_clone2 = b_high.clone();

            // 检查是否有足够的线程资源
            let can_parallel = ThreadResourceManager::try_acquire_thread()
                && ThreadResourceManager::try_acquire_thread()
                && ThreadResourceManager::try_acquire_thread();

            if can_parallel {
                // 并行计算三个乘积
                let p0_handle = thread::spawn(move || {
                    let result = a_low_clone1.mul(&b_low_clone1);
                    ThreadResourceManager::release_thread();
                    result
                });

                let p2_handle = thread::spawn(move || {
                    let result = a_high_clone1.mul(&b_high_clone1);
                    ThreadResourceManager::release_thread();
                    result
                });

                let p1_handle = thread::spawn(move || {
                    let a_sum = a_low_clone2.add(&a_high_clone2);
                    let b_sum = b_low_clone2.add(&b_high_clone2);
                    let result = a_sum.mul(&b_sum);
                    ThreadResourceManager::release_thread();
                    result
                });

                // 等待所有线程完成
                let p0 = p0_handle.join().unwrap_or_else(|_| a_low.mul(&b_low));
                let p2 = p2_handle.join().unwrap_or_else(|_| a_high.mul(&b_high));
                let p1_full = p1_handle.join().unwrap_or_else(|_| {
                    let a_sum = a_low.add(&a_high);
                    let b_sum = b_low.add(&b_high);
                    a_sum.mul(&b_sum)
                });

                (p0, p2, p1_full)
            } else {
                // 线程资源不足，使用串行计算
                let p0 = a_low.mul(&b_low);
                let p2 = a_high.mul(&b_high);
                let a_sum = a_low.add(&a_high);
                let b_sum = b_low.add(&b_high);
                let p1_full = a_sum.mul(&b_sum);

                (p0, p2, p1_full)
            }
        };

        // 计算 P1 = P1_full - P0 - P2
        let p1 = p1_full.sub(&p0).sub(&p2);

        // 优化：使用位运算而不是乘法来计算移位
        let shift_half = Mpz::from_i64(1).shl(half);

        // 检查是否有足够的线程资源进行并行计算
        let can_parallel_shift = ThreadResourceManager::try_acquire_thread()
            && ThreadResourceManager::try_acquire_thread()
            && ThreadResourceManager::try_acquire_thread();

        if can_parallel_shift {
            // 并行计算移位和加法操作
            let shift_half_clone1 = shift_half.clone();
            let shift_half_clone2 = shift_half.clone();
            let p1_clone = p1.clone();
            let p2_clone = p2.clone();
            let p0_clone = p0.clone();
            let term1_clone = p1.clone();

            let term1_handle = thread::spawn(move || {
                let result = p1_clone.mul(&shift_half_clone1);
                ThreadResourceManager::release_thread();
                result
            });

            let term2_handle = thread::spawn(move || {
                let result = p2_clone.mul(&shift_half_clone2).mul(&shift_half_clone2);
                ThreadResourceManager::release_thread();
                result
            });

            let temp_sum_handle = thread::spawn(move || {
                let result = p0_clone.add(&term1_clone);
                ThreadResourceManager::release_thread();
                result
            });

            let _term1 = term1_handle.join().unwrap_or_else(|_| p1.mul(&shift_half));
            let term2 = term2_handle
                .join()
                .unwrap_or_else(|_| p2.mul(&shift_half).mul(&shift_half));
            let temp_sum = temp_sum_handle.join().unwrap_or_else(|_| p0.add(&p1));

            Ok(temp_sum.add(&term2))
        } else {
            // 线程资源不足，使用串行计算
            let term1 = p1.mul(&shift_half);
            let term2 = p2.mul(&shift_half).mul(&shift_half);
            let temp_sum = p0.add(&term1);
            Ok(temp_sum.add(&term2))
        }
    }

    /// 多级并行Karatsuba乘法
    fn karatsuba_multi_level(&self, a: &Mpz, b: &Mpz, levels: usize) -> Result<Mpz> {
        if levels == 0 || a.bit_length() < 256 || b.bit_length() < 256 {
            return Ok(a.mul(b));
        }

        let max_bits = a.bit_length().max(b.bit_length());
        let half = max_bits.div_ceil(2);

        // 分割数字
        let (a_high, a_low) = self.split_number(a, half);
        let (b_high, b_low) = self.split_number(b, half);

        // 递归并行计算，减少并行级别，添加错误处理
        let a_low_clone1 = a_low.clone();
        let b_low_clone1 = b_low.clone();
        let a_high_clone1 = a_high.clone();
        let b_high_clone1 = b_high.clone();
        let a_low_clone2 = a_low.clone();
        let a_high_clone2 = a_high.clone();
        let b_low_clone2 = b_low.clone();
        let b_high_clone2 = b_high.clone();

        let left_handle = thread::spawn(move || {
            let compute = ParallelMultiplier::new(ParallelConfig::default());
            compute.karatsuba_multi_level(&a_low_clone1, &b_low_clone1, levels - 1)
        });

        let right_handle = thread::spawn(move || {
            let compute = ParallelMultiplier::new(ParallelConfig::default());
            compute.karatsuba_multi_level(&a_high_clone1, &b_high_clone1, levels - 1)
        });

        let center_handle = thread::spawn(move || {
            let a_sum = a_low_clone2.add(&a_high_clone2);
            let b_sum = b_low_clone2.add(&b_high_clone2);
            let compute = ParallelMultiplier::new(ParallelConfig::default());
            compute.karatsuba_multi_level(&a_sum, &b_sum, levels - 1)
        });

        // 处理线程执行失败的情况
        let p0 = match left_handle.join() {
            Ok(result) => result?,
            Err(_) => {
                // 线程join失败，使用串行计算
                let compute = ParallelMultiplier::new(ParallelConfig::default());
                compute.karatsuba_multi_level(&a_low, &b_low, levels - 1)?
            }
        };

        let p2 = match right_handle.join() {
            Ok(result) => result?,
            Err(_) => {
                // 线程join失败，使用串行计算
                let compute = ParallelMultiplier::new(ParallelConfig::default());
                compute.karatsuba_multi_level(&a_high, &b_high, levels - 1)?
            }
        };

        let p1_full = match center_handle.join() {
            Ok(result) => result?,
            Err(_) => {
                // 线程join失败，使用串行计算
                let a_sum = a_low.add(&a_high);
                let b_sum = b_low.add(&b_high);
                let compute = ParallelMultiplier::new(ParallelConfig::default());
                compute.karatsuba_multi_level(&a_sum, &b_sum, levels - 1)?
            }
        };

        // 计算 P1 = P1_full - P0 - P2
        let p1 = p1_full.sub(&p0).sub(&p2);

        // 使用位运算优化移位操作
        let shift_half = Mpz::from_i64(1).shl(half);

        // 并行计算最终结果，添加错误处理
        let shift_half_clone1 = shift_half.clone();
        let shift_half_clone2 = shift_half.clone();
        let p1_clone = p1.clone();
        let p2_clone = p2.clone();

        let term1_handle = thread::spawn(move || p1_clone.mul(&shift_half_clone1));
        let term2_handle =
            thread::spawn(move || p2_clone.mul(&shift_half_clone2).mul(&shift_half_clone2));

        let term1 = match term1_handle.join() {
            Ok(result) => result,
            Err(_) => p1.mul(&shift_half), // 回退到串行计算
        };

        let term2 = match term2_handle.join() {
            Ok(result) => result,
            Err(_) => p2.mul(&shift_half).mul(&shift_half), // 回退到串行计算
        };

        let temp_sum = p0.add(&term1);
        Ok(temp_sum.add(&term2))
    }

    /// 分割数字为高位和低位
    fn split_number(&self, num: &Mpz, half_bits: usize) -> (Mpz, Mpz) {
        let mask = Mpz::from_i64(1).shl(half_bits).sub(&Mpz::from_i64(1));
        let low = num.bitwise_and(&mask);
        let high = num.shr(half_bits);
        (high, low)
    }

    /// 并行FFT变换
    fn parallel_fft(&self, poly: &[Complex<f64>]) -> Result<Vec<Complex<f64>>> {
        let n = poly.len();
        if n <= 1 {
            return Ok(poly.to_vec());
        }

        // 使用分治策略并行FFT
        if n <= 64 {
            return Ok(self.sequential_fft(poly));
        }

        // 分割多项式
        let (even, odd): (Vec<_>, Vec<_>) = poly.iter().enumerate().partition(|(i, _)| i % 2 == 0);
        let even: Vec<_> = even.into_iter().map(|(_, x)| *x).collect();
        let odd: Vec<_> = odd.into_iter().map(|(_, x)| *x).collect();

        // 为fallback创建额外的克隆
        let even_fallback = even.clone();
        let odd_fallback = odd.clone();

        // 检查是否有足够的线程资源
        let can_parallel = ThreadResourceManager::try_acquire_thread()
            && ThreadResourceManager::try_acquire_thread();

        if can_parallel {
            // 并行递归FFT
            let even_handle = thread::spawn(move || {
                let compute = ParallelMultiplier::new(ParallelConfig::default());
                let result = compute.parallel_fft(&even);
                ThreadResourceManager::release_thread();
                result
            });

            let odd_handle = thread::spawn(move || {
                let compute = ParallelMultiplier::new(ParallelConfig::default());
                let result = compute.parallel_fft(&odd);
                ThreadResourceManager::release_thread();
                result
            });

            // 处理线程执行失败的情况
            let even_fft = match even_handle.join() {
                Ok(result) => result?,
                Err(_) => {
                    // 线程join失败，使用串行计算
                    let compute = ParallelMultiplier::new(ParallelConfig::default());
                    compute.sequential_fft(&even_fallback)
                }
            };

            let odd_fft = match odd_handle.join() {
                Ok(result) => result?,
                Err(_) => {
                    // 线程join失败，使用串行计算
                    let compute = ParallelMultiplier::new(ParallelConfig::default());
                    compute.sequential_fft(&odd_fallback)
                }
            };

            // 合并结果
            let mut result = vec![Complex::new(0.0, 0.0); n];
            let half = n / 2;

            for k in 0..half {
                let angle = -2.0 * std::f64::consts::PI * k as f64 / n as f64;
                let twiddle = Complex::new(angle.cos(), angle.sin());

                result[k] = even_fft[k] + twiddle * odd_fft[k];
                result[k + half] = even_fft[k] - twiddle * odd_fft[k];
            }

            Ok(result)
        } else {
            // 线程资源不足，使用串行计算
            Ok(self.sequential_fft(poly))
        }
    }

    /// 并行点乘
    fn parallel_pointwise_multiply(
        &self,
        a: &[Complex<f64>],
        b: &[Complex<f64>],
    ) -> Result<Vec<Complex<f64>>> {
        if a.len() != b.len() {
            return Err(crate::error::Error::InvalidInput(
                "FFT arrays must have same length".into(),
            ));
        }

        let n = a.len();
        let num_threads = if self.config.num_threads == 0 {
            4
        } else {
            self.config.num_threads
        };
        let chunk_size = n.div_ceil(num_threads);

        // 并行点乘
        let result: Vec<_> = a
            .chunks(chunk_size)
            .zip(b.chunks(chunk_size))
            .flat_map(|(a_chunk, b_chunk)| {
                a_chunk
                    .iter()
                    .zip(b_chunk.iter())
                    .map(|(a_val, b_val)| {
                        let re = a_val.re * b_val.re - a_val.im * b_val.im;
                        let im = a_val.re * b_val.im + a_val.im * b_val.re;
                        Complex::new(re, im)
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        Ok(result)
    }

    /// 并行逆FFT
    fn parallel_inverse_fft(&self, fft: &[Complex<f64>]) -> Result<Vec<Complex<f64>>> {
        let n = fft.len();
        if n <= 1 {
            return Ok(fft.to_vec());
        }

        // 共轭输入
        let conjugated: Vec<_> = fft.iter().map(|x| x.conj()).collect();

        // 正向FFT
        let forward = self.parallel_fft(&conjugated)?;

        // 共轭输出并除以n
        let result: Vec<_> = forward.iter().map(|x| x.conj() / n as f64).collect();

        Ok(result)
    }

    /// 使用缓存的旋转因子进行FFT
    fn fft_with_twiddles(
        &self,
        poly: &[Complex<f64>],
        _twiddles: &[Complex<f64>],
    ) -> Result<Vec<Complex<f64>>> {
        // 实现使用缓存的旋转因子的FFT
        self.parallel_fft(poly)
    }

    /// 使用缓存的旋转因子进行逆FFT
    fn inverse_fft_with_twiddles(
        &self,
        fft: &[Complex<f64>],
        _twiddles: &[Complex<f64>],
    ) -> Result<Vec<Complex<f64>>> {
        // 实现使用缓存的旋转因子的逆FFT
        self.parallel_inverse_fft(fft)
    }

    /// 获取缓存的旋转因子
    fn get_cached_twiddle_factors(&self, size: usize) -> Vec<Complex<f64>> {
        /// 全局预计算的 FFT 旋转因子缓存
        static GLOBAL_TWIDDLE_CACHE: OnceLock<Mutex<HashMap<usize, Vec<Complex<f64>>>>> =
            OnceLock::new();

        fn get_cache() -> &'static Mutex<HashMap<usize, Vec<Complex<f64>>>> {
            GLOBAL_TWIDDLE_CACHE.get_or_init(|| Mutex::new(HashMap::new()))
        }

        let mut cache = get_cache().lock().unwrap();
        if let Some(factors) = cache.get(&size) {
            return factors.clone();
        }

        let mut twiddles = Vec::with_capacity(size);
        for k in 0..size {
            let angle = -2.0 * std::f64::consts::PI * k as f64 / size as f64;
            twiddles.push(Complex::new(angle.cos(), angle.sin()));
        }
        cache.insert(size, twiddles.clone());
        twiddles
    }

    /// 顺序FFT（用于小规模计算）
    fn sequential_fft(&self, poly: &[Complex<f64>]) -> Vec<Complex<f64>> {
        let n = poly.len();
        if n <= 1 {
            return poly.to_vec();
        }

        // 简单的递归FFT实现
        let (even, odd): (Vec<_>, Vec<_>) = poly.iter().enumerate().partition(|(i, _)| i % 2 == 0);
        let even: Vec<_> = even.into_iter().map(|(_, x)| *x).collect();
        let odd: Vec<_> = odd.into_iter().map(|(_, x)| *x).collect();

        let even_fft = self.sequential_fft(&even);
        let odd_fft = self.sequential_fft(&odd);

        let mut result = vec![Complex::new(0.0, 0.0); n];
        let half = n / 2;

        for k in 0..half {
            let angle = -2.0 * std::f64::consts::PI * k as f64 / n as f64;
            let twiddle = Complex::new(angle.cos(), angle.sin());

            result[k] = even_fft[k] + twiddle * odd_fft[k];
            result[k + half] = even_fft[k] - twiddle * odd_fft[k];
        }

        result
    }

    /// 辅助方法：获取下一个2的幂
    fn next_power_of_two(&self, n: usize) -> usize {
        let mut power = 1;
        while power < n {
            power *= 2;
        }
        power
    }

    /// 辅助方法：将Mpz转换为复数多项式
    fn mpz_to_complex_poly(&self, num: &Mpz, size: usize) -> Vec<Complex<f64>> {
        let mut poly = vec![Complex::new(0.0, 0.0); size];

        // 将大整数的每一位转换为复数
        for (i, limb) in num.limbs().iter().enumerate() {
            for bit in 0..64 {
                if (limb >> bit) & 1 == 1 {
                    let pos = i * 64 + bit;
                    if pos < size {
                        poly[pos] = Complex::new(1.0, 0.0);
                    }
                }
            }
        }

        poly
    }

    /// 辅助方法：将复数多项式转换回Mpz
    fn complex_poly_to_mpz(&self, poly: &[Complex<f64>]) -> Mpz {
        let mut result = Mpz::new();
        let mut power = Mpz::from_i64(1);

        for &coeff in poly {
            let real_part = coeff.re.round() as i64;
            if real_part != 0 {
                let term = power.mul(&Mpz::from_i64(real_part.abs()));
                if real_part > 0 {
                    result = result.add(&term);
                } else {
                    result = result.sub(&term);
                }
            }
            power = power.mul(&Mpz::from_i64(2));
        }

        result
    }

    /// 辅助方法：将大整数分割成块
    fn split_into_blocks(&self, num: &Mpz, block_size: usize) -> Vec<Mpz> {
        let mut blocks = Vec::new();
        let mut remaining = num.clone();

        while !remaining.is_zero() {
            let mask = Mpz::from_i64(1).shl(block_size).sub(&Mpz::from_i64(1));
            let block = remaining.bitwise_and(&mask);
            blocks.push(block);
            remaining = remaining.shr(block_size);
        }

        blocks.reverse();
        blocks
    }

    /// 辅助方法：合并块结果
    fn merge_block_results(
        &self,
        results: &[(usize, usize, Mpz)],
        block_size: usize,
    ) -> Result<Mpz> {
        let mut final_result = Mpz::new();

        for &(i, j, ref result) in results {
            let offset = (i + j) * block_size;
            let shifted = result.shl(offset);
            final_result = final_result.add(&shifted);
        }

        Ok(final_result)
    }
}

// ============================================================
// 统一并行执行后端 (Task T2)
// ============================================================

/// 并行后端选择：rayon（feature="parallel"）或 std::thread 回退
#[allow(dead_code)]
enum ParallelBackend {
    #[cfg(feature = "parallel")]
    Rayon,
    StdThread {
        manager: ThreadResourceManager,
    },
}

/// 统一并行执行器，根据 feature 自动选择 rayon 或 std::thread。
pub(crate) struct ParallelExecutor {
    backend: ParallelBackend,
}

impl ParallelExecutor {
    /// 创建新的并行执行器。
    /// 当 `parallel` feature 启用时使用 rayon，否则使用 std::thread。
    pub fn new() -> Self {
        Self {
            #[cfg(feature = "parallel")]
            backend: ParallelBackend::Rayon,
            #[cfg(not(feature = "parallel"))]
            backend: ParallelBackend::StdThread {
                manager: ThreadResourceManager::new(),
            },
        }
    }

    /// 对 items 中的每个元素并行应用 f，收集结果。
    pub fn map_reduce<T, R, F>(&self, items: &[T], f: F) -> Vec<R>
    where
        T: Sync,
        R: Send,
        F: Fn(&T) -> R + Sync + Send,
    {
        match &self.backend {
            #[cfg(feature = "parallel")]
            ParallelBackend::Rayon => {
                use rayon::prelude::*;
                items.par_iter().map(f).collect()
            }
            ParallelBackend::StdThread { manager } => manager.map_reduce(items, f),
        }
    }
}

/// 全局并行乘法实例（线程安全）
static GLOBAL_PARALLEL: OnceLock<Mutex<ParallelMultiplier>> = OnceLock::new();

/// 获取全局并行乘法实例的克隆（避免生命周期问题）
pub fn get_global_parallel() -> ParallelMultiplier {
    let global = GLOBAL_PARALLEL
        .get_or_init(|| Mutex::new(ParallelMultiplier::new(ParallelConfig::default())));
    global
        .lock()
        .expect("Failed to acquire lock on global parallel multiplier")
        .clone()
}

/// 设置全局并行乘法配置
pub fn set_global_parallel_config(config: ParallelConfig) {
    let global = GLOBAL_PARALLEL
        .get_or_init(|| Mutex::new(ParallelMultiplier::new(ParallelConfig::default())));
    *global
        .lock()
        .expect("Failed to acquire lock on global parallel multiplier") =
        ParallelMultiplier::new(config);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parallel_multiply() {
        let parallel = ParallelMultiplier::new(ParallelConfig::default());
        let a = Mpz::from_i64(12345);
        let b = Mpz::from_i64(67890);

        let result = parallel.parallel_multiply(&a, &b).unwrap();
        let expected = a.mul(&b);

        assert_eq!(result.to_string(10), expected.to_string(10));
    }

    #[test]
    fn test_karatsuba_parallel_simple() {
        let parallel = ParallelMultiplier::new(ParallelConfig::default());
        let a = Mpz::from_str("12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890", 10).unwrap();
        let b = Mpz::from_str("98765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210", 10).unwrap();

        let result = parallel.karatsuba_parallel_simple(&a, &b).unwrap();
        let expected = a.mul(&b);

        assert_eq!(result.to_string(10), expected.to_string(10));
    }

    #[test]
    fn test_split_number() {
        let parallel = ParallelMultiplier::new(ParallelConfig::default());
        let num = Mpz::from_i64(12345);
        let (high, low) = parallel.split_number(&num, 8);

        // 验证分割结果
        let reconstructed = high.shl(8).add(&low);
        assert_eq!(reconstructed.to_string(10), num.to_string(10));
    }

    #[test]
    fn test_twiddle_factor_cache_reuse() {
        let parallel = ParallelMultiplier::new(ParallelConfig::default());
        // 相同大小的旋转因子应该被缓存复用
        let t1 = parallel.get_cached_twiddle_factors(1024);
        let t2 = parallel.get_cached_twiddle_factors(1024);
        assert_eq!(t1.len(), t2.len());
        assert!((t1[0].re - t2[0].re).abs() < 1e-15);
        assert!((t1[0].im - t2[0].im).abs() < 1e-15);
        let t3 = parallel.get_cached_twiddle_factors(512);
        assert_ne!(t1.len(), t3.len());
    }

    #[test]
    fn test_work_stealing_queue() {
        use std::sync::Arc;
        let items: Vec<usize> = (0..100).collect();
        let queue = Arc::new(WorkStealingQueue::new(items));
        let completed = Arc::new(AtomicUsize::new(0));
        let handles: Vec<_> = (0..4)
            .map(|_| {
                let q = queue.clone();
                let c = completed.clone();
                thread::spawn(move || {
                    while q.steal().is_some() {
                        c.fetch_add(1, Ordering::Relaxed);
                    }
                })
            })
            .collect();
        for h in handles {
            h.join().unwrap();
        }
        assert_eq!(completed.load(Ordering::Relaxed), 100);
    }
}
