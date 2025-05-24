#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <pthread.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <thread> 
using namespace boost::multiprecision;
// 可以自行添加需要的头文件
cpp_int mod_pow_big(cpp_int a, cpp_int b, cpp_int p);
cpp_int mod_inverse_big(cpp_int a, cpp_int p);
void fRead(int *a, int *b, int *n, int *p, int input_id){
    // 数据输入函数
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin>>*n>>*p;
    for (int i = 0; i < *n; i++){
        fin>>a[i];
    }
    for (int i = 0; i < *n; i++){   
        fin>>b[i];
    }
}

void fCheck(int *ab, int n, int input_id){
    // 判断多项式乘法结果是否正确
    std::string str1 = "/nttdata/"; // 确保此路径前缀正确
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';

    std::ifstream fin;
    fin.open(data_path, std::ios::in);

    if (!fin.is_open()) {
        std::cout << "错误: 无法打开预期的输出文件 " << data_path << std::endl;
        return;
    }

    bool error_found = false;
    int error_count = 0;
    const int MAX_ERRORS_TO_PRINT = 5; // 限制打印的错误数量

    for (int i = 0; i < n * 2 - 1; i++){
        int expected_val;
        if (!(fin >> expected_val)) {
            std::cout << "错误: 无法从文件读取索引 " << i << " 处的预期值。" << std::endl;
            error_found = true;
            break; // 文件读取问题，终止检查
        }
        if(expected_val != ab[i]){
            if (error_count < MAX_ERRORS_TO_PRINT) {
                std::cout << "多项式乘法结果错误于索引 " << i << ". 预期值: " << expected_val << ", 实际值: " << ab[i] << std::endl;
            }
            error_found = true;
            error_count++;
        }
    }
    fin.close();

    if (!error_found) {
        std::cout << "多项式乘法结果正确" << std::endl;
    } else {
        std::cout << "多项式乘法检测完毕，共发现 " << error_count << " 处错误。" << std::endl;
        if (error_count > MAX_ERRORS_TO_PRINT) {
            std::cout << "(只显示了前 " << MAX_ERRORS_TO_PRINT << " 个错误)" << std::endl;
        }
    }
    return;
}

void fWrite(int *ab, int n, int input_id){
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++){
        fout<<ab[i]<<'\n';
    }
}
// 朴素多项式乘法
void poly_multiply(int *a, int *b, int *ab, int n, int p){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            ab[i+j]=(1LL * a[i] * b[j] % p + ab[i+j]) % p;
        }
    }
}
//============================ 基础函数 ============================

// 快速幂 a^b mod p
long long mod_pow(long long a, long long b, int p) {
    long long res = 1;
    a %= p;
    while (b > 0) {
        if (b & 1) res = res * a % p;
        a = a * a % p;
        b >>= 1;
    }
    return res;
}

// 求逆元 a^(p-2) mod p
long long mod_inverse(long long a, int p) {
    return mod_pow(a, p - 2, p);
}

// 计算大于等于n的最小的2的幂
int get_min_power_of_2(int n) {
    int res = 1;
    while (res < n) res <<= 1;
    return res;
}

// 位逆序置换
void bit_reverse(int *a, int len) {
    for (int i = 1, j = len / 2, k; i < len - 1; i++) {
        if (i < j) std::swap(a[i], a[j]);
        for (k = len / 2; j >= k; j -= k, k /= 2);
        j += k;
    }
}
//============================ 朴素NTT ============================

// 朴素NTT
void ntt(int *a, int len, int type, int p) {
    bit_reverse(a, len);
    
    for (int h = 2; h <= len; h <<= 1) {
        int gn = mod_pow(3, (p - 1) / h, p);  // 原根
        if (type == -1) gn = mod_inverse(gn, p);  // INTT
        
        for (int j = 0; j < len; j += h) {
            long long g = 1;
            for (int k = j; k < j + h / 2; k++) {
                long long u = a[k];
                long long v = g * a[k + h / 2] % p;
                a[k] = (u + v) % p;
                a[k + h / 2] = (u - v + p) % p;
                g = g * gn % p;
            }
        }
    }
    
    if (type == -1) {
        long long inv_len = mod_inverse(len, p);
        for (int i = 0; i < len; i++) {
            a[i] = a[i] * inv_len % p;
        }
    }
}

// 序列填充和NTT多项式乘法
void ntt_multiply(int *a, int *b, int *c, int n, int p) {
    int len = get_min_power_of_2(2 * n - 1);
    
    // 创建临时数组并填充0
    int *ta = new int[len];
    int *tb = new int[len];
    std::fill(ta, ta + len, 0);
    std::fill(tb, tb + len, 0);
    
    // 复制原始数据
    for (int i = 0; i < n; i++) {
        ta[i] = a[i];
        tb[i] = b[i];
    }
    
    // 进行NTT变换
    ntt(ta, len, 1, p);
    ntt(tb, len, 1, p);
    
    // 点值乘法
    for (int i = 0; i < len; i++) {
        ta[i] = 1LL * ta[i] * tb[i] % p;
    }
    
    // 进行INTT变换
    ntt(ta, len, -1, p);
    
    // 结果赋值给c
    for (int i = 0; i < 2 * n - 1; i++) {
        c[i] = ta[i];
    }
    
    delete[] ta;
    delete[] tb;
}

//============================ pthread多线程实现 ============================
// pthread参数结构体
struct ThreadArgs {
    int *a;             // 数组指针
    int len;            // 数组长度
    int h;              // 当前层级(2^h) or m for four_div
    int p;              // 模数
    long long gn;       // 当前单位原根 or type for four_div or inv_len for scaling
    pthread_barrier_t *barrier; // 同步屏障
    int thread_id;      // 线程ID (can be repurposed as start_index for some tasks)
    int num_threads;    // 总线程数 (can be repurposed as count/step for some tasks)
    // 为四分NTT和特定任务添加新成员
    int start_index;    // 通用起始索引
    int work_amount;    // 通用工作量或步长
};

// pthread朴素NTT线程函数
void* pthread_ntt_worker(void *arg) {
    ThreadArgs *args = (ThreadArgs*)arg;
    int *a = args->a;
    int len = args->len;
    int h = args->h;
    int p = args->p;
    long long gn = args->gn;
    int thread_id = args->thread_id;
    int num_threads = args->num_threads;

    // 每个线程处理 len/h 个蝶形运算块中的一部分
    // 每个块大小为 h
    // 总共有 len/h 个块
    // 每个线程负责处理大约 (len/h) / num_threads 个块
    // 或者更简单地，让每个线程处理其分配到的 j 值
    for (int j = thread_id * h; j < len; j += num_threads * h) {
        long long g = 1;
        for (int k = j; k < j + h / 2; k++) {
            if (k + h/2 >= len) continue; // 边界检查
            long long u = a[k];
            long long v = g * a[k + h / 2] % p;
            a[k] = (u + v) % p;
            a[k + h / 2] = (u - v + p) % p;
            g = g * gn % p;
        }
    }
    
    // 等待所有线程完成当前级别处理
    if (args->barrier) { // 只有在正向/逆向变换的主循环中才使用barrier
        pthread_barrier_wait(args->barrier);
    }
    return NULL;
}

// pthread朴素NTT实现
void pthread_ntt(int *a, int len, int type, int p, int num_threads_requested) {
    if (num_threads_requested <= 0) num_threads_requested = 1;
    int num_threads = std::min(num_threads_requested, len / 2); // 限制线程数，避免过多线程处理小任务
    if (num_threads < 1) num_threads = 1;


    bit_reverse(a, len);
    
    pthread_t* threads = new pthread_t[num_threads];
    ThreadArgs* args = new ThreadArgs[num_threads];
    pthread_barrier_t barrier;
    
    bool barrier_initialized = false;
    if (num_threads > 1) {
        if (pthread_barrier_init(&barrier, NULL, num_threads) != 0) {
            std::cerr << "pthread_barrier_init failed" << std::endl;
            delete[] threads;
            delete[] args;
            // Fallback to serial NTT if barrier fails for multithreading
            ntt(a, len, type, p);
            return;
        }
        barrier_initialized = true;
    }
    
    for (int h = 2; h <= len; h <<= 1) {
        long long gn_h = mod_pow(3, (p - 1) / h, p);
        if (type == -1) gn_h = mod_inverse(gn_h, p);
        
        if (num_threads > 1) {
            for (int t = 0; t < num_threads; t++) {
                args[t] = {a, len, h, p, gn_h, &barrier, t, num_threads};
                if (pthread_create(&threads[t], NULL, pthread_ntt_worker, &args[t]) != 0) {
                     std::cerr << "Failed to create thread " << t << " in main loop" << std::endl;
                     // Handle thread creation failure, e.g., cleanup and fallback
                     for(int i=0; i<t; ++i) pthread_join(threads[i], NULL);
                     if(barrier_initialized) pthread_barrier_destroy(&barrier);
                     delete[] threads;
                     delete[] args;
                     ntt(a,len,type,p); // Fallback
                     return;
                }
            }
            
            for (int t = 0; t < num_threads; t++) {
                pthread_join(threads[t], NULL);
            }
        } else { // Serial execution for this h level if num_threads is 1
            for (int j = 0; j < len; j += h) {
                long long g = 1;
                for (int k = j; k < j + h / 2; k++) {
                    long long u = a[k];
                    long long v = g * a[k + h / 2] % p;
                    a[k] = (u + v) % p;
                    a[k + h / 2] = (u - v + p) % p;
                    g = g * gn_h % p;
                }
            }
        }
    }
    
    if (type == -1) {
        long long inv_len_val = mod_inverse(len, p);
        if (num_threads > 1) {
            // Lambda for final scaling
            auto scaling_worker_lambda = [](void *arg_lambda) -> void* {
                ThreadArgs *args_lambda = (ThreadArgs*)arg_lambda;
                int *arr_lambda = args_lambda->a;
                int len_lambda = args_lambda->len;
                int p_lambda = args_lambda->p;
                long long inv_val_lambda = args_lambda->gn; // Reusing gn for inv_len
                int thread_id_lambda = args_lambda->thread_id;
                int num_threads_lambda = args_lambda->num_threads;

                int chunk_size = (len_lambda + num_threads_lambda - 1) / num_threads_lambda;
                int start_idx = thread_id_lambda * chunk_size;
                int end_idx = std::min(start_idx + chunk_size, len_lambda);

                for (int i = start_idx; i < end_idx; i++) {
                    arr_lambda[i] = (1LL * arr_lambda[i] * inv_val_lambda) % p_lambda;
                    if (arr_lambda[i] < 0) arr_lambda[i] += p_lambda;
                }
                return NULL;
            };

            for (int t = 0; t < num_threads; t++) {
                // Pass inv_len_val via gn field. Set barrier to NULL as it's not used here.
                args[t] = {a, len, 0, p, inv_len_val, NULL, t, num_threads};
                if(pthread_create(&threads[t], NULL, scaling_worker_lambda, &args[t])!=0){
                    std::cerr << "Failed to create thread " << t << " in scaling" << std::endl;
                    for(int i=0; i<t; ++i) pthread_join(threads[i], NULL);
                    if(barrier_initialized) pthread_barrier_destroy(&barrier);
                    delete[] threads;
                    delete[] args;
                    // Fallback for scaling part
                    for(int i=0; i<len; ++i) a[i] = (1LL*a[i]*inv_len_val)%p;
                    return;
                }
            }
            
            for (int t = 0; t < num_threads; t++) {
                pthread_join(threads[t], NULL);
            }
        } else { // Serial scaling if num_threads is 1
            for (int i = 0; i < len; i++) {
                a[i] = (1LL * a[i] * inv_len_val) % p;
                 if (a[i] < 0) a[i] += p;
            }
        }
    }
    
    if (barrier_initialized) {
        pthread_barrier_destroy(&barrier);
    }
    delete[] threads;
    delete[] args;
}

// pthread朴素NTT多项式乘法
void pthread_ntt_multiply(int *a_orig, int *b_orig, int *c_final, int n, int p, int num_threads_requested) {
    if (num_threads_requested <= 0) num_threads_requested = 1;
    int len = get_min_power_of_2(2 * n - 1);
    
    int *ta = new int[len];
    int *tb = new int[len];
    // Initialize ta and tb in parallel if num_threads > 1
    if (num_threads_requested > 1) {
        std::thread t_init_a([&]() {
            std::fill(ta, ta + len, 0);
            std::copy(a_orig, a_orig + n, ta);
        });
        std::thread t_init_b([&]() {
            std::fill(tb, tb + len, 0);
            std::copy(b_orig, b_orig + n, tb);
        });
        t_init_a.join();
        t_init_b.join();
    } else {
        std::fill(ta, ta + len, 0);
        std::copy(a_orig, a_orig + n, ta);
        std::fill(tb, tb + len, 0);
        std::copy(b_orig, b_orig + n, tb);
    }
    
    pthread_ntt(ta, len, 1, p, num_threads_requested);
    pthread_ntt(tb, len, 1, p, num_threads_requested);
    
    // Pointwise multiplication
    if (num_threads_requested > 1) {
        std::vector<std::thread> pointwise_threads;
        int num_actual_threads_pm = std::min(num_threads_requested, len);
        if(num_actual_threads_pm < 1) num_actual_threads_pm = 1;

        for (int t = 0; t < num_actual_threads_pm; t++) {
            pointwise_threads.emplace_back([t, num_actual_threads_pm, len, ta, tb, p]() {
                int chunk_size = (len + num_actual_threads_pm - 1) / num_actual_threads_pm;
                int start_idx = t * chunk_size;
                int end_idx = std::min(start_idx + chunk_size, len);
                for (int i = start_idx; i < end_idx; i++) {
                    ta[i] = (1LL * ta[i] * tb[i]) % p;
                }
            });
        }
        for (auto& th : pointwise_threads) {
            th.join();
        }
    } else {
        for (int i = 0; i < len; i++) {
            ta[i] = (1LL * ta[i] * tb[i]) % p;
        }
    }
    
    pthread_ntt(ta, len, -1, p, num_threads_requested);
    
    std::copy(ta, ta + (2 * n - 1), c_final);
    
    delete[] ta;
    delete[] tb;
}
//============================ OpenMP多线程实现 ============================
// OpenMP朴素NTT实现
void openmp_ntt(int *a, int len, int type, int p, int num_threads) {
    bit_reverse(a, len);
    
    for (int h = 2; h <= len; h <<= 1) {
        int gn = mod_pow(3, (p - 1) / h, p);
        if (type == -1) gn = mod_inverse(gn, p);
        
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (int j = 0; j < len; j += h) {
            long long g = 1;
            for (int k = j; k < j + h / 2; k++) {
                long long u = a[k];
                long long v = g * a[k + h / 2] % p;
                a[k] = (u + v) % p;
                a[k + h / 2] = (u - v + p) % p;
                g = g * gn % p;
            }
        }
        // OpenMP隐式同步点 - 所有线程完成当前层处理后才能进入下一层
    }
    
    if (type == -1) {
        long long inv_len = mod_inverse(len, p);
        #pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < len; i++) {
            a[i] = a[i] * inv_len % p;
        }
    }
}

// OpenMP朴素NTT多项式乘法
void openmp_ntt_multiply(int *a, int *b, int *c, int n, int p, int num_threads) {
    int len = get_min_power_of_2(2 * n - 1);
    
    int *ta = new int[len];
    int *tb = new int[len];
    
    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp sections
        {
            #pragma omp section
            {
                std::fill(ta, ta + len, 0);
                for (int i = 0; i < n; i++) {
                    ta[i] = a[i];
                }
            }
            #pragma omp section
            {
                std::fill(tb, tb + len, 0);
                for (int i = 0; i < n; i++) {
                    tb[i] = b[i];
                }
            }
        }
    }
    
    // NTT变换
    openmp_ntt(ta, len, 1, p, num_threads);
    openmp_ntt(tb, len, 1, p, num_threads);
    
    // 点值乘法
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < len; i++) {
        ta[i] = 1LL * ta[i] * tb[i] % p;
    }
    
    // INTT变换
    openmp_ntt(ta, len, -1, p, num_threads);
    
    // 结果赋值给c
    for (int i = 0; i < 2 * n - 1; i++) {
        c[i] = ta[i];
    }
    
    delete[] ta;
    delete[] tb;
}
//============================ 四分NTT实现 ============================
// 判断是否为4的幂
bool is_power_of_4(int n) {
    if (n <= 0) return false;
    if ((n & (n - 1)) != 0) return false; // 先确保是2的幂
    // A power of 2, n = 2^k. n is power of 4 if k is even.
    // Check if (n & 0xAAAAAAAA) == 0.
    // 0xAAAAAAAA is 10101010... in binary. If n has a bit set here, k is odd (e.g. 2, 8, 32...).
    // So, if it's a power of 2, and (n & 0xAAAAAAAA) == 0, it means the 1 is at an even position (0, 2, 4...).
    // Example: 1 (0001), (1 & 0xAA...) = 0. Correct.
    // Example: 4 (0100), (4 & 0xAA...) = 0. Correct.
    // Example: 16 (00010000), (16 & 0xAA...) = 0. Correct.
    // Example: 2 (0010), (2 & 0xAA...) = 2. Not power of 4. Correct.
    return (n & 0xAAAAAAAA) == 0;
}

// 获取大于等于n的最小的4的幂
int get_min_power_of_4(int n) {
    int res = 1;
    while (res < n) res <<= 2;
    return res;
}

// 四分NTT的迭代实现(串行)
// 修改四分NTT蝴蝶操作中的计算
// 修改四分NTT蝴蝶操作中的计算
// 改进的四分NTT函数
// 改进的四分NTT函数
void four_div_ntt(int *a, int len, int type, int p) {
    // Ensure len is a power of 4, handled by four_div_ntt_multiply
    bit_reverse(a, len);
    
    for (int m = 4; m <= len; m <<= 2) {
        int m4 = m / 4;
        long long W_m_root = mod_pow(3, (p - 1) / m, p); // W_m (m-th principal root)
        if (type == -1) {
            W_m_root = mod_inverse(W_m_root, p);
        }

        for (int j = 0; j < len; j += m) {
            long long w_pow_k = 1; // This will be W_m^k, starts at W_m^0
            for (int k = 0; k < m4; k++) { 
                long long w_pow_2k = (w_pow_k * w_pow_k) % p;
                long long w_pow_3k = (w_pow_2k * w_pow_k) % p;

                int idx0 = j + k;
                int idx1 = j + k + m4;
                int idx2 = j + k + 2 * m4;
                int idx3 = j + k + 3 * m4;

                long long x0 = a[idx0];
                long long x1 = a[idx1];
                long long x2 = a[idx2];
                long long x3 = a[idx3];

                long long t0 = x0;
                long long t1 = (x1 * w_pow_k) % p;  
                long long t2 = (x2 * w_pow_2k) % p; 
                long long t3 = (x3 * w_pow_3k) % p; 
                
                long long X0 = (t0 + t1 + t2 + t3);
                long long X1 = (t0 - t1 + t2 - t3); 
                long long X2 = (t0 + t1 - t2 - t3); 
                long long X3 = (t0 - t1 - t2 + t3); 

                a[idx0] = (X0 % p + p) % p;
                a[idx1] = (X1 % p + p) % p; 
                a[idx2] = (X2 % p + p) % p;
                a[idx3] = (X3 % p + p) % p;
                
                w_pow_k = (w_pow_k * W_m_root) % p; 
            }
        }
    }
    
    if (type == -1) {
        long long inv_len = mod_inverse(len, p);
        for (int i = 0; i < len; i++) {
            a[i] = (1LL * a[i] * inv_len) % p;
            if (a[i] < 0) a[i] += p; 
        }
    }
}
// 四分NTT多项式乘法
void four_div_ntt_multiply(int *a_orig, int *b_orig, int *c_final, int n, int p) {
    int len = get_min_power_of_4(2 * n - 1); 
    
    int *ta = new int[len];
    int *tb = new int[len];
    std::fill(ta, ta + len, 0);
    std::fill(tb, tb + len, 0);
    
    std::copy(a_orig, a_orig + n, ta);
    std::copy(b_orig, b_orig + n, tb);
    
    four_div_ntt(ta, len, 1, p);
    four_div_ntt(tb, len, 1, p);
    
    for (int i = 0; i < len; i++) {
        ta[i] = (1LL * ta[i] * tb[i]) % p;
    }
    
    four_div_ntt(ta, len, -1, p);
    
    std::copy(ta, ta + (2 * n - 1), c_final);
    
    delete[] ta;
    delete[] tb;
}


// pthread四分NTT线程函数
void* pthread_four_div_ntt_worker(void *arg) {
    ThreadArgs *args = (ThreadArgs*)arg;
    int *a = args->a;
    int len = args->len;
    int m = args->h; // h is repurposed as m (current DFT size)
    int p_val = args->p; // Renamed p to p_val
    pthread_barrier_t *barrier = args->barrier;
    // Using start_index and work_amount for k-loop division
    int k_start_loop = args->start_index; // Renamed k_start to k_start_loop
    int k_count = args->work_amount;
    int k_end_loop = k_start_loop + k_count; // Renamed k_end to k_end_loop
    
    int m4 = m / 4;
    if (m4 == 0) { // Should not happen if m starts at 4
        if (barrier) pthread_barrier_wait(barrier);
        return NULL;
    }

    long long wm_root_val = mod_pow(3, (p_val - 1) / m, p_val); 
    // args->gn is repurposed to pass 'type'
    if (args->gn == -1) wm_root_val = mod_inverse(wm_root_val, p_val); 
    
    // Ensure k_end_loop does not exceed m4
    k_end_loop = std::min(k_end_loop, m4);

    if (k_start_loop >= k_end_loop) { // No work for this thread in k-loop
        if (barrier) pthread_barrier_wait(barrier);
        return NULL;
    }
    
    // Calculate initial w_pow_k for this thread's k_start_loop
    long long current_w_pow_k = mod_pow(wm_root_val, k_start_loop, p_val);

    for (int k_iter = k_start_loop; k_iter < k_end_loop; k_iter++) { // Renamed k to k_iter
        long long w_pow_2k = (current_w_pow_k * current_w_pow_k) % p_val;
        long long w_pow_3k = (w_pow_2k * current_w_pow_k) % p_val;
        
        for (int j = 0; j < len; j += m) {
            // No need for j + m > len check, outer loop ensures j < len
            // No need for j + m4*3 >= len check if len is multiple of m (which it is)
            
            int idx0 = j + k_iter;
            int idx1 = j + k_iter + m4;
            int idx2 = j + k_iter + 2 * m4;
            int idx3 = j + k_iter + 3 * m4;

            // Boundary checks for safety, though theoretically not needed if len is a power of 4 and m divides len
            if (idx3 >= len) continue;


            long long x0 = a[idx0];
            long long x1 = a[idx1];
            long long x2 = a[idx2];
            long long x3 = a[idx3];

            long long t0 = x0;
            long long t1 = (x1 * current_w_pow_k) % p_val;  
            long long t2 = (x2 * w_pow_2k) % p_val; 
            long long t3 = (x3 * w_pow_3k) % p_val; 
            
            long long X0 = (t0 + t1 + t2 + t3);
            long long X1 = (t0 - t1 + t2 - t3); 
            long long X2 = (t0 + t1 - t2 - t3); 
            long long X3 = (t0 - t1 - t2 + t3); 

            a[idx0] = (X0 % p_val + p_val) % p_val;
            a[idx1] = (X1 % p_val + p_val) % p_val; 
            a[idx2] = (X2 % p_val + p_val) % p_val;
            a[idx3] = (X3 % p_val + p_val) % p_val;
        }
        current_w_pow_k = (current_w_pow_k * wm_root_val) % p_val; // Update for next k_iter
    }
    
    if (barrier) {
        pthread_barrier_wait(barrier);
    }
    return NULL;
}
// pthread四分NTT实现
void pthread_four_div_ntt(int *a_ptr, int len, int type, int p_val, int num_threads_req) { // Renamed a, p, num_threads
    if (num_threads_req <= 0) num_threads_req = 1;
    int actual_num_threads = std::min(num_threads_req, len / 4); // Heuristic for four_div
    if (actual_num_threads < 1) actual_num_threads = 1;

    // Check if len is a power of 2, required for bit_reverse and NTT structure
    if ((len & (len - 1)) != 0 || len == 0) {
        std::cout << "错误: pthread_four_div_ntt 输入长度不是2的幂次或为0，使用普通pthread_ntt" << std::endl;
        pthread_ntt(a_ptr, len, type, p_val, actual_num_threads);
        return;
    }
    // Four_div_ntt typically expects len to be a power of 4 for optimal structure.
    // If not, it might need a preliminary radix-2 stage or careful handling.
    // The serial four_div_ntt assumes len is already a power of 4.
    // For simplicity here, we'll proceed, but this could be a source of issues if len is 2, 8, 32 etc.
    // A robust version would call get_min_power_of_4 for len in the multiply function.

    bit_reverse(a_ptr, len);
    
    pthread_t *threads_arr = new pthread_t[actual_num_threads]; // Renamed threads
    ThreadArgs *args_arr = new ThreadArgs[actual_num_threads];   // Renamed args
    pthread_barrier_t barrier_obj; // Renamed barrier
    
    bool barrier_active = false; // Renamed barrier_initialized
    if (actual_num_threads > 1) {
        if (pthread_barrier_init(&barrier_obj, NULL, actual_num_threads) != 0) {
            std::cerr << "pthread_four_div_ntt: 创建barrier失败!" << std::endl;
            delete[] threads_arr;
            delete[] args_arr;
            four_div_ntt(a_ptr, len, type, p_val); // Fallback to serial four_div
            return;
        }
        barrier_active = true;
    }
    
    // The initial radix-2 stage if len is 2 * (power of 4) is complex to parallelize correctly
    // with the subsequent radix-4 stages. The original serial four_div_ntt does not have this.
    // We will assume len is a power of 4 for the parallel radix-4 stages.
    // If four_div_ntt_multiply ensures len is power of 4, this is fine.

    for (int m_loop = 4; m_loop <= len; m_loop <<= 2) { // Renamed m to m_loop
        if (actual_num_threads > 1) {
            int m4 = m_loop / 4;
            if (m4 == 0) continue; // Should not happen
            int k_total_iterations = m4;
            int k_per_thread = (k_total_iterations + actual_num_threads - 1) / actual_num_threads;

            for (int t = 0; t < actual_num_threads; t++) {
                args_arr[t].a = a_ptr;
                args_arr[t].len = len;
                args_arr[t].h = m_loop; // h is m (current DFT size)
                args_arr[t].p = p_val;
                args_arr[t].gn = type; // Pass type via gn
                args_arr[t].barrier = &barrier_obj;
                args_arr[t].thread_id = t; // Keep original meaning
                args_arr[t].num_threads = actual_num_threads; // Keep original meaning
                args_arr[t].start_index = t * k_per_thread;
                args_arr[t].work_amount = k_per_thread; // Each thread attempts to do this many k iterations
                                                       // Worker will cap at m4

                if (pthread_create(&threads_arr[t], NULL, pthread_four_div_ntt_worker, &args_arr[t]) != 0) {
                    std::cerr << "pthread_four_div_ntt: 创建线程 " << t << " 失败!" << std::endl;
                    for (int i_cancel = 0; i_cancel < t; i_cancel++) { // Renamed i to i_cancel
                        // pthread_cancel is generally not safe with mutexes/barriers
                        pthread_join(threads_arr[i_cancel], NULL); // Wait for already started threads
                    }
                    if (barrier_active) pthread_barrier_destroy(&barrier_obj);
                    delete[] threads_arr;
                    delete[] args_arr;
                    four_div_ntt(a_ptr, len, type, p_val); // Fallback
                    return;
                }
            }
            
            for (int t = 0; t < actual_num_threads; t++) {
                pthread_join(threads_arr[t], NULL);
            }
        } else { // Serial execution for this m_loop level
            int m4 = m_loop / 4;
            long long W_m_root_serial = mod_pow(3, (p_val - 1) / m_loop, p_val);
            if (type == -1) W_m_root_serial = mod_inverse(W_m_root_serial, p_val);

            for (int j_serial = 0; j_serial < len; j_serial += m_loop) { // Renamed j to j_serial
                long long w_pow_k_serial = 1;
                for (int k_serial = 0; k_serial < m4; k_serial++) { // Renamed k to k_serial
                    long long w_pow_2k_s = (w_pow_k_serial * w_pow_k_serial) % p_val;
                    long long w_pow_3k_s = (w_pow_2k_s * w_pow_k_serial) % p_val;

                    int idx0 = j_serial + k_serial;
                    int idx1 = j_serial + k_serial + m4;
                    int idx2 = j_serial + k_serial + 2 * m4;
                    int idx3 = j_serial + k_serial + 3 * m4;

                    long long x0 = a_ptr[idx0]; long long x1 = a_ptr[idx1];
                    long long x2 = a_ptr[idx2]; long long x3 = a_ptr[idx3];

                    long long t0 = x0;
                    long long t1 = (x1 * w_pow_k_serial) % p_val;
                    long long t2 = (x2 * w_pow_2k_s) % p_val;
                    long long t3 = (x3 * w_pow_3k_s) % p_val;

                    long long X0 = (t0 + t1 + t2 + t3);
                    long long X1 = (t0 - t1 + t2 - t3);
                    long long X2 = (t0 + t1 - t2 - t3);
                    long long X3 = (t0 - t1 - t2 + t3);

                    a_ptr[idx0] = (X0 % p_val + p_val) % p_val;
                    a_ptr[idx1] = (X1 % p_val + p_val) % p_val;
                    a_ptr[idx2] = (X2 % p_val + p_val) % p_val;
                    a_ptr[idx3] = (X3 % p_val + p_val) % p_val;
                    
                    w_pow_k_serial = (w_pow_k_serial * W_m_root_serial) % p_val;
                }
            }
        }
    }
      
    if (type == -1) {
        long long inv_len_val = mod_inverse(len, p_val); // Renamed inv_len
        if (actual_num_threads > 1) {
            // Lambda for final scaling (reusing ThreadArgs and pthread_t arrays)
            auto scaling_worker_lambda_fd = [](void *arg_lambda_fd) -> void* { // Renamed arg and lambda
                ThreadArgs *args_lambda_fd = (ThreadArgs*)arg_lambda_fd;
                int *arr_fd = args_lambda_fd->a; // Renamed arr
                // args_lambda_fd->start_index and args_lambda_fd->work_amount are used here
                int start_idx_fd = args_lambda_fd->start_index; // Renamed start_idx
                int end_idx_fd = args_lambda_fd->start_index + args_lambda_fd->work_amount; // Renamed end_idx
                int p_val_fd = args_lambda_fd->p; // Renamed p_val
                long long inv_val_fd = args_lambda_fd->gn; // Renamed inv_val (gn is repurposed for inv_len)
                
                for (int i_scale = start_idx_fd; i_scale < end_idx_fd; i_scale++) { // Renamed i to i_scale
                    arr_fd[i_scale] = (int)(1LL * arr_fd[i_scale] * inv_val_fd % p_val_fd);
                    if (arr_fd[i_scale] < 0) arr_fd[i_scale] += p_val_fd;
                }
                return NULL;
            };
            
            int chunk_size_scale = (len + actual_num_threads - 1) / actual_num_threads;
            for (int t = 0; t < actual_num_threads; t++) {
                args_arr[t].a = a_ptr;
                args_arr[t].len = len; 
                args_arr[t].start_index = t * chunk_size_scale;
                args_arr[t].work_amount = std::min(chunk_size_scale, len - args_arr[t].start_index);
                args_arr[t].h = 0; // Not used
                args_arr[t].p = p_val;
                args_arr[t].gn = inv_len_val; // Pass inv_len
                args_arr[t].barrier = NULL; // No barrier for scaling
                args_arr[t].thread_id = t; // Keep original meaning
                args_arr[t].num_threads = actual_num_threads; // Keep original meaning
                
                if (pthread_create(&threads_arr[t], NULL, scaling_worker_lambda_fd, &args_arr[t]) != 0) {
                     std::cerr << "pthread_four_div_ntt: 创建缩放线程 " << t << " 失败!" << std::endl;
                     for(int i_cancel=0; i_cancel<t; ++i_cancel) pthread_join(threads_arr[i_cancel], NULL);
                     if(barrier_active) pthread_barrier_destroy(&barrier_obj);
                     delete[] threads_arr;
                     delete[] args_arr;
                     // Fallback scaling
                     for(int i_fb_scale=0; i_fb_scale<len; ++i_fb_scale) 
                        a_ptr[i_fb_scale] = (1LL*a_ptr[i_fb_scale]*inv_len_val)%p_val;
                     return;
                }
            }
            
            for (int t = 0; t < actual_num_threads; t++) {
                pthread_join(threads_arr[t], NULL);
            }
        } else { // Serial scaling
            for (int i_s = 0; i_s < len; i_s++) { // Renamed i to i_s
                a_ptr[i_s] = (int)(1LL * a_ptr[i_s] * inv_len_val % p_val);
                if (a_ptr[i_s] < 0) a_ptr[i_s] += p_val;
            }
        }
    }
    
    if (barrier_active) {
        pthread_barrier_destroy(&barrier_obj);
    }
    delete[] threads_arr;
    delete[] args_arr;
}

// pthread四分NTT多项式乘法
// pthread四分NTT多项式乘法的修复版本
void pthread_four_div_ntt_multiply(int *a_orig, int *b_orig, int *c_final, int n_val, int p_val, int num_threads_req) { // Renamed n,p,num_threads
    if (num_threads_req <= 0) num_threads_req = 1;
    int len = get_min_power_of_4(2 * n_val - 1); // Ensure len is power of 4 for four_div_ntt

    int *ta = new int[len];
    int *tb = new int[len];
    
    // Parallel initialization using std::thread for simplicity here
    // Could also use pthreads if preferred for consistency
    if (num_threads_req > 1) {
        std::vector<std::thread> init_threads_fd; // Renamed init_threads
        init_threads_fd.emplace_back([&]() {
            std::fill(ta, ta + len, 0);
            std::copy(a_orig, a_orig + n_val, ta);
        });
        init_threads_fd.emplace_back([&]() {
            std::fill(tb, tb + len, 0);
            std::copy(b_orig, b_orig + n_val, tb);
        });
        for(auto& th : init_threads_fd) th.join();
    } else {
        std::fill(ta, ta + len, 0);
        std::copy(a_orig, a_orig + n_val, ta);
        std::fill(tb, tb + len, 0);
        std::copy(b_orig, b_orig + n_val, tb);
    }

    pthread_four_div_ntt(ta, len, 1, p_val, num_threads_req);
    pthread_four_div_ntt(tb, len, 1, p_val, num_threads_req);

    // Pointwise multiplication (parallelized using std::thread for simplicity)
    if (num_threads_req > 1) {
        std::vector<std::thread> pointwise_threads_fd; // Renamed pointwise_threads
        int actual_num_threads_pm_fd = std::min(num_threads_req, len); // Renamed num_actual_threads_pm
        if(actual_num_threads_pm_fd < 1) actual_num_threads_pm_fd = 1;

        for (int t = 0; t < actual_num_threads_pm_fd; t++) {
            pointwise_threads_fd.emplace_back([t, actual_num_threads_pm_fd, len, ta, tb, p_val]() {
                int chunk_size = (len + actual_num_threads_pm_fd - 1) / actual_num_threads_pm_fd;
                int start_idx = t * chunk_size;
                int end_idx = std::min(start_idx + chunk_size, len);
                for (int i_pm = start_idx; i_pm < end_idx; i_pm++) { // Renamed i to i_pm
                    ta[i_pm] = (1LL * ta[i_pm] * tb[i_pm]) % p_val;
                }
            });
        }
        for (auto& th : pointwise_threads_fd) {
            th.join();
        }
    } else {
        for (int i_pm_s = 0; i_pm_s < len; i_pm_s++) { // Renamed i to i_pm_s
            ta[i_pm_s] = (1LL * ta[i_pm_s] * tb[i_pm_s]) % p_val;
        }
    }
    
    pthread_four_div_ntt(ta, len, -1, p_val, num_threads_req);

    std::copy(ta, ta + (2 * n_val - 1), c_final);

    delete[] ta;
    delete[] tb;
}

// OpenMP四分NTT实现
void openmp_four_div_ntt(int *a_ptr, int len, int type, int p_val, int num_threads_req) { // Renamed a,p,num_threads
    if (num_threads_req <=0) num_threads_req = 1;
    // Ensure len is power of 4 for this specific four_div_ntt logic
    // The calling multiply function should ensure this.
    if (!is_power_of_4(len) || (len & (len-1)) != 0 ) { // Stricter check
        std::cout << "错误: openmp_four_div_ntt 输入长度不是4的幂次，使用普通OpenMP NTT" << std::endl;
        openmp_ntt(a_ptr, len, type, p_val, num_threads_req);
        return;
    }
    
    bit_reverse(a_ptr, len);
    
    // No initial radix-2 stage in this version, assumes len is power of 4.
    
    for (int m_loop = 4; m_loop <= len; m_loop <<= 2) { // Renamed m to m_loop
        int m4 = m_loop / 4;
        long long wm_root_omp = mod_pow(3, (p_val - 1) / m_loop, p_val); // Renamed wm
        if (type == -1) wm_root_omp = mod_inverse(wm_root_omp, p_val);
        
        #pragma omp parallel for num_threads(num_threads_req) schedule(dynamic)
        for (int k_iter = 0; k_iter < m4; k_iter++) { // Renamed k to k_iter
            long long w_k_omp = mod_pow(wm_root_omp, k_iter, p_val); // Renamed w
            long long w2_k_omp = (w_k_omp * w_k_omp) % p_val; // Renamed w2
            long long w3_k_omp = (w2_k_omp * w_k_omp) % p_val; // Renamed w3
            
            for (int j_iter = 0; j_iter < len; j_iter += m_loop) { // Renamed j to j_iter
                int idx0 = j_iter + k_iter;
                int idx1 = j_iter + k_iter + m4;
                int idx2 = j_iter + k_iter + 2 * m4;
                int idx3 = j_iter + k_iter + 3 * m4;

                // Boundary check though theoretically not needed if len is power of 4 and m_loop divides len
                if (idx3 >= len) continue;

                long long x0 = a_ptr[idx0]; long long x1 = a_ptr[idx1];
                long long x2 = a_ptr[idx2]; long long x3 = a_ptr[idx3];
                
                long long t0 = x0;
                long long t1 = (x1 * w_k_omp) % p_val;
                long long t2 = (x2 * w2_k_omp) % p_val;
                long long t3 = (x3 * w3_k_omp) % p_val;
                
                long long X0 = (t0 + t1 + t2 + t3);
                long long X1 = (t0 - t1 + t2 - t3);
                long long X2 = (t0 + t1 - t2 - t3);
                long long X3 = (t0 - t1 - t2 + t3);
                
                a_ptr[idx0] = (X0 % p_val + p_val) % p_val;
                a_ptr[idx1] = (X1 % p_val + p_val) % p_val;
                a_ptr[idx2] = (X2 % p_val + p_val) % p_val;
                a_ptr[idx3] = (X3 % p_val + p_val) % p_val;
            }
        }
    }
    
    if (type == -1) {
        long long inv_len_val = mod_inverse(len, p_val); // Renamed inv_len
        #pragma omp parallel for num_threads(num_threads_req)
        for (int i_s = 0; i_s < len; i_s++) { // Renamed i to i_s
            a_ptr[i_s] = (int)(1LL * a_ptr[i_s] * inv_len_val % p_val);
            if (a_ptr[i_s] < 0) a_ptr[i_s] += p_val;
        }
    }
}
// OpenMP四分NTT多项式乘法
// OpenMP四分NTT多项式乘法
void openmp_four_div_ntt_multiply(int *a_orig, int *b_orig, int *c_final, int n_val, int p_val, int num_threads_req) { // Renamed n,p,num_threads
    if (num_threads_req <=0) num_threads_req = 1;
    int len = get_min_power_of_4(2 * n_val - 1);
    
    int *ta = new int[len];
    int *tb = new int[len];
    
    #pragma omp parallel sections num_threads(std::min(num_threads_req,2))
    {
        #pragma omp section
        {
            std::fill(ta, ta + len, 0);
            std::copy(a_orig, a_orig + n_val, ta);
        }
        #pragma omp section
        {
            std::fill(tb, tb + len, 0);
            std::copy(b_orig, b_orig + n_val, tb);
        }
    }
    
    openmp_four_div_ntt(ta, len, 1, p_val, num_threads_req);
    openmp_four_div_ntt(tb, len, 1, p_val, num_threads_req);
    
    #pragma omp parallel for num_threads(num_threads_req)
    for (int i_pm = 0; i_pm < len; i_pm++) { // Renamed i to i_pm
        ta[i_pm] = (1LL * ta[i_pm] * tb[i_pm]) % p_val;
    }
    
    openmp_four_div_ntt(ta, len, -1, p_val, num_threads_req);
    
    std::copy(ta, ta + (2 * n_val - 1), c_final);
    
    delete[] ta;
    delete[] tb;
}
//============================ CRT多线程实现 ============================
// 预定义的互质模数
const int MOD_COUNT = 3;
const int MODULI[MOD_COUNT] = {469762049, 998244353, 167772161};  // 三个互质模数

// CRT合并多个模数下的结果
struct CrtResult {
    int *a;
    int len;
};

// CRT线程函数参数
struct CrtThreadArgs {
    int *a;
    int *b;
    int n;
    int p;
    int *c;
};

// 预计算CRT常量结构体
struct CrtConstants {
    cpp_int M;
    cpp_int M_divs[MOD_COUNT]; // 大小应与 MOD_COUNT 一致
    cpp_int M_invs[MOD_COUNT]; // 大小应与 MOD_COUNT 一致
    
    CrtConstants() {
        M = 1;
        for (int i = 0; i < MOD_COUNT; i++) {
            M *= MODULI[i];
        }
        for (int i = 0; i < MOD_COUNT; i++) {
            M_divs[i] = M / MODULI[i];
            // Ensure M_divs[i] is taken modulo MODULI[i] before inverse if M_divs[i] can be > MODULI[i]
            cpp_int M_div_mod_mi = M_divs[i] % cpp_int(MODULI[i]);
            if (M_div_mod_mi < 0) M_div_mod_mi += MODULI[i]; // Ensure positive before inverse
            M_invs[i] = mod_inverse_big(M_div_mod_mi, cpp_int(MODULI[i]));
        }
    }
};
const CrtConstants& get_crt_constants() {
    static CrtConstants instance;
    return instance;
}


// pthread CRT线程函数
// pthread CRT线程函数 (可以进一步优化内存操作)
void* pthread_crt_ntt_worker(void *arg) {
    CrtThreadArgs *args_ptr = (CrtThreadArgs*)arg; // Renamed args to args_ptr
    int *a_input = args_ptr->a; // Renamed a to a_input
    int *b_input = args_ptr->b; // Renamed b to b_input
    int n_val = args_ptr->n;    // Renamed n to n_val
    int p_mod = args_ptr->p;    // Renamed p to p_mod
    int *c_output = args_ptr->c; // Renamed c to c_output
    
    int len = get_min_power_of_2(2 * n_val - 1);
    
    // 一次性分配内存并填充
    // 使用 std::vector 可以自动管理内存，但 new/delete 也可以
    int *ta = new int[len];
    int *tb = new int[len];
    
    // 初始化为0不是必须的，因为后续会复制或计算每个元素
    // std::fill(ta, ta + len, 0); 
    // std::fill(tb, tb + len, 0);

    for (int i = 0; i < n_val; i++) {
        long long val_a = a_input[i];
        long long val_b = b_input[i];
        // Ensure values are within [0, p_mod-1] before NTT
        ta[i] = static_cast<int>((val_a % p_mod + p_mod) % p_mod);
        tb[i] = static_cast<int>((val_b % p_mod + p_mod) % p_mod);
    }
    // Pad the rest with zeros
    std::fill(ta + n_val, ta + len, 0);
    std::fill(tb + n_val, tb + len, 0);
    
    // 使用最稳定的NTT实现 (您当前的 ntt 函数)
    ntt(ta, len, 1, p_mod);
    ntt(tb, len, 1, p_mod);
    
    // 高效点值乘法
    for (int i = 0; i < len; i++) {
        ta[i] = static_cast<int>((1LL * ta[i] * tb[i]) % p_mod);
        // No need for if (ta[i] < 0) ta[i] += p_mod; if p_mod is positive and % gives result in [-p_mod+1, p_mod-1]
        // However, to be absolutely safe for all C++ % behaviors with negative numbers:
        if (ta[i] < 0) ta[i] += p_mod;
    }
    
    ntt(ta, len, -1, p_mod);
    
    // 直接复制结果到输出数组 c_output
    // memcpy is generally faster for plain old data types
    memcpy(c_output, ta, (2 * n_val - 1) * sizeof(int));
    
    delete[] ta;
    delete[] tb;
    return NULL;
}

// 修改CRT合并函数，使用大整数处理
// 优化的crt_combine函数
void crt_combine(int *results[], int n_orig, int p_target, int *c_final) {
    const CrtConstants& consts = get_crt_constants(); // 获取预计算的常量
    
    int result_len = 2 * n_orig - 1;

    // 如果 result_len 很大，可以考虑并行化这个外层循环
    // #pragma omp parallel for // 如果启用，需要考虑线程安全和数据竞争
    for (int i = 0; i < result_len; i++) {
        cpp_int current_sum_mod_M = 0;

        for (int j = 0; j < MOD_COUNT; j++) {
            cpp_int val_j = results[j][i]; 
            // val_j 已经是模 MODULI[j] 的结果，但可能为负，确保为正
            val_j = (val_j % MODULI[j] + MODULI[j]) % MODULI[j];

            // CRT项计算: term_j = val_j * M_invs[j] * M_divs[j]
            // (val_j * consts.M_invs[j]) % MODULI[j] * consts.M_divs[j] - this is not standard CRT formula structure
            // Standard: sum ( a_i * (M/m_i) * inv((M/m_i), m_i) ) mod M
            
            cpp_int term = val_j;
            term *= consts.M_invs[j];       // term = val_j * inv((M/m_j) mod m_j, m_j)
            term %= consts.M;               // Reduce intermediate product if M_invs[j] is large
                                            // This modulo is not strictly part of the standard formula here,
                                            // but can help keep numbers smaller if M_invs[j] is large.
                                            // More accurately: term = val_j * M_invs[j] (this product is mod MODULI[j] implicitly by M_invs)
                                            // then multiply by M_divs[j]
            
            // Corrected term calculation:
            // term = val_j * consts.M_invs[j] (this is val_j * inv(M/m_j mod m_j, m_j) )
            // This product should be taken modulo MODULI[j] if we want to be precise about steps,
            // but since M_invs[j] is already inv( (M/m_j) mod m_j, m_j),
            // (val_j * consts.M_invs[j]) % MODULI[j] effectively gives (val_j * ( (M/m_j) mod m_j )^-1) mod m_j
            // Then this is multiplied by M_divs[j] = M/m_j.
            // The sum is: sum_j ( val_j * M_invs[j] * M_divs[j] )
            // Let's use the direct formula: (a_i * (M/m_i) * inv((M/m_i) mod m_i, m_i))

            cpp_int single_term = val_j;                  // a_i
            single_term *= consts.M_divs[j];          // a_i * (M/m_i)
            // single_term %= consts.M; // Optional: reduce before multiplying by inverse, but M_invs is inv of (M_divs[j] % MODULI[j])
            single_term *= consts.M_invs[j];          // a_i * (M/m_i) * inv((M/m_i) mod m_i, m_i)
                                                    // This product is implicitly modulo M because M_invs[j] is inv mod MODULI[j]
                                                    // and M_divs[j] contains all other moduli.
                                                    // The sum of these terms will be taken modulo M.
            current_sum_mod_M += single_term;
            // current_sum_mod_M %= consts.M; // Defer final modulo M to after the loop for sum
        }

        current_sum_mod_M %= consts.M; // Final sum modulo M
        if (current_sum_mod_M < 0) {
            current_sum_mod_M += consts.M;
        }
        
        // 最后，将结果对目标模数 p_target 取模
        cpp_int final_value_mod_p = current_sum_mod_M % cpp_int(p_target); // Ensure p_target is cpp_int for %
        if (final_value_mod_p < 0) {
            final_value_mod_p += p_target;
        }
        
        c_final[i] = static_cast<int>(final_value_mod_p);
    }
}
// pthread CRT NTT多项式乘法
// 优化的pthread CRT NTT多项式乘法
void pthread_crt_ntt_multiply(int *a, int *b, int *c, int n, int p_target) { // Renamed p to p_target
    pthread_t threads[MOD_COUNT];
    CrtThreadArgs args[MOD_COUNT];
    int *results[MOD_COUNT]; // Array of pointers to results for each modulus
    
    // std::cout << "   使用优化的" << MOD_COUNT << "线程CRT计算" << std::endl; // More generic message
    
    // 一次性分配所有内存
    int result_len_for_alloc = get_min_power_of_2(2 * n - 1); // Allocate enough for NTT length

    for (int i = 0; i < MOD_COUNT; i++) {
        // Allocate memory for results[i] to hold the output of ntt which is of 'len'
        // but crt_combine will only use up to 2*n-1 elements.
        // pthread_crt_ntt_worker writes 2*n-1 elements.
        results[i] = new int[2 * n - 1]; // Only need 2*n-1 for final result storage
        args[i] = {a, b, n, MODULI[i], results[i]}; // Pass correct result buffer
        if (pthread_create(&threads[i], NULL, pthread_crt_ntt_worker, &args[i]) != 0) {
            std::cerr << "pthread_crt_ntt_multiply: Failed to create thread " << i << std::endl;
            // Handle error: cleanup already created threads and allocated memory
            for (int k=0; k<i; ++k) pthread_join(threads[k], NULL);
            for (int k=0; k<=i; ++k) delete[] results[k]; // results[i] was allocated
            // Potentially re-throw or exit
            return;
        }
    }
    
    // 等待所有线程完成
    for (int i = 0; i < MOD_COUNT; i++) {
        pthread_join(threads[i], NULL);
    }
    
    // 使用优化的CRT合并
    crt_combine(results, n, p_target, c);
    
    // 释放内存
    for (int i = 0; i < MOD_COUNT; i++) {
        delete[] results[i];
    }
}

// OpenMP CRT NTT多项式乘法
void openmp_crt_ntt_multiply(int *a, int *b, int *c, int n, int p_target) { // Renamed p to p_target
    int *results[MOD_COUNT];
    int output_len = 2 * n - 1; // Length of the final polynomial product

    for (int i = 0; i < MOD_COUNT; i++) {
        results[i] = new int[output_len]; // Allocate for final product length
    }
    
    // std::cout << "   使用" << MOD_COUNT << "个固定线程进行CRT计算" << std::endl;
    
    #pragma omp parallel for num_threads(MOD_COUNT)
    for (int i = 0; i < MOD_COUNT; i++) {
        int len_ntt = get_min_power_of_2(output_len); // NTT length
        int mod_current = MODULI[i]; // Renamed mod to mod_current
        
        // Thread-local temporary arrays for NTT
        int *ta = new int[len_ntt];
        int *tb = new int[len_ntt];
        // No need to fill with 0 if copy/calc covers all used elements
        // std::fill(ta, ta + len_ntt, 0);
        // std::fill(tb, tb + len_ntt, 0);
        
        for (int j = 0; j < n; j++) {
            long long val_a = a[j];
            long long val_b = b[j];
            ta[j] = static_cast<int>((val_a % mod_current + mod_current) % mod_current);
            tb[j] = static_cast<int>((val_b % mod_current + mod_current) % mod_current);
        }
        std::fill(ta + n, ta + len_ntt, 0);
        std::fill(tb + n, tb + len_ntt, 0);

        ntt(ta, len_ntt, 1, mod_current);
        ntt(tb, len_ntt, 1, mod_current);
        
        for (int j = 0; j < len_ntt; j++) {
            ta[j] = static_cast<int>((1LL * ta[j] * tb[j]) % mod_current);
            if (ta[j] < 0) ta[j] += mod_current;
        }
        
        ntt(ta, len_ntt, -1, mod_current);
        
        // Copy only the necessary part to results[i]
        memcpy(results[i], ta, output_len * sizeof(int));
        
        delete[] ta;
        delete[] tb;
    }
    
    // 使用CRT合并结果
    crt_combine(results, n, p_target, c);
    
    for (int i = 0; i < MOD_COUNT; i++) {
        delete[] results[i];
    }
}
//============================ 测试函数 ============================
// 测试各种实现的性能
// 测试各种实现的性能
void test_performance(int *a, int *b, int *c, int n, int p, int current_input_id) { // 添加 current_input_id 参数
    int num_threads_list[] = {1, 2, 4, 8};
    int num_threads_count = sizeof(num_threads_list) / sizeof(int);
    
    std::chrono::high_resolution_clock::time_point start, end;
    long long duration;

    std::cout << "\n===== NTT性能测试 (n=" << n << ", p=" << p << ", input_id=" << current_input_id << ") =====" << std::endl;
    
    // 0. 朴素多项式乘法
    std::cout << "0. 朴素多项式乘法:" << std::endl;
    memset(c, 0, (2 * n - 1) * sizeof(int)); // 注意 memset 的大小参数
    start = std::chrono::high_resolution_clock::now();
    poly_multiply(a, b, c, n, p);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "   耗时: " << duration << " us" << std::endl;
    fCheck(c, n, current_input_id); // 使用 current_input_id
    
    // 1. 串行NTT
    std::cout << "1. 串行NTT:" << std::endl;
    memset(c, 0, (2 * n - 1) * sizeof(int));
    start = std::chrono::high_resolution_clock::now();
    ntt_multiply(a, b, c, n, p);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "   耗时: " << duration << " us" << std::endl;
    fCheck(c, n, current_input_id); // 使用 current_input_id

    // 2. pthread朴素NTT
    std::cout << "2. pthread朴素NTT:" << std::endl;
    for (int i = 0; i < num_threads_count; i++) {
        int num_threads = num_threads_list[i];
        memset(c, 0, (2 * n - 1) * sizeof(int));
        start = std::chrono::high_resolution_clock::now();
        pthread_ntt_multiply(a, b, c, n, p, num_threads);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "   " << num_threads << "线程耗时: " << duration << " us" << std::endl;
        fCheck(c, n, current_input_id); // 使用 current_input_id
    }
    
    // 3. OpenMP朴素NTT
    std::cout << "3. OpenMP朴素NTT:" << std::endl;
    for (int i = 0; i < num_threads_count; i++) {
        int num_threads = num_threads_list[i];
        memset(c, 0, (2 * n - 1) * sizeof(int));
        start = std::chrono::high_resolution_clock::now();
        openmp_ntt_multiply(a, b, c, n, p, num_threads);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "   " << num_threads << "线程耗时: " << duration << " us" << std::endl;
        fCheck(c, n, current_input_id); // 使用 current_input_id
    }
    
    // 4. 串行四分NTT
    std::cout << "4. 串行四分NTT:" << std::endl;
    memset(c, 0, (2 * n - 1) * sizeof(int));
    start = std::chrono::high_resolution_clock::now();
    four_div_ntt_multiply(a, b, c, n, p);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "   耗时: " << duration << " us" << std::endl;
    fCheck(c, n, current_input_id); // 使用 current_input_id
    
    // 5. pthread四分NTT
    std::cout << "5. pthread四分NTT:" << std::endl;
    for (int i = 0; i < num_threads_count; i++) {
        int num_threads = num_threads_list[i];
        memset(c, 0, (2 * n - 1) * sizeof(int));
        start = std::chrono::high_resolution_clock::now();
        pthread_four_div_ntt_multiply(a, b, c, n, p, num_threads);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "   " << num_threads << "线程耗时: " << duration << " us" << std::endl;
        fCheck(c, n, current_input_id); // 使用 current_input_id
    }
    
    // 6. OpenMP四分NTT
    std::cout << "6. OpenMP四分NTT:" << std::endl;
    for (int i = 0; i < num_threads_count; i++) {
        int num_threads = num_threads_list[i];
        memset(c, 0, (2 * n - 1) * sizeof(int));
        start = std::chrono::high_resolution_clock::now();
        openmp_four_div_ntt_multiply(a, b, c, n, p, num_threads);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "   " << num_threads << "线程耗时: " << duration << " us" << std::endl;
        fCheck(c, n, current_input_id); // 使用 current_input_id
    }
    
    // 7. pthread CRT NTT
    std::cout << "7. pthread CRT NTT:" << std::endl;
    memset(c, 0, (2 * n - 1) * sizeof(int));
    start = std::chrono::high_resolution_clock::now();
    pthread_crt_ntt_multiply(a, b, c, n, p);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "   耗时: " << duration << " us" << std::endl;
    fCheck(c, n, current_input_id); // 使用 current_input_id
    
    // 8. OpenMP CRT NTT
    std::cout << "8. OpenMP CRT NTT:" << std::endl;
    memset(c, 0, (2 * n - 1) * sizeof(int));
    start = std::chrono::high_resolution_clock::now();
    openmp_crt_ntt_multiply(a, b, c, n, p);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "   耗时: " << duration << " us" << std::endl;
    fCheck(c, n, current_input_id); // 使用 current_input_id
}
// 为大整数类型添加快速幂函数
cpp_int mod_pow_big(cpp_int a, cpp_int b, cpp_int p) {
    cpp_int res = 1;
    a %= p;
    while (b > 0) {
        if (b & 1) res = res * a % p;
        a = a * a % p;
        b >>= 1;
    }
    return res;
}

// 为大整数类型添加求逆元函数
cpp_int mod_inverse_big(cpp_int a, cpp_int p) {
    return mod_pow_big(a, p - 2, p);
}
// 添加大模数NTT的支持函数
void big_mod_ntt_multiply(int *a, int *b, int *c, int n, cpp_int p_target_big) { // Renamed p to p_target_big
    int *results[MOD_COUNT];
    int output_len = 2 * n - 1;

    for (int i = 0; i < MOD_COUNT; i++) {
        results[i] = new int[output_len];
    }
    // std::cout << "   使用" << MOD_COUNT << "个固定线程进行大整数CRT计算" << std::endl;
    #pragma omp parallel for num_threads(MOD_COUNT)
    for (int i = 0; i < MOD_COUNT; i++) {
        int len_ntt = get_min_power_of_2(output_len);
        int mod_current = MODULI[i];
        
        int *ta = new int[len_ntt];
        int *tb = new int[len_ntt];
        
        for (int j = 0; j < n; j++) {
            long long val_a = a[j];
            long long val_b = b[j];
            ta[j] = static_cast<int>((val_a % mod_current + mod_current) % mod_current);
            tb[j] = static_cast<int>((val_b % mod_current + mod_current) % mod_current);
        }
        std::fill(ta + n, ta + len_ntt, 0);
        std::fill(tb + n, tb + len_ntt, 0);
        
        ntt(ta, len_ntt, 1, mod_current);
        ntt(tb, len_ntt, 1, mod_current);
        
        for (int j = 0; j < len_ntt; j++) {
            ta[j] = static_cast<int>((1LL * ta[j] * tb[j]) % mod_current);
             if (ta[j] < 0) ta[j] += mod_current;
        }
        
        ntt(ta, len_ntt, -1, mod_current);
        
        memcpy(results[i], ta, output_len * sizeof(int));
        
        delete[] ta;
        delete[] tb;
    }
    
    // 使用CRT合并不同模数下的结果, 与 crt_combine 函数逻辑一致
    const CrtConstants& consts = get_crt_constants();

    // #pragma omp parallel for // If output_len is very large
    for (int i = 0; i < output_len; i++) {
        cpp_int current_sum_mod_M = 0;
        for (int j = 0; j < MOD_COUNT; j++) {
            cpp_int val_j = results[j][i];
            val_j = (val_j % MODULI[j] + MODULI[j]) % MODULI[j];
            
            cpp_int single_term = val_j;
            single_term *= consts.M_divs[j];
            single_term *= consts.M_invs[j];
            current_sum_mod_M += single_term;
        }
        current_sum_mod_M %= consts.M;
        if (current_sum_mod_M < 0) {
            current_sum_mod_M += consts.M;
        }
        
        cpp_int final_val = current_sum_mod_M % p_target_big;
        if (final_val < 0) {
            final_val += p_target_big;
        }
        c[i] = static_cast<int>(final_val); // Assumes final result fits in int
    }
    
    for (int i = 0; i < MOD_COUNT; i++) {
        delete[] results[i];
    }
}

int a[300000], b[300000], ab[300000];
int main(int argc, char *argv[]) {
    // 保留原有的初始化代码
    int test_begin = 0;
    int test_end = 1;
    int thread_count = 8;
    
    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 1337006139375617（注意修正的大模数）
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    
    // 参数处理部分保持不变
    if (argc > 1) {
        test_begin = atoi(argv[1]);
        test_end = test_begin;
    }
    if (argc > 2) {
        thread_count = atoi(argv[2]);
    }
    
    std::cout << "运行测试 " << test_begin << " 到 " << test_end 
              << " 使用 " << thread_count << " 个线程" << std::endl;
    
    // 对每个测试用例运行所有算法
for(int i = test_begin; i <= test_end; ++i) {
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        
        std::cout << "测试用例 " << i << "：n = " << n_ << ", p = " << p_ << std::endl;
        
        memset(ab, 0, (2*n_ -1) * sizeof(int)); // 清理ab以准备接收结果
        test_performance(a, b, ab, n_, p_, i); // 传递 i 作为 current_input_id
       
        // 如果是大模数，额外测试大模数算法
        if (p_ > INT_MAX || p_ == 1337006139375617LL) {
            std::cout << "\n大模数专用算法测试:" << std::endl;
            memset(ab, 0, sizeof(ab));
            
            auto start = std::chrono::high_resolution_clock::now();
            cpp_int big_p = p_;
            if (p_ == 1337006139375617LL) {
                big_p = cpp_int("1337006139375617");
            }
            big_mod_ntt_multiply(a, b, ab, n_, big_p);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            
            std::cout << "   大模数NTT耗时: " << duration << " us" << std::endl;
            fCheck(ab, n_, i);
        }
        
        // 修改此处，使用最后一个测试的算法结果作为最终结果
        // 不再单独执行某个算法
        std::cout << "\n使用最后一次测试的算法结果作为最终结果" << std::endl;
        
        // test_performance中最后测试的是OpenMP CRT NTT，其结果已经存储在ab中
        // 不需要重新计算，直接检查并写入
        fCheck(ab, n_, i);
        fWrite(ab, n_, i);
    }
    
    return 0;
}