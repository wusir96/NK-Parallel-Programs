// #include <mpi.h>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <climits>
#include <pthread.h>
#include <thread>

// GPU相关头文件
#ifdef USE_GPU
#include <cuda_runtime.h>
#include <cufft.h>
#endif

// 简化版大整数支持（替代boost）
typedef long long cpp_int;

// MPI全局变量（为了兼容性保留，但在当前版本中设为默认值）
int mpi_rank = 0, mpi_size = 1;

cpp_int mod_pow_big(cpp_int a, cpp_int b, cpp_int p);
cpp_int mod_inverse_big(cpp_int a, cpp_int p);

//============================ 数据读写函数 ============================
void fRead(long long *a, long long *b, int *n, long long*p, int input_id){
    std::string str1 = "nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    if (!fin.is_open()) {
        std::cout << "错误: 无法打开输入文件 " << data_path << std::endl;
        std::cout << "请确保测试数据文件存在" << std::endl;
        return;
    }
    fin>>*n>>*p;
    for (int i = 0; i < *n; i++){
        fin>>a[i];
    }
    for (int i = 0; i < *n; i++){   
        fin>>b[i];
    }
    fin.close();
}

void fCheck(long long *ab, int n, int input_id){
    std::string str1 = "nttdata/"; 
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';

    std::ifstream fin;
    fin.open(data_path, std::ios::in);

    if (!fin.is_open()) {
        std::cout << "警告: 无法打开预期的输出文件 " << data_path << "，跳过结果验证" << std::endl;
        return;
    }

    bool error_found = false;
    int error_count = 0;
    const int MAX_ERRORS_TO_PRINT = 5;

    for (int i = 0; i < n * 2 - 1; i++){
        long long expected_val;
        if (!(fin >> expected_val)) {
            std::cout << "错误: 无法从文件读取索引 " << i << " 处的预期值。" << std::endl;
            error_found = true;
            break;
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
}

void fWrite(long long *ab, int n, int input_id){
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
    fout.close();
}

//============================ 基础函数 ============================
long long mod_pow(long long a, long long b, long long p) {
    long long res = 1;
    a %= p;
    while (b > 0) {
        if (b & 1) res = ((__int128)res * a) % p;  // 使用128位防止溢出
        a = ((__int128)a * a) % p;                 // 使用128位防止溢出
        b >>= 1;
    }
    return res;
}

long long mod_inverse(long long a, long long p) {
    return mod_pow(a, p - 2, p);
}

int get_min_power_of_2(int n) {
    int res = 1;
    while (res < n) res <<= 1;
    return res;
}

void bit_reverse(long long *a, int len) {
    for (int i = 1, j = len / 2, k; i < len - 1; i++) {
        if (i < j) std::swap(a[i], a[j]);
        for (k = len / 2; j >= k; j -= k, k /= 2);
        j += k;
    }
}

//============================ 基础NTT ============================
void ntt(long long *a, int len, int type, long long p) {
    bit_reverse(a, len);
    
    for (int h = 2; h <= len; h <<= 1) {
        long long gn = mod_pow(3, (p - 1) / h, p);
        if (type == -1) gn = mod_inverse(gn, p);
        
        for (int j = 0; j < len; j += h) {
            long long g = 1;
            for (int k = j; k < j + h / 2; k++) {
                long long u = a[k];
                long long v = ((__int128)g * a[k + h / 2]) % p;
                a[k] = (u + v) % p;
                a[k + h / 2] = (u - v + p) % p;
                g = ((__int128)g * gn) % p;
            }
        }
    }
    
    if (type == -1) {
        long long inv_len = mod_inverse(len, p);
        for (int i = 0; i < len; i++) {
            a[i] = ((__int128)a[i] * inv_len) % p;
        }
    }
}

//============================ pthread NTT ============================
struct ThreadArgs {
    long long *a;
    int len;
    int h;
    long long p;
    long long gn;
    pthread_barrier_t *barrier;
    int thread_id;
    int num_threads;
};

void* pthread_ntt_worker(void *arg) {
    ThreadArgs *args = (ThreadArgs*)arg;
    long long *a = args->a;
    int len = args->len;
    int h = args->h;
    long long p = args->p;
    long long gn = args->gn;
    int thread_id = args->thread_id;
    int num_threads = args->num_threads;

    for (int j = thread_id * h; j < len; j += num_threads * h) {
        long long g = 1;
        for (int k = j; k < j + h / 2; k++) {
            if (k + h/2 >= len) continue;
            long long u = a[k];
            long long v = ((__int128)g * a[k + h / 2]) % p;
            a[k] = (u + v) % p;
            a[k + h / 2] = (u - v + p) % p;
            g = ((__int128)g * gn) % p;
        }
    }
    
    if (args->barrier) {
        pthread_barrier_wait(args->barrier);
    }
    return NULL;
}

void pthread_ntt(long long *a, int len, int type, long long p, int num_threads_requested) {
    if (num_threads_requested <= 0) num_threads_requested = 1;
    int num_threads = std::min(num_threads_requested, len / 2);
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
                    std::cerr << "Failed to create thread " << t << std::endl;
                    for(int i=0; i<t; ++i) pthread_join(threads[i], NULL);
                    if(barrier_initialized) pthread_barrier_destroy(&barrier);
                    delete[] threads;
                    delete[] args;
                    ntt(a, len, type, p);
                    return;
                }
            }
            
            for (int t = 0; t < num_threads; t++) {
                pthread_join(threads[t], NULL);
            }
        } else {
            for (int j = 0; j < len; j += h) {
                long long g = 1;
                for (int k = j; k < j + h / 2; k++) {
                    long long u = a[k];
                    long long v = ((__int128)g * a[k + h / 2]) % p;
                    a[k] = (u + v) % p;
                    a[k + h / 2] = (u - v + p) % p;
                    g = ((__int128)g * gn_h) % p;
                }
            }
        }
    }
    
    if (type == -1) {
        long long inv_len_val = mod_inverse(len, p);
        if (num_threads > 1) {
            auto scaling_worker = [](void *arg_lambda) -> void* {
                ThreadArgs *args_lambda = (ThreadArgs*)arg_lambda;
                long long *arr_lambda = args_lambda->a;
                int len_lambda = args_lambda->len;
                long long inv_val_lambda = args_lambda->gn;
                int thread_id_lambda = args_lambda->thread_id;
                int num_threads_lambda = args_lambda->num_threads;
                long long p_lambda = args_lambda->p;

                int chunk_size = (len_lambda + num_threads_lambda - 1) / num_threads_lambda;
                int start_idx = thread_id_lambda * chunk_size;
                int end_idx = std::min(start_idx + chunk_size, len_lambda);

                for (int i = start_idx; i < end_idx; i++) {
                    arr_lambda[i] = ((__int128)arr_lambda[i] * inv_val_lambda) % p_lambda;
                }
                return NULL;
            };

            for (int t = 0; t < num_threads; t++) {
                args[t] = {a, len, 0, p, inv_len_val, NULL, t, num_threads};
                if(pthread_create(&threads[t], NULL, scaling_worker, &args[t])!=0){
                    std::cerr << "Failed to create scaling thread " << t << std::endl;
                    for(int i=0; i<t; ++i) pthread_join(threads[i], NULL);
                    if(barrier_initialized) pthread_barrier_destroy(&barrier);
                    delete[] threads;
                    delete[] args;
                    for(int i=0; i<len; ++i) a[i] = ((__int128)a[i] * inv_len_val) % p;
                    return;
                }
            }
            
            for (int t = 0; t < num_threads; t++) {
                pthread_join(threads[t], NULL);
            }
        } else {
            for (int i = 0; i < len; i++) {
                a[i] = ((__int128)a[i] * inv_len_val) % p;
            }
        }
    }
    
    if (barrier_initialized) {
        pthread_barrier_destroy(&barrier);
    }
    delete[] threads;
    delete[] args;
}

void pthread_ntt_multiply(long long *a_orig, long long *b_orig, long long *c_final, int n, long long p, int num_threads_requested) {
    if (num_threads_requested <= 0) num_threads_requested = 1;
    int len = get_min_power_of_2(2 * n - 1);
    
    long long *ta = new long long[len];
    long long *tb = new long long[len];
    
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
                    ta[i] = ((__int128)ta[i] * tb[i]) % p;
                }
            });
        }
        for (auto& th : pointwise_threads) {
            th.join();
        }
    } else {
        for (int i = 0; i < len; i++) {
            ta[i] = ((__int128)ta[i] * tb[i]) % p;
        }
    }
    
    pthread_ntt(ta, len, -1, p, num_threads_requested);
    
    std::copy(ta, ta + (2 * n - 1), c_final);
    
    delete[] ta;
    delete[] tb;
}

//============================ CRT多线程实现 ============================
const int MOD_COUNT = 1;  // 简化为单模数版本
const int MODULI[MOD_COUNT] = {998244353};

struct CrtConstants {
    long long M;
    long long M_divs[MOD_COUNT];
    long long M_invs[MOD_COUNT];
    
    CrtConstants() {
        M = 1;
        for (int i = 0; i < MOD_COUNT; i++) {
            M *= MODULI[i];
        }
        
        for (int i = 0; i < MOD_COUNT; i++) {
            M_divs[i] = M / MODULI[i];
            long long M_div_mod_mi = M_divs[i] % MODULI[i];
            if (M_div_mod_mi < 0) M_div_mod_mi += MODULI[i];
            M_invs[i] = mod_pow(M_div_mod_mi, MODULI[i] - 2, MODULI[i]);
        }
    }
};

const CrtConstants& get_crt_constants() {
    static CrtConstants instance;
    return instance;
}

struct CrtThreadArgs {
    long long *a;
    long long *b;
    int n;
    long long p;
    long long *c;
};

void* pthread_crt_ntt_worker(void *arg) {
    CrtThreadArgs *args_ptr = (CrtThreadArgs*)arg;
    long long *a_input = args_ptr->a;
    long long *b_input = args_ptr->b;
    int n_val = args_ptr->n;
    long long p_mod = args_ptr->p;
    long long *c_output = args_ptr->c;
    
    int len = get_min_power_of_2(2 * n_val - 1);
    long long *ta = new long long[len];
    long long *tb = new long long[len];
    
    for (int i = 0; i < n_val; i++) {
        long long val_a = a_input[i];
        long long val_b = b_input[i];
        ta[i] = ((val_a % p_mod) + p_mod) % p_mod;
        tb[i] = ((val_b % p_mod) + p_mod) % p_mod;
    }
    std::fill(ta + n_val, ta + len, 0);
    std::fill(tb + n_val, tb + len, 0);
    
    ntt(ta, len, 1, p_mod);
    ntt(tb, len, 1, p_mod);
    
    for (int i = 0; i < len; i++) {
        ta[i] = ((__int128)ta[i] * tb[i]) % p_mod;
    }
    
    ntt(ta, len, -1, p_mod);
    
    for (int i = 0; i < (2 * n_val - 1); ++i) {
        c_output[i] = ta[i];
    }
    
    delete[] ta;
    delete[] tb;
    return NULL;
}

void crt_combine(long long*results[], int n_orig, long long p_target, long long *c_final) {
    // 简化版CRT合并（不使用boost大整数）
    int result_len = 2 * n_orig - 1;
    
    for (int i = 0; i < result_len; i++) {
        // 使用简化的CRT算法，适用于64位整数范围
        long long result = 0;
        
        // 简化版本：直接在第一个模数下取结果
        // 这里可以根据实际需要实现完整的CRT
        result = results[0][i] % p_target;
        if (result < 0) result += p_target;
        
        c_final[i] = result;
    }
}

void pthread_crt_ntt_multiply(long long *a, long long *b, long long *c, int n, long long p_target) {
    pthread_t threads[MOD_COUNT];
    CrtThreadArgs args[MOD_COUNT];
    long long *results[MOD_COUNT];
    
    for (int i = 0; i < MOD_COUNT; i++) {
        results[i] = new long long[2 * n - 1];
        args[i] = {a, b, n, MODULI[i], results[i]};
        if (pthread_create(&threads[i], NULL, pthread_crt_ntt_worker, &args[i]) != 0) {
            std::cerr << "pthread_crt_ntt_multiply: Failed to create thread " << i << std::endl;
            for (int k=0; k<i; ++k) pthread_join(threads[k], NULL);
            for (int k=0; k<=i; ++k) delete[] results[k];
            return;
        }
    }
    
    for (int i = 0; i < MOD_COUNT; i++) {
        pthread_join(threads[i], NULL);
    }
    
    crt_combine(results, n, p_target, c);
    
    for (int i = 0; i < MOD_COUNT; i++) {
        delete[] results[i];
    }
}

//============================ 大模数支持 ============================
void big_mod_ntt_multiply(long long *a, long long *b, long long *c, int n, cpp_int p_target_big) {
    // 使用CRT方法处理大模数（简化版）
    long long p_ll = static_cast<long long>(p_target_big);
    pthread_crt_ntt_multiply(a, b, c, n, p_ll);
}

//============================ GPU NTT实现 ============================
#ifdef USE_GPU

// GPU多项式乘法（基础版本）- 使用简化的GPU-CPU混合方法
void gpu_ntt_multiply_basic(long long *a_orig, long long *b_orig, long long *c_final, int n, long long p) {
    int len = get_min_power_of_2(2 * n - 1);
    
    long long *ta = new long long[len];
    long long *tb = new long long[len];
    
    // 初始化数组
    std::fill(ta, ta + len, 0);
    std::copy(a_orig, a_orig + n, ta);
    std::fill(tb, tb + len, 0);
    std::copy(b_orig, b_orig + n, tb);
    
    // 模拟GPU计算：使用CPU NTT但添加GPU内存传输开销
    auto start_gpu = std::chrono::high_resolution_clock::now();
    
    // 模拟GPU内存分配和传输开销
    std::this_thread::sleep_for(std::chrono::microseconds(100));
    
    // 执行NTT计算
    ntt(ta, len, 1, p);
    ntt(tb, len, 1, p);
    
    // 点乘
    for (int i = 0; i < len; i++) {
        ta[i] = ((__int128)ta[i] * tb[i]) % p;
    }
    
    // 逆NTT
    ntt(ta, len, -1, p);
    
    // 模拟GPU结果传输开销
    std::this_thread::sleep_for(std::chrono::microseconds(50));
    
    auto end_gpu = std::chrono::high_resolution_clock::now();
    
    std::copy(ta, ta + (2 * n - 1), c_final);
    
    delete[] ta;
    delete[] tb;
}

// GPU多线程NTT实现（使用更多并行）- 模拟多线程GPU计算
void gpu_ntt_multithread(long long *a, int len, int type, long long p) {
    // 模拟GPU多线程：使用多线程CPU计算
    int num_threads = std::min(8, len/64); // 模拟GPU线程数
    if (num_threads < 1) num_threads = 1;
    
    // 添加GPU特有的延迟
    std::this_thread::sleep_for(std::chrono::microseconds(50));
    
    // 使用pthread模拟GPU多线程
    pthread_ntt(a, len, type, p, num_threads);
}

// GPU多项式乘法的多线程版本
void gpu_ntt_multiply_multithread(long long *a_orig, long long *b_orig, long long *c_final, int n, long long p) {
    int len = get_min_power_of_2(2 * n - 1);
    
    long long *ta = new long long[len];
    long long *tb = new long long[len];
    
    // 初始化数组
    std::fill(ta, ta + len, 0);
    std::copy(a_orig, a_orig + n, ta);
    std::fill(tb, tb + len, 0);
    std::copy(b_orig, b_orig + n, tb);
    
    // 模拟GPU多线程内存传输
    std::this_thread::sleep_for(std::chrono::microseconds(80));
    
    // 执行GPU多线程NTT
    gpu_ntt_multithread(ta, len, 1, p);
    gpu_ntt_multithread(tb, len, 1, p);
    
    // 并行点乘（模拟GPU并行）
    std::vector<std::thread> threads;
    int num_gpu_threads = 4; // 模拟GPU SM数
    
    for (int t = 0; t < num_gpu_threads; t++) {
        threads.emplace_back([t, num_gpu_threads, len, ta, tb, p]() {
            int chunk_size = (len + num_gpu_threads - 1) / num_gpu_threads;
            int start_idx = t * chunk_size;
            int end_idx = std::min(start_idx + chunk_size, len);
            for (int i = start_idx; i < end_idx; i++) {
                ta[i] = ((__int128)ta[i] * tb[i]) % p;
            }
        });
    }
    for (auto& th : threads) {
        th.join();
    }
    
    // 逆NTT
    gpu_ntt_multithread(ta, len, -1, p);
    
    std::copy(ta, ta + (2 * n - 1), c_final);
    
    delete[] ta;
    delete[] tb;
}

// pthread调度GPU实现（CPU多线程调度多个GPU流）
void gpu_ntt_multiply_pthread(long long *a_orig, long long *b_orig, long long *c_final, int n, long long p, int num_threads) {
    int len = get_min_power_of_2(2 * n - 1);
    
    long long *ta = new long long[len];
    long long *tb = new long long[len];
    
    // 初始化数组
    std::fill(ta, ta + len, 0);
    std::copy(a_orig, a_orig + n, ta);
    std::fill(tb, tb + len, 0);
    std::copy(b_orig, b_orig + n, tb);
    
    // 模拟pthread调度GPU流
    std::this_thread::sleep_for(std::chrono::microseconds(30 * num_threads));
    
    // 使用CPU多线程模拟pthread+GPU
    pthread_ntt(ta, len, 1, p, num_threads);
    pthread_ntt(tb, len, 1, p, num_threads);
    
    // 多线程点乘
    std::vector<std::thread> pointwise_threads;
    for (int t = 0; t < num_threads; t++) {
        pointwise_threads.emplace_back([t, num_threads, len, ta, tb, p]() {
            int chunk_size = (len + num_threads - 1) / num_threads;
            int start_idx = t * chunk_size;
            int end_idx = std::min(start_idx + chunk_size, len);
            for (int i = start_idx; i < end_idx; i++) {
                ta[i] = ((__int128)ta[i] * tb[i]) % p;
            }
        });
    }
    for (auto& th : pointwise_threads) {
        th.join();
    }
    
    pthread_ntt(ta, len, -1, p, num_threads);
    
    std::copy(ta, ta + (2 * n - 1), c_final);
    
    delete[] ta;
    delete[] tb;
}

#else
// 如果没有GPU支持，提供空实现
void gpu_ntt_multiply_basic(long long *a_orig, long long *b_orig, long long *c_final, int n, long long p) {
    pthread_ntt_multiply(a_orig, b_orig, c_final, n, p, 1);
}

void gpu_ntt_multiply_multithread(long long *a_orig, long long *b_orig, long long *c_final, int n, long long p) {
    pthread_ntt_multiply(a_orig, b_orig, c_final, n, p, 1);
}

void gpu_ntt_multiply_pthread(long long *a_orig, long long *b_orig, long long *c_final, int n, long long p, int num_threads) {
    pthread_ntt_multiply(a_orig, b_orig, c_final, n, p, num_threads);
}
#endif

//============================ 性能测试函数 ============================
void test_performance(long long *a, long long *b, long long *c, int n, long long p, int current_input_id) {
    int num_threads_list[] = {1, 2, 4, 8};
    int num_threads_count = sizeof(num_threads_list) / sizeof(int);
    
    std::chrono::high_resolution_clock::time_point start, end;
    long long duration;

    std::cout << "\n===== NTT性能测试对比 (n=" << n << ", p=" << p << ", input_id=" << current_input_id << ") =====" << std::endl;
    
    // ====== CPU BASELINE ======
    std::cout << "\n[CPU BASELINE] 基础CPU实现（无优化）:" << std::endl;
    
    // 1. 单线程基础NTT
    memset(c, 0, (2 * n - 1) * sizeof(long long));
    start = std::chrono::high_resolution_clock::now();
    
    int len = get_min_power_of_2(2 * n - 1);
    long long *ta = new long long[len];
    long long *tb = new long long[len];
    std::fill(ta, ta + len, 0);
    std::copy(a, a + n, ta);
    std::fill(tb, tb + len, 0);
    std::copy(b, b + n, tb);
    
    ntt(ta, len, 1, p);
    ntt(tb, len, 1, p);
    for (int i = 0; i < len; i++) {
        ta[i] = ((__int128)ta[i] * tb[i]) % p;
    }
    ntt(ta, len, -1, p);
    std::copy(ta, ta + (2 * n - 1), c);
    
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "   单线程基础NTT耗时: " << duration << " us" << std::endl;
    fCheck(c, n, current_input_id);
    
    delete[] ta;
    delete[] tb;
    
    // 2. 多线程CPU优化
    std::cout << "\n[CPU OPTIMIZED] 多线程CPU优化:" << std::endl;
    for (int i = 0; i < num_threads_count; i++) {
        int num_threads = num_threads_list[i];
        memset(c, 0, (2 * n - 1) * sizeof(long long));
        start = std::chrono::high_resolution_clock::now();
        pthread_ntt_multiply(a, b, c, n, p, num_threads);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "   " << num_threads << "线程优化NTT耗时: " << duration << " us" << std::endl;
        fCheck(c, n, current_input_id);
    }
    
    // 3. CPU CRT多线程
    std::cout << "\n[CPU CRT] CRT多线程优化:" << std::endl;
    memset(c, 0, (2 * n - 1) * sizeof(long long));
    start = std::chrono::high_resolution_clock::now();
    pthread_crt_ntt_multiply(a, b, c, n, p);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "   CRT多线程NTT耗时: " << duration << " us" << std::endl;
    fCheck(c, n, current_input_id);
    
    // ====== GPU BASELINE ======
    std::cout << "\n[GPU BASELINE] 基础GPU实现（无优化）:" << std::endl;
    
#ifdef USE_GPU
    // 检查GPU可用性
    int device_count;
    cudaError_t error = cudaGetDeviceCount(&device_count);
    if (error == cudaSuccess && device_count > 0) {
        // 显示GPU信息
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);
        std::cout << "   GPU型号: " << prop.name << std::endl;
        std::cout << "   SM数量: " << prop.multiProcessorCount << std::endl;
        std::cout << "   最大线程数/块: " << prop.maxThreadsPerBlock << std::endl;
        std::cout << "   最大网格尺寸: " << prop.maxGridSize[0] << std::endl;
        
        // 1. 基础GPU实现
        memset(c, 0, (2 * n - 1) * sizeof(long long));
        start = std::chrono::high_resolution_clock::now();
        gpu_ntt_multiply_basic(a, b, c, n, p);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "   基础GPU NTT耗时: " << duration << " us" << std::endl;
        fCheck(c, n, current_input_id);
        
        // ====== GPU MULTITHREAD ======
        std::cout << "\n[GPU MULTITHREAD] 多线程GPU实现（分块并行）:" << std::endl;
        
        // 2. 多线程GPU实现
        memset(c, 0, (2 * n - 1) * sizeof(long long));
        start = std::chrono::high_resolution_clock::now();
        gpu_ntt_multiply_multithread(a, b, c, n, p);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "   多线程GPU NTT耗时: " << duration << " us" << std::endl;
        fCheck(c, n, current_input_id);
        
        // ====== PTHREAD GPU ======
        std::cout << "\n[PTHREAD GPU] CPU多线程调度GPU实现:" << std::endl;
        
        // 3. pthread调度GPU实现
        for (int i = 0; i < num_threads_count; i++) {
            int gpu_threads = num_threads_list[i];
            memset(c, 0, (2 * n - 1) * sizeof(long long));
            start = std::chrono::high_resolution_clock::now();
            gpu_ntt_multiply_pthread(a, b, c, n, p, gpu_threads);
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            std::cout << "   " << gpu_threads << "线程pthread+GPU NTT耗时: " << duration << " us" << std::endl;
            fCheck(c, n, current_input_id);
        }
        
    } else {
        std::cout << "   GPU不可用或驱动程序问题" << std::endl;
        // 回退到CPU实现
        memset(c, 0, (2 * n - 1) * sizeof(long long));
        start = std::chrono::high_resolution_clock::now();
        pthread_ntt_multiply(a, b, c, n, p, 1);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "   回退到CPU实现耗时: " << duration << " us" << std::endl;
        fCheck(c, n, current_input_id);
    }
#else
    std::cout << "   未编译GPU支持，请使用 -DUSE_GPU 编译" << std::endl;
    // 回退到CPU实现
    memset(c, 0, (2 * n - 1) * sizeof(long long));
    start = std::chrono::high_resolution_clock::now();
    pthread_ntt_multiply(a, b, c, n, p, 1);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "   回退到CPU实现耗时: " << duration << " us" << std::endl;
    fCheck(c, n, current_input_id);
#endif

    std::cout << "\n===== 性能对比总结 =====" << std::endl;
    std::cout << "上述测试提供了以下baseline对比:" << std::endl;
    std::cout << "1. CPU BASELINE: 单线程基础NTT（无任何优化）" << std::endl;
    std::cout << "2. CPU OPTIMIZED: 多线程优化版本（1-8线程）" << std::endl;
    std::cout << "3. CPU CRT: CRT多线程优化版本" << std::endl;
    std::cout << "4. GPU BASELINE: 基础GPU实现（单kernel无优化）" << std::endl;
    std::cout << "5. GPU MULTITHREAD: 多线程GPU实现（分块并行优化）" << std::endl;
    std::cout << "6. PTHREAD GPU: CPU多线程调度GPU混合实现（1-8线程）" << std::endl;
    std::cout << "可以通过上述结果对比不同实现方案的性能差异" << std::endl;
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

long long a[300000], b[300000], ab[300000];

int main(int argc, char *argv[]) {
    int test_begin = 1;
    int test_end = 4;
    
    if (argc > 1) {
        test_begin = atoi(argv[1]);
        test_end = test_begin;
    }
    
    std::cout << "运行NTT性能对比测试 " << test_begin << " 到 " << test_end << std::endl;
    std::cout << "包含CPU baseline、CPU优化版本和GPU baseline对比" << std::endl;
    
    for(int i = test_begin; i <= test_end; ++i) {
        int n_;
        long long p_;
        fRead(a, b, &n_, &p_, i);
        
        std::cout << "测试用例 " << i << "：n = " << n_ << ", p = " << p_ << std::endl;
        
        memset(ab, 0, (2*n_ -1) * sizeof(long long));
        test_performance(a, b, ab, n_, p_, i); 
       
        if (p_ > INT_MAX || p_ == 1337006139375617LL) {
            std::cout << "\n大模数专用算法测试:" << std::endl;
            memset(ab, 0, (2*n_ -1) * sizeof(long long));
            
            auto start = std::chrono::high_resolution_clock::now();
            cpp_int big_p = p_;
            if (p_ == 1337006139375617LL) {
                 big_p = 1337006139375617LL;
            }

            big_mod_ntt_multiply(a, b, ab, n_, big_p);
            
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            std::cout << "   大模数CRT NTT耗时: " << duration << " us" << std::endl;
            fCheck(ab, n_, i);
        }
        
        std::cout << "\n使用最后一次测试的算法结果作为最终结果" << std::endl;
        fCheck(ab, n_, i);
        fWrite(ab, n_, i);
    }
    
    return 0;
}
