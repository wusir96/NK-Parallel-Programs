#include <mpi.h>
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
#include <pthread.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <thread>
using namespace boost::multiprecision;

// MPI全局变量
int mpi_rank, mpi_size;

cpp_int mod_pow_big(cpp_int a, cpp_int b, cpp_int p);
cpp_int mod_inverse_big(cpp_int a, cpp_int p);

//============================ 数据读写函数 ============================
void fRead(long long *a, long long *b, int *n, long long*p, int input_id){
    if (mpi_rank == 0) {
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
        fin.close();
    }
    
    // 广播数据到所有进程
    MPI_Bcast(n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(p, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(a, *n, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, *n, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
}

void fCheck(long long *ab, int n, int input_id){
    if (mpi_rank != 0) return;
    
    std::string str1 = "/nttdata/"; 
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
    if (mpi_rank != 0) return;
    
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

//============================ Barrett规约优化 ============================
struct BarrettParams {
    long long q;          // 模数
    __uint128_t r;        // 预计算的 r = floor(2^64 / q)
    int k;                // k = 64
    
    BarrettParams(long long modulus) : q(modulus), k(64) {
        r = ((__uint128_t(1) << 64) / q);
    }
};

// Barrett模乘优化版本
inline long long barrett_reduce(long long x, const BarrettParams& params) {
    if (x < params.q) return x;
    
    __uint128_t tmp = (__uint128_t)x * params.r;
    long long quotient = tmp >> 64;
    long long result = x - quotient * params.q;
    
    if (result >= params.q) result -= params.q;
    if (result < 0) result += params.q;
    
    return result;
}

// Barrett模乘
inline long long barrett_mul_mod(long long a, long long b, const BarrettParams& params) {
    __uint128_t product = (__uint128_t)a * b;
    __uint128_t tmp = (product * params.r) >> 64;
    long long result = product - tmp * params.q;
    
    if (result >= params.q) result -= params.q;
    if (result < 0) result += params.q;
    
    return result;
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

//============================ Barrett优化的NTT ============================
void barrett_ntt(long long *a, int len, int type, long long p) {
    BarrettParams params(p);
    bit_reverse(a, len);
    
    for (int h = 2; h <= len; h <<= 1) {
        long long gn = mod_pow(3, (p - 1) / h, p);
        if (type == -1) gn = mod_inverse(gn, p);
        
        for (int j = 0; j < len; j += h) {
            long long g = 1;
            for (int k = j; k < j + h / 2; k++) {
                long long u = a[k];
                long long v = barrett_mul_mod(g, a[k + h / 2], params);
                a[k] = barrett_reduce(u + v, params);
                a[k + h / 2] = barrett_reduce(u - v + p, params);
                g = barrett_mul_mod(g, gn, params);
            }
        }
    }
    
    if (type == -1) {
        long long inv_len = mod_inverse(len, p);
        for (int i = 0; i < len; i++) {
            a[i] = barrett_mul_mod(a[i], inv_len, params);
        }
    }
}

//============================ pthread Barrett NTT ============================
struct ThreadArgs {
    long long *a;
    int len;
    int h;
    long long p;
    long long gn;
    pthread_barrier_t *barrier;
    int thread_id;
    int num_threads;
    BarrettParams *params;
};

void* pthread_barrett_ntt_worker(void *arg) {
    ThreadArgs *args = (ThreadArgs*)arg;
    long long *a = args->a;
    int len = args->len;
    int h = args->h;
    long long p = args->p;
    long long gn = args->gn;
    int thread_id = args->thread_id;
    int num_threads = args->num_threads;
    BarrettParams *params = args->params;

    for (int j = thread_id * h; j < len; j += num_threads * h) {
        long long g = 1;
        for (int k = j; k < j + h / 2; k++) {
            if (k + h/2 >= len) continue;
            long long u = a[k];
            long long v = barrett_mul_mod(g, a[k + h / 2], *params);
            a[k] = barrett_reduce(u + v, *params);
            a[k + h / 2] = barrett_reduce(u - v + p, *params);
            g = barrett_mul_mod(g, gn, *params);
        }
    }
    
    if (args->barrier) {
        pthread_barrier_wait(args->barrier);
    }
    return NULL;
}

void pthread_barrett_ntt(long long *a, int len, int type, long long p, int num_threads_requested) {
    if (num_threads_requested <= 0) num_threads_requested = 1;
    int num_threads = std::min(num_threads_requested, len / 2);
    if (num_threads < 1) num_threads = 1;

    BarrettParams params(p);
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
            barrett_ntt(a, len, type, p);
            return;
        }
        barrier_initialized = true;
    }
    
    for (int h = 2; h <= len; h <<= 1) {
        long long gn_h = mod_pow(3, (p - 1) / h, p);
        if (type == -1) gn_h = mod_inverse(gn_h, p);
        
        if (num_threads > 1) {
            for (int t = 0; t < num_threads; t++) {
                args[t] = {a, len, h, p, gn_h, &barrier, t, num_threads, &params};
                if (pthread_create(&threads[t], NULL, pthread_barrett_ntt_worker, &args[t]) != 0) {
                    std::cerr << "Failed to create thread " << t << std::endl;
                    for(int i=0; i<t; ++i) pthread_join(threads[i], NULL);
                    if(barrier_initialized) pthread_barrier_destroy(&barrier);
                    delete[] threads;
                    delete[] args;
                    barrett_ntt(a, len, type, p);
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
                    long long v = barrett_mul_mod(g, a[k + h / 2], params);
                    a[k] = barrett_reduce(u + v, params);
                    a[k + h / 2] = barrett_reduce(u - v + p, params);
                    g = barrett_mul_mod(g, gn_h, params);
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
                BarrettParams *params_lambda = args_lambda->params;

                int chunk_size = (len_lambda + num_threads_lambda - 1) / num_threads_lambda;
                int start_idx = thread_id_lambda * chunk_size;
                int end_idx = std::min(start_idx + chunk_size, len_lambda);

                for (int i = start_idx; i < end_idx; i++) {
                    arr_lambda[i] = barrett_mul_mod(arr_lambda[i], inv_val_lambda, *params_lambda);
                }
                return NULL;
            };

            for (int t = 0; t < num_threads; t++) {
                args[t] = {a, len, 0, p, inv_len_val, NULL, t, num_threads, &params};
                if(pthread_create(&threads[t], NULL, scaling_worker, &args[t])!=0){
                    std::cerr << "Failed to create scaling thread " << t << std::endl;
                    for(int i=0; i<t; ++i) pthread_join(threads[i], NULL);
                    if(barrier_initialized) pthread_barrier_destroy(&barrier);
                    delete[] threads;
                    delete[] args;
                    for(int i=0; i<len; ++i) a[i] = barrett_mul_mod(a[i], inv_len_val, params);
                    return;
                }
            }
            
            for (int t = 0; t < num_threads; t++) {
                pthread_join(threads[t], NULL);
            }
        } else {
            for (int i = 0; i < len; i++) {
                a[i] = barrett_mul_mod(a[i], inv_len_val, params);
            }
        }
    }
    
    if (barrier_initialized) {
        pthread_barrier_destroy(&barrier);
    }
    delete[] threads;
    delete[] args;
}

void pthread_barrett_ntt_multiply(long long *a_orig, long long *b_orig, long long *c_final, int n, long long p, int num_threads_requested) {
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
    
    pthread_barrett_ntt(ta, len, 1, p, num_threads_requested);
    pthread_barrett_ntt(tb, len, 1, p, num_threads_requested);
    
    BarrettParams params(p);
    if (num_threads_requested > 1) {
        std::vector<std::thread> pointwise_threads;
        int num_actual_threads_pm = std::min(num_threads_requested, len);
        if(num_actual_threads_pm < 1) num_actual_threads_pm = 1;

        for (int t = 0; t < num_actual_threads_pm; t++) {
            pointwise_threads.emplace_back([t, num_actual_threads_pm, len, ta, tb, &params]() {
                int chunk_size = (len + num_actual_threads_pm - 1) / num_actual_threads_pm;
                int start_idx = t * chunk_size;
                int end_idx = std::min(start_idx + chunk_size, len);
                for (int i = start_idx; i < end_idx; i++) {
                    ta[i] = barrett_mul_mod(ta[i], tb[i], params);
                }
            });
        }
        for (auto& th : pointwise_threads) {
            th.join();
        }
    } else {
        for (int i = 0; i < len; i++) {
            ta[i] = barrett_mul_mod(ta[i], tb[i], params);
        }
    }
    
    pthread_barrett_ntt(ta, len, -1, p, num_threads_requested);
    
    std::copy(ta, ta + (2 * n - 1), c_final);
    
    delete[] ta;
    delete[] tb;
}

//============================ CRT多线程实现 ============================
const int MOD_COUNT = 4;
const int MODULI[MOD_COUNT] = {469762049, 998244353, 167772161, 1004535809};

struct CrtConstants {
    cpp_int M;
    cpp_int M_divs[MOD_COUNT];
    cpp_int M_invs[MOD_COUNT];
    
    CrtConstants() {
        M = 1;
        for (int i = 0; i < MOD_COUNT; i++) {
            M *= MODULI[i];
        }
        
        for (int i = 0; i < MOD_COUNT; i++) {
            M_divs[i] = M / MODULI[i];
            cpp_int M_div_mod_mi = M_divs[i] % cpp_int(MODULI[i]);
            if (M_div_mod_mi < 0) M_div_mod_mi += MODULI[i];
            M_invs[i] = mod_inverse_big(M_div_mod_mi, cpp_int(MODULI[i]));
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

void* pthread_crt_barrett_ntt_worker(void *arg) {
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
    
    barrett_ntt(ta, len, 1, p_mod);
    barrett_ntt(tb, len, 1, p_mod);
    
    BarrettParams params(p_mod);
    for (int i = 0; i < len; i++) {
        ta[i] = barrett_mul_mod(ta[i], tb[i], params);
    }
    
    barrett_ntt(ta, len, -1, p_mod);
    
    for (int i = 0; i < (2 * n_val - 1); ++i) {
        c_output[i] = ta[i];
    }
    
    delete[] ta;
    delete[] tb;
    return NULL;
}

void crt_combine(long long*results[], int n_orig, long long p_target, long long *c_final) {
    const CrtConstants& consts = get_crt_constants(); 
    int result_len = 2 * n_orig - 1;

    for (int i = 0; i < result_len; i++) {
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
        
        cpp_int final_value_mod_p = current_sum_mod_M % cpp_int(p_target);
        if (final_value_mod_p < 0) {
            final_value_mod_p += p_target;
        }
        
        c_final[i] = static_cast<long long>(final_value_mod_p);
    }
}

void pthread_crt_barrett_ntt_multiply(long long *a, long long *b, long long *c, int n, long long p_target) {
    pthread_t threads[MOD_COUNT];
    CrtThreadArgs args[MOD_COUNT];
    long long *results[MOD_COUNT];
    
    for (int i = 0; i < MOD_COUNT; i++) {
        results[i] = new long long[2 * n - 1];
        args[i] = {a, b, n, MODULI[i], results[i]};
        if (pthread_create(&threads[i], NULL, pthread_crt_barrett_ntt_worker, &args[i]) != 0) {
            std::cerr << "pthread_crt_barrett_ntt_multiply: Failed to create thread " << i << std::endl;
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

//============================ MPI CRT实现 ============================
void mpi_crt_barrett_ntt_multiply(long long *a, long long *b, long long *c, int n, long long p_target) {
    int mods_per_process = (MOD_COUNT + mpi_size - 1) / mpi_size;
    int start_mod = mpi_rank * mods_per_process;
    int end_mod = std::min(start_mod + mods_per_process, MOD_COUNT);
    
    // 每个进程处理分配给它的模数
    std::vector<long long*> local_results;
    for (int i = start_mod; i < end_mod; i++) {
        long long *result = new long long[2 * n - 1];
        local_results.push_back(result);
        
        // 在当前模数下进行Barrett NTT运算
        int len = get_min_power_of_2(2 * n - 1);
        long long *ta = new long long[len];
        long long *tb = new long long[len];
        
        for (int j = 0; j < n; j++) {
            long long val_a = a[j];
            long long val_b = b[j];
            ta[j] = ((val_a % MODULI[i]) + MODULI[i]) % MODULI[i];
            tb[j] = ((val_b % MODULI[i]) + MODULI[i]) % MODULI[i];
        }
        std::fill(ta + n, ta + len, 0);
        std::fill(tb + n, tb + len, 0);
        
        barrett_ntt(ta, len, 1, MODULI[i]);
        barrett_ntt(tb, len, 1, MODULI[i]);
        
        BarrettParams params(MODULI[i]);
        for (int j = 0; j < len; j++) {
            ta[j] = barrett_mul_mod(ta[j], tb[j], params);
        }
        
        barrett_ntt(ta, len, -1, MODULI[i]);
        
        for (int k = 0; k < (2 * n - 1); ++k) {
            result[k] = ta[k];
        }
        
        delete[] ta;
        delete[] tb;
    }
    
    // 收集所有结果到rank 0
    long long *all_results[MOD_COUNT];
    for (int i = 0; i < MOD_COUNT; i++) {
        all_results[i] = new long long[2 * n - 1];
    }
    
    // 每个进程将其结果发送到rank 0
    for (int proc = 0; proc < mpi_size; proc++) {
        int proc_start_mod = proc * mods_per_process;
        int proc_end_mod = std::min(proc_start_mod + mods_per_process, MOD_COUNT);
        
        for (int mod_idx = proc_start_mod; mod_idx < proc_end_mod; mod_idx++) {
            if (mpi_rank == proc) {
                // 发送数据
                int local_idx = mod_idx - start_mod;
                if (proc == 0) {
                    // rank 0直接复制
                    memcpy(all_results[mod_idx], local_results[local_idx], (2 * n - 1) * sizeof(long long));
                } else {
                    MPI_Send(local_results[local_idx], 2 * n - 1, MPI_LONG_LONG, 0, mod_idx, MPI_COMM_WORLD);
                }
            } else if (mpi_rank == 0) {
                // rank 0接收数据
                MPI_Recv(all_results[mod_idx], 2 * n - 1, MPI_LONG_LONG, proc, mod_idx, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    
    // rank 0进行CRT合并
    if (mpi_rank == 0) {
        crt_combine(all_results, n, p_target, c);
    }
    
    // 将结果广播到所有进程
    MPI_Bcast(c, 2 * n - 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    
    // 清理内存
    for (auto& result : local_results) {
        delete[] result;
    }
    for (int i = 0; i < MOD_COUNT; i++) {
        delete[] all_results[i];
    }
}

//============================ 大模数支持 ============================
void big_mod_barrett_ntt_multiply(long long *a, long long *b, long long *c, int n, cpp_int p_target_big) {
    // 使用MPI CRT方法处理大模数
    long long p_ll = static_cast<long long>(p_target_big);
    mpi_crt_barrett_ntt_multiply(a, b, c, n, p_ll);
}

//============================ 测试函数 ============================
void test_performance(long long *a, long long *b, long long *c, int n, long long p, int current_input_id) {
    int num_threads_list[] = {1, 2, 4, 8};
    int num_threads_count = sizeof(num_threads_list) / sizeof(int);
    
    std::chrono::high_resolution_clock::time_point start, end;
    long long duration;

    if (mpi_rank == 0) {
        std::cout << "\n===== Barrett规约优化MPI+pthread NTT性能测试 (n=" << n << ", p=" << p << ", input_id=" << current_input_id << ") =====" << std::endl;
        std::cout << "使用 " << mpi_size << " 个MPI进程" << std::endl;
    }
    
    // 1. pthread Barrett NTT (单进程测试)
    if (mpi_rank == 0) {
        std::cout << "1. pthread Barrett NTT:" << std::endl;
        for (int i = 0; i < num_threads_count; i++) {
            int num_threads = num_threads_list[i];
            memset(c, 0, (2 * n - 1) * sizeof(long long));
            start = std::chrono::high_resolution_clock::now();
            pthread_barrett_ntt_multiply(a, b, c, n, p, num_threads);
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            std::cout << "   " << num_threads << "线程耗时: " << duration << " us" << std::endl;
            fCheck(c, n, current_input_id);
        }
    }
    
    // 2. pthread CRT Barrett NTT (单进程测试)
    if (mpi_rank == 0) {
        std::cout << "2. pthread CRT Barrett NTT:" << std::endl;
        memset(c, 0, (2 * n - 1) * sizeof(long long));
        start = std::chrono::high_resolution_clock::now();
        pthread_crt_barrett_ntt_multiply(a, b, c, n, p);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "   耗时: " << duration << " us" << std::endl;
        fCheck(c, n, current_input_id);
    }
    
    // 3. MPI CRT Barrett NTT (多进程测试)
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        std::cout << "3. MPI CRT Barrett NTT:" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }
    
    memset(c, 0, (2 * n - 1) * sizeof(long long));
    mpi_crt_barrett_ntt_multiply(a, b, c, n, p);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "   耗时: " << duration << " us" << std::endl;
        fCheck(c, n, current_input_id);
    }
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
    // 初始化MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    int test_begin = 0;
    int test_end = 0;
    
    if (argc > 1) {
        test_begin = atoi(argv[1]);
        test_end = test_begin;
    }
    
    if (mpi_rank == 0) {
        std::cout << "运行Barrett规约优化MPI测试 " << test_begin << " 到 " << test_end 
                  << " 使用 " << mpi_size << " 个MPI进程" << std::endl;
    }
    
    for(int i = test_begin; i <= test_end; ++i) {
        int n_;
        long long p_;
        fRead(a, b, &n_, &p_, i);
        
        if (mpi_rank == 0) {
            std::cout << "测试用例 " << i << "：n = " << n_ << ", p = " << p_ << std::endl;
        }
        
        memset(ab, 0, (2*n_ -1) * sizeof(long long));
        test_performance(a, b, ab, n_, p_, i); 
       
        if (p_ > INT_MAX || p_ == 1337006139375617LL) {
            if (mpi_rank == 0) {
                std::cout << "\n大模数专用Barrett算法测试:" << std::endl;
            }
            memset(ab, 0, (2*n_ -1) * sizeof(long long));
            
            MPI_Barrier(MPI_COMM_WORLD);
            auto start = std::chrono::high_resolution_clock::now();
            cpp_int big_p = p_;
            if (p_ == 1337006139375617LL) {
                 big_p = cpp_int("1337006139375617");
            }

            big_mod_barrett_ntt_multiply(a, b, ab, n_, big_p);
            MPI_Barrier(MPI_COMM_WORLD);
            
            if (mpi_rank == 0) {
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
                std::cout << "   大模数MPI Barrett NTT耗时: " << duration << " us" << std::endl;
                fCheck(ab, n_, i);
            }
        }
        
        if (mpi_rank == 0) {
            std::cout << "\n使用最后一次测试的算法结果作为最终结果" << std::endl;
            fCheck(ab, n_, i);
            fWrite(ab, n_, i);
        }
    }
    
    MPI_Finalize();
    return 0;
}