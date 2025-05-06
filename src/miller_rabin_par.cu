#include "miller_rabin.h"
#include "utils.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <gmp.h>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <random>

#define BLOCK_SIZE 1024
#define MAX_BITS 32768
#define DEFAULT_ITERATIONS 50

using namespace std;


typedef struct {
    double total_time_ms;
    double base_generation_time_ms;
    double decomposition_time_ms;
    double conversion_time_ms;
    double memory_htod_time_ms;
    double kernel_time_ms;
    double memory_dtoh_time_ms;
    double verification_time_ms;
} PerformanceMetrics;

// hehe big int
typedef struct {
    unsigned int data[MAX_BITS/32 + 1];
    int bits;
} CudaBigInt;

// GMP mpz_t to CUDA format
CudaBigInt mpz_to_cuda_bigint(mpz_t num) {
    CudaBigInt result;
    result.bits = mpz_sizeinbase(num, 2);
    memset(result.data, 0, sizeof(result.data));
    size_t count = (result.bits + 31) / 32;
    mpz_export(result.data, &count, -1, sizeof(unsigned int), 0, 0, num);
    
    return result;
}

// CUDA format back to GMP mpz_t
void cuda_bigint_to_mpz(CudaBigInt* cbi, mpz_t result) {
    size_t count = (cbi->bits + 31) / 32;
    mpz_import(result, count, -1, sizeof(unsigned int), 0, 0, cbi->data);
}

__device__ void cuda_mod_mul(CudaBigInt* result, const CudaBigInt* a, const CudaBigInt* b, const CudaBigInt* modulus) {
    memset(result->data, 0, sizeof(unsigned int) * (MAX_BITS/32 + 1));
    result->bits = 0;
    
    unsigned long long carry = 0;
    for (int i = 0; i < (a->bits + 31) / 32; i++) {
        for (int j = 0; j < (b->bits + 31) / 32; j++) {
            unsigned long long product = (unsigned long long)a->data[i] * b->data[j] + carry;
            for (int k = i + j; k < (MAX_BITS/32 + 1); k++) {
                unsigned long long sum = (unsigned long long)result->data[k] + (product & 0xFFFFFFFF);
                result->data[k] = (unsigned int)(sum & 0xFFFFFFFF);
                product >>= 32;
                carry = (sum >> 32) + product;
                if (carry == 0) break;
            }
        }
    }
    
    unsigned int q[MAX_BITS/32 + 1];
    memset(q, 0, sizeof(q));
    
    int r_bits = 0;
    for (int i = MAX_BITS/32; i >= 0; i--) {
        if (result->data[i] > 0) {
            r_bits = i * 32 + 32;
            while ((result->data[i] & (1 << (r_bits % 32 - 1))) == 0) {
                r_bits--;
            }
            break;
        }
    }
    
    if (r_bits <= modulus->bits) {
        result->bits = r_bits;
        return;
    }
    
    for (int shift = r_bits - modulus->bits; shift >= 0; shift--) {
        bool can_subtract = true;
        for (int i = (modulus->bits + shift + 31) / 32 - 1; i >= 0; i--) {
            unsigned int m_val = 0;
            int m_idx = i - shift / 32;
            int bit_shift = shift % 32;
            
            if (m_idx >= 0 && m_idx < (modulus->bits + 31) / 32) {
                m_val = modulus->data[m_idx] << bit_shift;
                if (bit_shift > 0 && m_idx + 1 < (modulus->bits + 31) / 32) {
                    m_val |= modulus->data[m_idx + 1] >> (32 - bit_shift);
                }
            }
            
            if (result->data[i] < m_val) {
                can_subtract = false;
                break;
            } else if (result->data[i] > m_val) {
                break;
            }
        }
        
        if (can_subtract) {
            q[shift / 32] |= 1 << (shift % 32);
            unsigned int borrow = 0;
            for (int i = 0; i < (r_bits + 31) / 32; i++) {
                unsigned int m_val = 0;
                int m_idx = i - shift / 32;
                int bit_shift = shift % 32;
                
                if (m_idx >= 0 && m_idx < (modulus->bits + 31) / 32) {
                    m_val = modulus->data[m_idx] << bit_shift;
                    if (bit_shift > 0 && m_idx + 1 < (modulus->bits + 31) / 32) {
                        m_val |= modulus->data[m_idx + 1] >> (32 - bit_shift);
                    }
                }
                
                unsigned int diff;
                if (result->data[i] >= m_val + borrow) {
                    diff = result->data[i] - m_val - borrow;
                    borrow = 0;
                } else {
                    diff = result->data[i] + (1ULL << 32) - m_val - borrow;
                    borrow = 1;
                }
                result->data[i] = diff;
            }
        }
    }
    
    for (int i = MAX_BITS/32; i >= 0; i--) {
        if (result->data[i] > 0) {
            result->bits = i * 32 + 32;
            while ((result->data[i] & (1 << (result->bits % 32 - 1))) == 0) {
                result->bits--;
            }
            break;
        }
    }
    if (result->bits == 0) result->bits = 1;
}

// modular exponentiation on CUDA
__device__ void cuda_mod_pow(CudaBigInt* result, const CudaBigInt* base, const CudaBigInt* exponent, const CudaBigInt* modulus) {
    memset(result->data, 0, sizeof(unsigned int) * (MAX_BITS/32 + 1));
    result->data[0] = 1;
    result->bits = 1;
    
    bool is_zero = true;
    for (int i = 0; i < (exponent->bits + 31) / 32; i++) {
        if (exponent->data[i] != 0) {
            is_zero = false;
            break;
        }
    }
    if (is_zero) return;
    
    CudaBigInt base_copy;
    memcpy(&base_copy, base, sizeof(CudaBigInt));
    
    for (int bit_pos = exponent->bits - 1; bit_pos >= 0; bit_pos--) {
        CudaBigInt squared;
        cuda_mod_mul(&squared, result, result, modulus);
        memcpy(result, &squared, sizeof(CudaBigInt));
        
        int word_idx = bit_pos / 32;
        int bit_idx = bit_pos % 32;
        if (exponent->data[word_idx] & (1 << bit_idx)) {
            CudaBigInt product;
            cuda_mod_mul(&product, result, &base_copy, modulus);
            memcpy(result, &product, sizeof(CudaBigInt));
        }
    }
}

__global__ void miller_rabin_kernel(CudaBigInt n, CudaBigInt* bases, bool* results, CudaBigInt d, CudaBigInt r, unsigned int iterations) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < iterations) {
        CudaBigInt a = bases[idx];
        bool is_prime = true;
        
        CudaBigInt x;
        cuda_mod_pow(&x, &a, &d, &n);
        
        bool is_one = true;
        for (int i = 1; i < (x.bits + 31) / 32; i++) {
            if (x.data[i] != 0) {
                is_one = false;
                break;
            }
        }
        if (x.data[0] != 1) is_one = false;
        
        CudaBigInt n_minus_1;
        memcpy(&n_minus_1, &n, sizeof(CudaBigInt));
        n_minus_1.data[0]--; 
        for (int i = 0; i < (n_minus_1.bits + 31) / 32; i++) {
            if (n_minus_1.data[i] != 0xFFFFFFFF) break;
            n_minus_1.data[i] = 0;
            if (i + 1 < (n_minus_1.bits + 31) / 32) n_minus_1.data[i + 1]--;
        }
        
        bool is_n_minus_1 = true;
        for (int i = 0; i < (n.bits + 31) / 32; i++) {
            if (x.data[i] != n_minus_1.data[i]) {
                is_n_minus_1 = false;
                break;
            }
        }
        
        if (!is_one && !is_n_minus_1) {
            bool found_minus_1 = false;
            for (unsigned int j = 0; j < r.data[0] - 1; j++) {
                CudaBigInt squared;
                cuda_mod_mul(&squared, &x, &x, &n);
                memcpy(&x, &squared, sizeof(CudaBigInt));
                
                bool is_one = true;
                for (int i = 1; i < (x.bits + 31) / 32; i++) {
                    if (x.data[i] != 0) {
                        is_one = false;
                        break;
                    }
                }
                if (x.data[0] != 1) is_one = false;
                
                if (is_one) {
                    is_prime = false;
                    break;
                }
                
                bool is_n_minus_1 = true;
                for (int i = 0; i < (n.bits + 31) / 32; i++) {
                    if (x.data[i] != n_minus_1.data[i]) {
                        is_n_minus_1 = false;
                        break;
                    }
                }
                
                if (is_n_minus_1) {
                    found_minus_1 = true;
                    break;
                }
            }
            
            if (!found_minus_1) {
                is_prime = false;
            }
        }
        
        results[idx] = is_prime;
    }
}

void print_performance_metrics(const PerformanceMetrics& metrics) {
    printf("===== Performance Metrics =====\n");
    printf("Total execution time:      %.3f ms\n", metrics.total_time_ms);
    printf("Base generation time:      %.3f ms (%.1f%%)\n", metrics.base_generation_time_ms, 
           (metrics.base_generation_time_ms / metrics.total_time_ms) * 100);
    printf("Decomposition time:        %.3f ms (%.1f%%)\n", metrics.decomposition_time_ms, 
           (metrics.decomposition_time_ms / metrics.total_time_ms) * 100);
    printf("GMP to CUDA conversion:    %.3f ms (%.1f%%)\n", metrics.conversion_time_ms, 
           (metrics.conversion_time_ms / metrics.total_time_ms) * 100);
    printf("Host to device transfer:   %.3f ms (%.1f%%)\n", metrics.memory_htod_time_ms, 
           (metrics.memory_htod_time_ms / metrics.total_time_ms) * 100);
    printf("Kernel execution time:     %.3f ms (%.1f%%)\n", metrics.kernel_time_ms, 
           (metrics.kernel_time_ms / metrics.total_time_ms) * 100);
    printf("Device to host transfer:   %.3f ms (%.1f%%)\n", metrics.memory_dtoh_time_ms, 
           (metrics.memory_dtoh_time_ms / metrics.total_time_ms) * 100);
    printf("Result verification time:  %.3f ms (%.1f%%)\n", metrics.verification_time_ms, 
           (metrics.verification_time_ms / metrics.total_time_ms) * 100);
    printf("===================================\n");
}

bool miller_rabin_par(mpz_t n, unsigned int iterations, PerformanceMetrics* metrics = nullptr) {
    chrono::high_resolution_clock::time_point start_total, end_total;
    chrono::high_resolution_clock::time_point start_temp, end_temp;
    
    if (metrics) {
        metrics->total_time_ms = 0;
        metrics->base_generation_time_ms = 0;
        metrics->decomposition_time_ms = 0;
        metrics->conversion_time_ms = 0;
        metrics->memory_htod_time_ms = 0;
        metrics->kernel_time_ms = 0;
        metrics->memory_dtoh_time_ms = 0;
        metrics->verification_time_ms = 0;
        
        start_total = chrono::high_resolution_clock::now();
    }
    
    if (mpz_cmp_ui(n, 1) <= 0) return false;
    if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0) return true;
    if (mpz_divisible_ui_p(n, 2)) return false;
    
    int bits = mpz_sizeinbase(n, 2);
    int batch_size = iterations;
    
    if (bits > 8192) {
        batch_size = min(iterations, 32);
    }
    else if (bits > 4096) {
        batch_size = min(iterations, 64);
    }
    
    if (bits > MAX_BITS) {
        printf("Number has %d bits, which exceeds maximum of %d bits for CUDA implementation.\n", bits, MAX_BITS);
        printf("Falling back to CPU implementation.\n");
        
        if (metrics) start_temp = chrono::high_resolution_clock::now();
        bool cpu_result = miller_rabin_cpu(n, iterations);
        if (metrics) {
             end_temp = chrono::high_resolution_clock::now();
             metrics->kernel_time_ms = chrono::duration<double, milli>(end_temp - start_temp).count();
             
             end_total = chrono::high_resolution_clock::now();
             metrics->total_time_ms = chrono::duration<double, milli>(end_total - start_total).count();
             
             metrics->base_generation_time_ms = 0;
             metrics->decomposition_time_ms = 0;
             metrics->conversion_time_ms = 0;
             metrics->memory_htod_time_ms = 0;
             metrics->memory_dtoh_time_ms = 0;
             metrics->verification_time_ms = 0;
        }
        return cpu_result;
    }
    
    if (metrics) start_temp = chrono::high_resolution_clock::now();
    vector<mpz_t> bases = generate_random_bases(n, iterations);
    if (metrics) {
        end_temp = chrono::high_resolution_clock::now();
        metrics->base_generation_time_ms = chrono::duration<double, milli>(end_temp - start_temp).count();
    }
    
    if (metrics) start_temp = chrono::high_resolution_clock::now();
    MRDecomposition decomp = decompose(n);
    if (metrics) {
        end_temp = chrono::high_resolution_clock::now();
        metrics->decomposition_time_ms = chrono::duration<double, milli>(end_temp - start_temp).count();
    }
    
    bool is_prime = true;
    CudaBigInt cuda_n = mpz_to_cuda_bigint(n);
    CudaBigInt cuda_d = mpz_to_cuda_bigint(decomp.d);
    CudaBigInt cuda_r = mpz_to_cuda_bigint(decomp.r);
    
    if (metrics) start_temp = chrono::high_resolution_clock::now();
    
    cudaEvent_t kernel_start, kernel_end;
    if (metrics) {
        cudaEventCreate(&kernel_start);
        cudaEventCreate(&kernel_end);
    }
    
    for (unsigned int batch_start = 0; batch_start < iterations && is_prime; batch_start += batch_size) {
        unsigned int current_batch_size = min(batch_size, iterations - batch_start);
        
        CudaBigInt* host_cuda_bases = new CudaBigInt[current_batch_size];
        
        #pragma omp parallel for
        for (unsigned int i = 0; i < current_batch_size; i++) {
            host_cuda_bases[i] = mpz_to_cuda_bigint(bases[batch_start + i]);
        }
        
        CudaBigInt* device_bases;
        bool* device_results;
        bool* host_results = new bool[current_batch_size];
        
        cudaMalloc((void**)&device_bases, current_batch_size * sizeof(CudaBigInt));
        cudaMalloc((void**)&device_results, current_batch_size * sizeof(bool));
        
        cudaMemcpy(device_bases, host_cuda_bases, current_batch_size * sizeof(CudaBigInt), cudaMemcpyHostToDevice);
        
        float kernel_time = 0;
        
        if (metrics) {
            cudaEventRecord(kernel_start);
        }
        
        int num_blocks = (current_batch_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
        miller_rabin_kernel<<<num_blocks, BLOCK_SIZE>>>(cuda_n, device_bases, device_results, cuda_d, cuda_r, current_batch_size);
        cudaGetLastError();
        
        if (metrics) {
            cudaEventRecord(kernel_end);
            cudaEventSynchronize(kernel_end);
            cudaEventElapsedTime(&kernel_time, kernel_start, kernel_end);
            metrics->kernel_time_ms += kernel_time;
        }
        else {
            cudaDeviceSynchronize();
        }
        
        cudaMemcpy(host_results, device_results, current_batch_size * sizeof(bool), cudaMemcpyDeviceToHost);
        
        for (unsigned int i = 0; i < current_batch_size; i++) {
            if (!host_results[i]) {
                is_prime = false;
                break;
            }
        }
        
        cudaFree(device_bases);
        cudaFree(device_results);
        delete[] host_cuda_bases;
        delete[] host_results;
        
        if (!is_prime) break;
    }
    
    if (metrics) {
        cudaEventDestroy(kernel_start);
        cudaEventDestroy(kernel_end);
        
        end_temp = chrono::high_resolution_clock::now();
        metrics->conversion_time_ms = chrono::duration<double, milli>(end_temp - start_temp).count() - metrics->kernel_time_ms;
        // Since we've included memory operations in the conversion time, set these to 0
        metrics->memory_htod_time_ms = 0;
        metrics->memory_dtoh_time_ms = 0;
        metrics->verification_time_ms = 0;
    }
    
    for (unsigned int i = 0; i < iterations; i++) {
        mpz_clear(bases[i]);
    }
    
    mpz_clear(decomp.r);
    mpz_clear(decomp.d);
    
    if (metrics) {
        end_total = chrono::high_resolution_clock::now();
        metrics->total_time_ms = chrono::duration<double, milli>(end_total - start_total).count();
    }
    
    return is_prime;
}

int main(int argc, char** argv) {
    mpz_t number;
    unsigned int iterations;
    
    if (!parse_arguments(argc, argv, number, iterations)) {
        mpz_clear(number);
        return 1;
    }
    
    cout << "Testing if ";
    print_mpz(number);
    cout << " is prime using " << iterations << " parallel iterations..." << endl;
    cout << "Using CPU (OpenMP) for base generation and GPU (CUDA) for computation" << endl;
    
    PerformanceMetrics metrics;
    bool is_prime = miller_rabin_par(number, iterations, &metrics);
    
    cout << "Total execution time: " << metrics.total_time_ms << " ms" << endl;
    print_result(number, is_prime, metrics.total_time_ms);
    print_performance_metrics(metrics);
    
    mpz_clear(number);
    
    return 0;
} 