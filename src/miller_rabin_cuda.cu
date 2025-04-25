#include "miller_rabin.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <gmp.h>
#include <curand_kernel.h>

#define BLOCK_SIZE 256
#define MAX_BITS 2048

#define cudaCheckError(ans) { cudaAssert((ans), __FILE__, __LINE__); }
inline void cudaAssert(cudaError_t code, const char *file, int line) {
    if (code != cudaSuccess) {
        fprintf(stderr, "CUDA Error: %s at %s:%d\n", cudaGetErrorString(code), file, line);
        exit(code);
    }
}

typedef struct {
    unsigned int data[MAX_BITS/32 + 1];
    int bits;
} CudaBigInt;

CudaBigInt mpz_to_cuda_bigint(mpz_t num) {
    CudaBigInt result;
    result.bits = mpz_sizeinbase(num, 2);
    memset(result.data, 0, sizeof(result.data));
    size_t count = (result.bits + 31) / 32;
    mpz_export(result.data, &count, -1, sizeof(unsigned int), 0, 0, num);
    
    return result;
}

void cuda_bigint_to_mpz(CudaBigInt* cbi, mpz_t result) {
    size_t count = (cbi->bits + 31) / 32;
    mpz_import(result, count, -1, sizeof(unsigned int), 0, 0, cbi->data);
}


// Device function for modular multiplication on CUDA
__device__ void cuda_mod_mul(CudaBigInt* result, const CudaBigInt* a, const CudaBigInt* b, const CudaBigInt* modulus) {
    // Reset result
    memset(result->data, 0, sizeof(unsigned int) * (MAX_BITS/32 + 1));
    result->bits = 0;
    
    // Compute a * b
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
    
    // Compute result mod modulus using division
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
        return; // Result is already smaller than modulus
    }
    
    // Perform modulo operation
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
    
    // Calculate result bits
    for (int i = MAX_BITS/32; i >= 0; i--) {
        if (result->data[i] > 0) {
            result->bits = i * 32 + 32;
            while ((result->data[i] & (1 << (result->bits % 32 - 1))) == 0) {
                result->bits--;
            }
            break;
        }
    }
    if (result->bits == 0) result->bits = 1; // At minimum 1 bit for zero
}

// Device function for modular exponentiation on CUDA
__device__ void cuda_mod_pow(CudaBigInt* result, const CudaBigInt* base, const CudaBigInt* exponent, const CudaBigInt* modulus) {
    // Initialize result to 1
    memset(result->data, 0, sizeof(unsigned int) * (MAX_BITS/32 + 1));
    result->data[0] = 1;
    result->bits = 1;
    
    // Handle base case of exponent = 0
    bool is_zero = true;
    for (int i = 0; i < (exponent->bits + 31) / 32; i++) {
        if (exponent->data[i] != 0) {
            is_zero = false;
            break;
        }
    }
    if (is_zero) return;
    
    // Create a copy of base for squaring
    CudaBigInt base_copy;
    memcpy(&base_copy, base, sizeof(CudaBigInt));
    
    // Square-and-multiply algorithm
    for (int bit_pos = exponent->bits - 1; bit_pos >= 0; bit_pos--) {
        // Square result
        CudaBigInt squared;
        cuda_mod_mul(&squared, result, result, modulus);
        memcpy(result, &squared, sizeof(CudaBigInt));
        
        // Multiply by base if current bit is set
        int word_idx = bit_pos / 32;
        int bit_idx = bit_pos % 32;
        if (exponent->data[word_idx] & (1 << bit_idx)) {
            CudaBigInt product;
            cuda_mod_mul(&product, result, &base_copy, modulus);
            memcpy(result, &product, sizeof(CudaBigInt));
        }
    }
}

// CUDA kernel for Miller-Rabin test with one base per thread
__global__ void miller_rabin_kernel(CudaBigInt n, CudaBigInt* bases, bool* results, 
                                   CudaBigInt d, CudaBigInt r, unsigned int iterations) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < iterations) {
        CudaBigInt a = bases[idx];
        bool is_prime = true;
        
        // Miller-Rabin test with a single base
        CudaBigInt x;
        cuda_mod_pow(&x, &a, &d, &n);
        
        // Check if x is 1 or n-1
        bool is_one = true;
        for (int i = 1; i < (x.bits + 31) / 32; i++) {
            if (x.data[i] != 0) {
                is_one = false;
                break;
            }
        }
        if (x.data[0] != 1) is_one = false;
        
        // Calculate n-1 for comparison
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
            // Square x for r-1 iterations and check if it equals n-1
            bool found_minus_1 = false;
            for (unsigned int j = 0; j < r.data[0] - 1; j++) {
                // Compute x = x^2 mod n
                CudaBigInt squared;
                cuda_mod_mul(&squared, &x, &x, &n);
                memcpy(&x, &squared, sizeof(CudaBigInt));
                
                // Check if x is 1 (composite)
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
                
                // Check if x is n-1
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
            
            // If we never found x = n-1, n is composite
            if (!found_minus_1) {
                is_prime = false;
            }
        }
        
        results[idx] = is_prime;
    }
}

// Host function to perform Miller-Rabin test on GPU with GMP
bool miller_rabin_parallel(mpz_t n, unsigned int iterations) {
    // Handle small cases on the host
    if (mpz_cmp_ui(n, 1) <= 0) return false;
    if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0) return true;
    if (mpz_divisible_ui_p(n, 2)) return false;
    
    // For very large numbers, we'll use the CPU implementation
    // as a full CUDA implementation for arbitrary precision is complex
    int bits = mpz_sizeinbase(n, 2);
    if (bits > MAX_BITS) {
        printf("Number has %d bits, which exceeds maximum of %d bits for CUDA implementation.\n", 
               bits, MAX_BITS);
        printf("Falling back to CPU implementation.\n");
        return miller_rabin_cpu(n, iterations);
    }
    
    // Generate random bases on the CPU using OpenMP
    std::vector<mpz_t> bases = generate_random_bases(n, iterations);
    
    // Decompose n-1 = 2^r * d
    MRDecomposition decomp = decompose(n);
    
    // Convert GMP numbers to CUDA-compatible format
    CudaBigInt cuda_n = mpz_to_cuda_bigint(n);
    CudaBigInt cuda_d = mpz_to_cuda_bigint(decomp.d);
    CudaBigInt cuda_r = mpz_to_cuda_bigint(decomp.r);
    
    // Convert bases to CUDA format
    CudaBigInt* host_cuda_bases = new CudaBigInt[iterations];
    for (unsigned int i = 0; i < iterations; i++) {
        host_cuda_bases[i] = mpz_to_cuda_bigint(bases[i]);
    }
    
    // Allocate device memory
    CudaBigInt* device_bases;
    bool* device_results;
    
    cudaCheckError(cudaMalloc((void**)&device_bases, iterations * sizeof(CudaBigInt)));
    cudaCheckError(cudaMalloc((void**)&device_results, iterations * sizeof(bool)));
    
    // Copy data to device
    cudaCheckError(cudaMemcpy(device_bases, host_cuda_bases, iterations * sizeof(CudaBigInt), cudaMemcpyHostToDevice));
    
    // Launch kernel
    int num_blocks = (iterations + BLOCK_SIZE - 1) / BLOCK_SIZE;
    miller_rabin_kernel<<<num_blocks, BLOCK_SIZE>>>(cuda_n, device_bases, device_results, cuda_d, cuda_r, iterations);
    cudaCheckError(cudaGetLastError());
    cudaCheckError(cudaDeviceSynchronize());
    
    // Copy results back to host
    bool* host_results = new bool[iterations];
    cudaCheckError(cudaMemcpy(host_results, device_results, iterations * sizeof(bool), cudaMemcpyDeviceToHost));
    
    // Check if all tests passed
    bool is_prime = true;
    for (unsigned int i = 0; i < iterations; i++) {
        if (!host_results[i]) {
            is_prime = false;
            break;
        }
    }
    
    // Clean up
    cudaCheckError(cudaFree(device_bases));
    cudaCheckError(cudaFree(device_results));
    
    for (unsigned int i = 0; i < iterations; i++) {
        mpz_clear(bases[i]);
    }
    
    mpz_clear(decomp.r);
    mpz_clear(decomp.d);
    
    delete[] host_cuda_bases;
    delete[] host_results;
    
    return is_prime;
} 