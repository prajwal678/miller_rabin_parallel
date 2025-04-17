#include "miller_rabin.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>

// CUDA kernel parameters
#define BLOCK_SIZE 256

// Error checking macro for CUDA calls
#define cudaCheckError(ans) { cudaAssert((ans), __FILE__, __LINE__); }
inline void cudaAssert(cudaError_t code, const char *file, int line) {
    if (code != cudaSuccess) {
        fprintf(stderr, "CUDA Error: %s at %s:%d\n", cudaGetErrorString(code), file, line);
        exit(code);
    }
}

// Device function for modular exponentiation on GPU
__device__ u64 mod_pow_gpu(u64 base, u64 exponent, u64 modulus) {
    if (modulus == 1) return 0;
    
    u64 result = 1;
    base = base % modulus;
    
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = (result * base) % modulus;
        }
        exponent >>= 1;
        base = (base * base) % modulus;
    }
    
    return result;
}

// Single Miller-Rabin test kernel for a specific base
__global__ void miller_rabin_test_kernel(u64 n, u64* bases, bool* results, u64 d, u64 r, unsigned int iterations) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < iterations) {
        u64 a = bases[idx];
        bool is_prime = true;
        
        // Handle small cases
        if (n <= 1 || n % 2 == 0 && n != 2) {
            is_prime = false;
        } else if (n <= 3) {
            is_prime = true;
        } else {
            // Compute a^d % n
            u64 x = mod_pow_gpu(a, d, n);
            
            // If a^d % n == 1 or a^d % n == n-1, n passes the test for this base
            if (x == 1 || x == n - 1) {
                is_prime = true;
            } else {
                // Square x, r-1 times
                bool found_witness = false;
                for (u64 i = 0; i < r - 1 && !found_witness; i++) {
                    x = mod_pow_gpu(x, 2, n);
                    if (x == n - 1) {
                        found_witness = true;
                        is_prime = true;
                    }
                }
                
                if (!found_witness) {
                    is_prime = false;
                }
            }
        }
        
        results[idx] = is_prime;
    }
}

// Host function to perform Miller-Rabin test on GPU
bool miller_rabin_parallel(u64 n, unsigned int iterations) {
    // Handle small cases on the host
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0) return false;
    
    // Decompose n-1 = 2^r * d on the host
    MRDecomposition decomp = decompose(n);
    
    // Generate random bases on the CPU using OpenMP
    std::vector<u64> host_bases = generate_random_bases(n, iterations);
    
    // Allocate device memory
    u64* device_bases;
    bool* device_results;
    
    cudaCheckError(cudaMalloc((void**)&device_bases, iterations * sizeof(u64)));
    cudaCheckError(cudaMalloc((void**)&device_results, iterations * sizeof(bool)));
    
    // Copy bases to device
    cudaCheckError(cudaMemcpy(device_bases, host_bases.data(), iterations * sizeof(u64), cudaMemcpyHostToDevice));
    
    // Calculate grid and block dimensions
    int num_blocks = (iterations + BLOCK_SIZE - 1) / BLOCK_SIZE;
    
    // Launch kernel
    miller_rabin_test_kernel<<<num_blocks, BLOCK_SIZE>>>(n, device_bases, device_results, decomp.d, decomp.r, iterations);
    cudaCheckError(cudaGetLastError());
    
    // Wait for kernel to finish
    cudaCheckError(cudaDeviceSynchronize());
    
    // Copy results back to host
    std::vector<bool> host_results(iterations);
    cudaCheckError(cudaMemcpy(host_results.data(), device_results, iterations * sizeof(bool), cudaMemcpyDeviceToHost));
    
    // Free device memory
    cudaCheckError(cudaFree(device_bases));
    cudaCheckError(cudaFree(device_results));
    
    // Check if all tests passed
    for (bool result : host_results) {
        if (!result) {
            return false;  // Definitely composite
        }
    }
    
    return true;  // Probably prime
} 