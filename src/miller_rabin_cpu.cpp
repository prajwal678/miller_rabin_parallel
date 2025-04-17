#include "miller_rabin.h"
#include "utils.h"
#include <iostream>

// Decompose n-1 as 2^r * d where d is odd
MRDecomposition decompose(u64 n) {
    MRDecomposition result;
    result.d = n - 1;
    result.r = 0;
    
    while (result.d % 2 == 0) {
        result.d /= 2;
        result.r++;
    }
    
    return result;
}

// Modular exponentiation: computes (base^exponent) % modulus
u64 mod_pow_cpu(u64 base, u64 exponent, u64 modulus) {
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

// Single Miller-Rabin test with base a
bool miller_rabin_test_cpu(u64 n, u64 a, MRDecomposition decomp) {
    if (n == 2 || n == 3) return true;
    if (n <= 1 || n % 2 == 0) return false;
    
    u64 d = decomp.d;
    u64 r = decomp.r;
    
    // Compute a^d % n
    u64 x = mod_pow_cpu(a, d, n);
    
    // If a^d % n == 1 or a^d % n == n-1, n passes the test
    if (x == 1 || x == n - 1) return true;
    
    // Square x, r-1 times
    for (u64 i = 0; i < r - 1; i++) {
        x = mod_pow_cpu(x, 2, n);
        if (x == n - 1) return true;
    }
    
    // If we get here, n is definitely composite
    return false;
}

// Full Miller-Rabin test with multiple iterations
bool miller_rabin_cpu(u64 n, unsigned int iterations) {
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0) return false;
    
    // Decompose n-1 = 2^r * d
    MRDecomposition decomp = decompose(n);
    
    // Generate random bases
    std::vector<u64> bases = generate_random_bases(n, iterations);
    
    // Test with each base
    for (u64 a : bases) {
        if (!miller_rabin_test_cpu(n, a, decomp)) {
            return false;  // Definitely composite
        }
    }
    
    return true;  // Probably prime
}

int main(int argc, char** argv) {
    u64 number;
    unsigned int iterations;
    
    if (!parse_arguments(argc, argv, number, iterations)) {
        return 1;
    }
    
    std::cout << "Testing if " << number << " is prime using " 
              << iterations << " iterations..." << std::endl;
    
    Timer timer("CPU Miller-Rabin test");
    bool is_prime = miller_rabin_cpu(number, iterations);
    double elapsed = timer.elapsed_ms();
    
    print_result(number, is_prime, elapsed);
    
    return 0;
} 