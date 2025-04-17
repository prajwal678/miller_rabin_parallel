#ifndef MILLER_RABIN_H
#define MILLER_RABIN_H

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <omp.h>
#include <stdint.h>

#define DEFAULT_ITERATIONS 10

typedef uint64_t u64;
typedef unsigned long long ull;

// Structure to hold decomposition of n-1 as 2^r * d
struct MRDecomposition {
    u64 r;
    u64 d;
};

bool miller_rabin_cpu(u64 n, unsigned int iterations);
bool miller_rabin_test_cpu(u64 n, u64 a, MRDecomposition decomp);
MRDecomposition decompose(u64 n);
u64 mod_pow_cpu(u64 base, u64 exponent, u64 modulus);

bool miller_rabin_parallel(u64 n, unsigned int iterations);

std::vector<u64> generate_random_bases(u64 n, unsigned int count);

#endif // MILLER_RABIN_H 