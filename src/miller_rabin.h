#ifndef MILLER_RABIN_H
#define MILLER_RABIN_H

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <omp.h>
#include <gmp.h>

#define DEFAULT_ITERATIONS 50


// decomposition of n-1 as 2^r * d
struct MRDecomposition {
    mpz_t r;
    mpz_t d;
};

MRDecomposition decompose(mpz_t n);
bool miller_rabin_test_cpu(mpz_t n, mpz_t a, MRDecomposition decomp);

std::vector<mpz_t> generate_random_bases(mpz_t n, unsigned int count); // Parallel (OpenMP)
std::vector<mpz_t> generate_random_bases_seq(mpz_t n, unsigned int count); // Sequential

bool miller_rabin_seq(mpz_t n, unsigned int iterations);
bool miller_rabin_cpu(mpz_t n, unsigned int iterations);
bool miller_rabin_par(mpz_t n, unsigned int iterations);

#endif // MILLER_RABIN_H 