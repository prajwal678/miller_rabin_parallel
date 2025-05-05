#include "miller_rabin.h"
#include "utils.h"
#include <iostream>
#include <gmp.h>
#include <chrono>
#include <stdexcept>
#include <random>

using namespace std;


// Single Miller-Rabin test with base a
bool miller_rabin_test_seq(mpz_t n, mpz_t a, MRDecomposition decomp) {
    if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0) return true;
    if (mpz_cmp_ui(n, 1) <= 0 || mpz_divisible_ui_p(n, 2)) return false;
    
    mpz_t x, n_minus_1;
    mpz_init(x);
    mpz_init(n_minus_1);
    mpz_sub_ui(n_minus_1, n, 1);
    
    mpz_powm(x, a, decomp.d, n);
    
    // if a^d % n == 1 or a^d % n == n-1, n passes the test
    if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_minus_1) == 0) {
        mpz_clear(x);
        mpz_clear(n_minus_1);
        return true;
    }
    
    mpz_t r_minus_1, i;
    mpz_init(r_minus_1);
    mpz_init_set_ui(i, 0);
    mpz_sub_ui(r_minus_1, decomp.r, 1);
    
    while (mpz_cmp(i, r_minus_1) < 0) {
        mpz_powm_ui(x, x, 2, n);  // x = x^2 % n
        if (mpz_cmp(x, n_minus_1) == 0) {
            mpz_clear(x);
            mpz_clear(n_minus_1);
            mpz_clear(r_minus_1);
            mpz_clear(i);
            return true;
        }
        mpz_add_ui(i, i, 1);
    }
    
    mpz_clear(x);
    mpz_clear(n_minus_1);
    mpz_clear(r_minus_1);
    mpz_clear(i);
    
    return false;
}

bool miller_rabin_seq(mpz_t n, unsigned int iterations) {
    if (mpz_cmp_ui(n, 1) <= 0) return false;
    if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0) return true;
    if (mpz_divisible_ui_p(n, 2)) return false;
    
    MRDecomposition decomp = decompose(n);
    vector<mpz_t> bases = generate_random_bases_seq(n, iterations);
    
    bool probable_prime = true;
    for (unsigned int i = 0; i < iterations; i++) {
        if (!miller_rabin_test_seq(n, bases[i], decomp)) {
            probable_prime = false;
            break;
        }
    }
    
    for (unsigned int i = 0; i < iterations; i++) {
        mpz_clear(bases[i]);
    }
    mpz_clear(decomp.r);
    mpz_clear(decomp.d);
    
    return probable_prime;
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
    cout << " is prime using " << iterations << " sequential iterations..." << endl;
    
    auto start_time = chrono::high_resolution_clock::now();
    bool is_prime = miller_rabin_seq(number, iterations);
    auto end_time = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double, milli>(end_time - start_time).count();
    
    cout << "Total execution time: " << elapsed << " ms" << endl;
    print_result(number, is_prime, elapsed);
    
    mpz_clear(number);
    
    return 0;
} 