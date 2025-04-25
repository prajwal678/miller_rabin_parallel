#include "miller_rabin.h"
#include "utils.h"
#include <iostream>
#include <gmp.h>

// n-1 as 2^r * d where d is odd
MRDecomposition decompose(mpz_t n) {
    MRDecomposition result;
    mpz_init(result.r);
    mpz_init(result.d);
    
    mpz_t n_minus_1;
    mpz_init(n_minus_1);
    mpz_sub_ui(n_minus_1, n, 1);
    
    // initialize d = n-1
    mpz_set(result.d, n_minus_1);
    mpz_set_ui(result.r, 0);
    
    // largest power of 2 that divides n-1
    mpz_t remainder;
    mpz_init(remainder);
    
    while (true) {
        mpz_mod_ui(remainder, result.d, 2);
        if (mpz_cmp_ui(remainder, 0) != 0) {
            break;  // d is odd
        }
        mpz_fdiv_q_ui(result.d, result.d, 2);  // d = d/2
        mpz_add_ui(result.r, result.r, 1);     // r = r+1
    }
    
    mpz_clear(n_minus_1);
    mpz_clear(remainder);
    
    return result;
}

// single Miller-Rabin test with base a
bool miller_rabin_test_cpu(mpz_t n, mpz_t a, MRDecomposition decomp) {
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

// full Miller-Rabin test with multiple iterations
bool miller_rabin_cpu(mpz_t n, unsigned int iterations) {
    if (mpz_cmp_ui(n, 1) <= 0) return false;
    if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0) return true;
    if (mpz_divisible_ui_p(n, 2)) return false;
    
    MRDecomposition decomp = decompose(n);
    std::vector<mpz_t> bases = generate_random_bases(n, iterations);
    
    for (unsigned int i = 0; i < iterations; i++) {
        if (!miller_rabin_test_cpu(n, bases[i], decomp)) {
            for (unsigned int j = 0; j < iterations; j++) {
                mpz_clear(bases[j]);
            }
            mpz_clear(decomp.r);
            mpz_clear(decomp.d);
            return false;
        }
    }
    
    for (unsigned int i = 0; i < iterations; i++) {
        mpz_clear(bases[i]);
    }
    mpz_clear(decomp.r);
    mpz_clear(decomp.d);
    
    return true;  // PROBABLY primy BRUV, dont take true as yes
}

int main(int argc, char** argv) {
    mpz_t number;
    unsigned int iterations;
    
    if (!parse_arguments(argc, argv, number, iterations)) {
        mpz_clear(number);
        return 1;
    }
    
    std::cout << "Testing if ";
    print_mpz(number);
    std::cout << " is prime using " << iterations << " iterations..." << std::endl;
    
    Timer timer("CPU Miller-Rabin test");
    bool is_prime = miller_rabin_cpu(number, iterations);
    double elapsed = timer.elapsed_ms();
    
    print_result(number, is_prime, elapsed);
    
    mpz_clear(number);
    
    return 0;
} 