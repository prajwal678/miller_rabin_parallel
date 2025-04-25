#include "miller_rabin.h"
#include "utils.h"
#include <iostream>
#include <gmp.h>

int main(int argc, char** argv) {
    mpz_t number;
    unsigned int iterations;
    
    if (!parse_arguments(argc, argv, number, iterations)) {
        mpz_clear(number);
        return 1;
    }
    
    std::cout << "Testing if ";
    print_mpz(number);
    std::cout << " is prime using " << iterations << " parallel iterations..." << std::endl;
    std::cout << "Using CPU (OpenMP) for base generation and GPU (CUDA) for computation" << std::endl;
    
    Timer timer("Parallel Miller-Rabin test");
    bool is_prime = miller_rabin_parallel(number, iterations);
    double elapsed = timer.elapsed_ms();
    
    print_result(number, is_prime, elapsed);
    
    mpz_clear(number);
    
    return 0;
} 