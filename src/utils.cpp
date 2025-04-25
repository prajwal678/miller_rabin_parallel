#include "utils.h"
#include <iostream>
#include <chrono>
#include <random>
#include <omp.h>
#include <gmp.h>

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <number_to_test> [iterations]" << std::endl;
    std::cerr << "  number_to_test: The number to test for primality" << std::endl;
    std::cerr << "  iterations: Number of iterations for the Miller-Rabin test (default: " 
              << DEFAULT_ITERATIONS << ")" << std::endl;
}

bool parse_arguments(int argc, char** argv, mpz_t number, unsigned int& iterations) {
    if (argc < 2 || argc > 3) {
        print_usage(argv[0]);
        return false;
    }

    try {
        mpz_init(number);
        if (mpz_set_str(number, argv[1], 10) != 0) {
            std::cerr << "Error: Invalid number format" << std::endl;
            return false;
        }
        
        if (mpz_cmp_ui(number, 2) < 0) {
            std::cerr << "Error: Number must be at least 2" << std::endl;
            return false;
        }

        if (argc == 3) {
            iterations = std::stoi(argv[2]);
            if (iterations < 1) {
                std::cerr << "Error: Iterations must be at least 1" << std::endl;
                return false;
            }
        } else {
            iterations = DEFAULT_ITERATIONS;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error parsing arguments: " << e.what() << std::endl;
        return false;
    }

    return true;
}

std::vector<mpz_t> generate_random_bases(mpz_t n, unsigned int count) {
    std::vector<mpz_t> bases(count);
    
    for (unsigned int i = 0; i < count; i++) {
        mpz_init(bases[i]);
    }
    
    #pragma omp parallel
    {
        gmp_randstate_t state;
        unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count() + omp_get_thread_num();
        gmp_randinit_default(state);
        gmp_randseed_ui(state, seed);
        
        mpz_t max_base;
        mpz_init(max_base);
        mpz_sub_ui(max_base, n, 3);  // n-3
        
        #pragma omp for
        for (unsigned int i = 0; i < count; i++) {
            mpz_urandomm(bases[i], state, max_base);
            mpz_add_ui(bases[i], bases[i], 2);  // shift range to [2, n-2]
        }
        
        mpz_clear(max_base);
        gmp_randclear(state);
    }
    
    return bases;
}

Timer::Timer(const std::string& name) : operation_name(name) {
    start_time = std::chrono::high_resolution_clock::now();
}

Timer::~Timer() {
    double duration = elapsed_ms();
    std::cout << operation_name << " took " << duration << " ms" << std::endl;
}

double Timer::elapsed_ms() {
    auto end_time = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(end_time - start_time).count();
}

void print_mpz(mpz_t num) {
    char* str = mpz_get_str(NULL, 10, num);
    std::cout << str;
    free(str);
}

void print_result(mpz_t number, bool is_prime, double elapsed_ms) {
    print_mpz(number);
    std::cout << " is " << (is_prime ? "probably prime" : "composite") 
              << " (determined in " << elapsed_ms << " ms)" << std::endl;
} 