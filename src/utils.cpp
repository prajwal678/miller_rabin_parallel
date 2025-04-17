#include "utils.h"
#include <iostream>
#include <chrono>
#include <random>
#include <omp.h>

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <number_to_test> [iterations]" << std::endl;
    std::cerr << "  number_to_test: The number to test for primality" << std::endl;
    std::cerr << "  iterations: Number of iterations for the Miller-Rabin test (default: " 
              << DEFAULT_ITERATIONS << ")" << std::endl;
}

bool parse_arguments(int argc, char** argv, u64& number, unsigned int& iterations) {
    if (argc < 2 || argc > 3) {
        print_usage(argv[0]);
        return false;
    }

    try {
        number = stoull(argv[1]);
        if (number < 2) {
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

// Generate random bases for Miller-Rabin test
std::vector<u64> generate_random_bases(u64 n, unsigned int count) {
    std::vector<u64> bases(count);
    
    #pragma omp parallel
    {
        // each thread gen a random number (CPU side)
        unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count() + omp_get_thread_num();
        std::mt19937_64 generator(seed);
        std::uniform_int_distribution<u64> distribution(2, n - 2);
        
        #pragma omp for
        for (unsigned int i = 0; i < count; i++) {
            bases[i] = distribution(generator);
        }
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

void print_result(u64 number, bool is_prime, double elapsed_ms) {
    std::cout << number << " is " << (is_prime ? "probably prime" : "composite") 
              << " (determined in " << elapsed_ms << " ms)" << std::endl;
} 