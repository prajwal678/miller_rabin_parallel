#include "utils.h"
#include <iostream>
#include <chrono>
#include <random>
#include <omp.h>
#include <gmp.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdexcept>
#include <string>

using namespace std;


void print_usage(const char* program_name) {
    cerr << "Usage: " << program_name << " <number_to_test> [iterations]" << endl;
    cerr << "  number_to_test: The number to test for primality (in base 10)" << endl;
    cerr << "  iterations: Number of iterations for the Miller-Rabin test (default: " 
              << DEFAULT_ITERATIONS << ")" << endl;
}

bool parse_arguments(int argc, char** argv, mpz_t number, unsigned int& iterations) {
    if (argc < 2 || argc > 3) {
        print_usage(argv[0]);
        return false;
    }

    try {
        mpz_init(number);
        if (mpz_set_str(number, argv[1], 10) != 0) {
            cerr << "Error: Invalid number format for <number_to_test>" << endl;
            mpz_clear(number);
            return false;
        }
        
        if (mpz_cmp_ui(number, 2) < 0) {
            cerr << "Error: Number must be at least 2" << endl;
            mpz_clear(number);
            return false;
        }

        if (argc == 3) {
            iterations = stoi(argv[2]);
            if (iterations < 1) {
                cerr << "Error: Iterations must be at least 1" << endl;
                mpz_clear(number);
                return false;
            }
        } else {
            iterations = DEFAULT_ITERATIONS;
        }
    } catch (const invalid_argument& e) {
        cerr << "Error: Invalid number format for iterations." << endl;
        mpz_clear(number);
        return false;
    } catch (const out_of_range& e) {
        cerr << "Error: Iterations value is out of range." << endl;
        mpz_clear(number);
        return false;
    } catch (const exception& e) {
        cerr << "Error parsing arguments: " << e.what() << endl;
        mpz_clear(number);
        return false;
    }

    return true;
}

void get_random_bytes(unsigned char* buffer, size_t nbytes) {
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd == -1) {
        cerr << "Failed to open /dev/urandom" << endl;
        exit(1);
    }
    
    size_t bytes_read = 0;
    while (bytes_read < nbytes) {
        ssize_t result = read(fd, buffer + bytes_read, nbytes - bytes_read);
        if (result <= 0) {
            cerr << "Error reading from /dev/urandom" << endl;
            close(fd);
            exit(1);
        }
        bytes_read += result;
    }
    
    close(fd);
}

vector<mpz_t> generate_random_bases(mpz_t n, unsigned int count) {
    vector<mpz_t> bases(count);
    
    for (unsigned int i = 0; i < count; i++) {
        mpz_init(bases[i]);
    }
    
    #pragma omp parallel
    {
        size_t nbits = mpz_sizeinbase(n, 2);
        size_t nbytes = (nbits + 7) / 8;
        unsigned char* random_data = new unsigned char[nbytes];
        
        mpz_t max_base;
        mpz_init(max_base);
        mpz_sub_ui(max_base, n, 3);
        
        #pragma omp for
        for (unsigned int i = 0; i < count; i++) {
            get_random_bytes(random_data, nbytes);
            
            mpz_t random_num;
            mpz_init(random_num);
            mpz_import(random_num, nbytes, 1, 1, 0, 0, random_data);
            
            mpz_mod(random_num, random_num, max_base);
            mpz_add_ui(bases[i], random_num, 2);
            mpz_clear(random_num);
        }
        
        mpz_clear(max_base);
        delete[] random_data;
    }
    
    return bases;
}

// seq version of random base generation that still uses /dev/urandom
vector<mpz_t> generate_random_bases_seq(mpz_t n, unsigned int count) {
    vector<mpz_t> bases(count);
    
    for (unsigned int i = 0; i < count; i++) {
        mpz_init(bases[i]);
    }
    
    size_t nbits = mpz_sizeinbase(n, 2);
    size_t nbytes = (nbits + 7) / 8;
    unsigned char* random_data = new unsigned char[nbytes];
    
    mpz_t max_base;
    mpz_init(max_base);
    mpz_sub_ui(max_base, n, 3);
    
    for (unsigned int i = 0; i < count; i++) {
        get_random_bytes(random_data, nbytes);
        
        mpz_t random_num;
        mpz_init(random_num);
        mpz_import(random_num, nbytes, 1, 1, 0, 0, random_data);
        
        mpz_mod(random_num, random_num, max_base);
        mpz_add_ui(bases[i], random_num, 2);
        mpz_clear(random_num);
    }
    
    mpz_clear(max_base);
    delete[] random_data;
    
    return bases;
}

Timer::Timer(const string& name) : operation_name(name) {
    start_time = chrono::high_resolution_clock::now();
}

Timer::~Timer() {
}

double Timer::elapsed_ms() {
    auto end_time = chrono::high_resolution_clock::now();
    return chrono::duration<double, milli>(end_time - start_time).count();
}

void print_mpz(mpz_t num) {
    char* str = mpz_get_str(NULL, 10, num);
    cout << str;
    free(str);
}

void print_result(mpz_t number, bool is_prime, double elapsed_ms) {
    print_mpz(number);
    cout << " is " << (is_prime ? "probably prime" : "composite") 
              << " (determined in " << elapsed_ms << " ms)" << endl;
}