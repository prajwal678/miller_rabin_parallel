#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <gmp.h>
#include "miller_rabin.h"

using namespace std;


class Timer {
private:
    chrono::high_resolution_clock::time_point start_time;
    string operation_name;

public:
    Timer(const string& name);
    ~Timer();
    double elapsed_ms();
};

void print_usage(const char* program_name);
bool parse_arguments(int argc, char** argv, mpz_t number, unsigned int& iterations);
void print_mpz(mpz_t num);
void print_result(mpz_t number, bool is_prime, double elapsed_ms);

#endif // UTILS_H 