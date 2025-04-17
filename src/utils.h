#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include "miller_rabin.h"

void print_usage(const char* program_name);

bool parse_arguments(int argc, char** argv, u64& number, unsigned int& iterations);

// Timer functions for benchmarking
class Timer {
private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::string operation_name;

public:
    Timer(const std::string& name);
    ~Timer();
    double elapsed_ms();
};

void print_result(u64 number, bool is_prime, double elapsed_ms);

#endif // UTILS_H 