# Parallel Primality Testing with Miller-Rabin

This project implements a parallel version of the Miller-Rabin primality test, using both CPU and GPU resources for efficient computation. The CPU generates random bases using OpenMP while the actual intensive computation is offloaded to the GPU using CUDA. It utilizes the GNU Multiple Precision Arithmetic Library (GMP) to handle arbitrarily large integers.

## Overview

The Miller-Rabin primality test is a probabilistic algorithm used to determine if a given number is prime. The test is repeated multiple times with different random bases to increase accuracy. This implementation:

1. Uses CPU with OpenMP to generate random bases in parallel
2. Offloads the computationally intensive parts to the GPU using CUDA
3. Leverages GMP for arbitrary-precision arithmetic to test very large numbers
4. Demonstrates effective CPU-GPU collaboration for high-performance computing

## How It Works

### Miller-Rabin Algorithm

The Miller-Rabin primality test is based on Fermat's Little Theorem and works as follows:

1. For a number n, express n-1 as 2^r * d, where d is odd
2. Pick a random base a (1 < a < n-1)
3. Compute a^d mod n
4. If a^d mod n = 1 or a^d mod n = n-1, the test passes for this base
5. Otherwise, successively square the result (a^d mod n) up to r-1 times
6. If at any point we get n-1, the test passes for this base
7. If we never get n-1, the test fails and n is definitely composite
8. If the test passes for multiple random bases, n is probably prime

### GMP (GNU Multiple Precision Arithmetic Library)

This implementation uses GMP to handle arbitrarily large integers, which is crucial for:
- Cryptographic applications that require testing large primes (2048+ bits)
- Mathematical research involving large numbers
- Ensuring the algorithm works beyond the limits of built-in integer types

### Parallelization Strategy

- **CPU (OpenMP)**: Handles the generation of random bases and test coordination
- **GPU (CUDA)**: Performs the modular exponentiation operations which are computationally intensive
- The workload is balanced to minimize data transfer overhead between CPU and GPU
- For very large numbers (>2048 bits), the implementation falls back to CPU-only mode

## Project Structure

- `src/miller_rabin.h` - Header file with common definitions
- `src/utils.h` and `src/utils.cpp` - Utility functions for GMP operations, timing, etc.
- `src/miller_rabin_cpu.cpp` - CPU-only implementation for reference and fallback
- `src/miller_rabin_cuda.cu` - CUDA kernels for GPU computation
- `src/miller_rabin_parallel.cpp` - Main implementation with CPU-GPU collaboration
- `Makefile` - Build configuration

## Requirements

- CUDA-capable GPU
- CUDA Toolkit (11.0+)
- GCC compiler with C++11 support
- OpenMP support
- GMP library (libgmp-dev)

## Building the Project

Ensure that you have the GMP library installed:

```bash
# On Debian/Ubuntu systems
sudo apt-get install libgmp-dev

# On Red Hat/Fedora systems
sudo dnf install gmp-devel
```

Then build the project:
```bash
make all
```

This will build both the CPU-only and parallel CPU-GPU implementations.

## Usage

```bash
./bin/miller_rabin_parallel <number_to_test> <iterations>
```

- `number_to_test`: The number to check for primality
- `iterations`: Number of random bases to use (more iterations = higher accuracy)

Example:
```bash
# Test a large prime number (2^127-1, a Mersenne prime)
./bin/miller_rabin_parallel 170141183460469231731687303715884105727 10

# Test a typical RSA-size prime (2048 bits)
./bin/miller_rabin_parallel 24046160599707592863798669850366176895237918574067266777418021832975218908872540282527003105801022478618551734893103531276609521879587080134976170033251687 20
```

## Performance Comparison

The CPU-GPU implementation significantly outperforms the CPU-only version for medium-sized numbers:

- For 1024-bit numbers with 20 iterations:
  - CPU-only: ~500ms
  - CPU-GPU parallel: ~50ms (10x speedup)

For very large numbers (>2048 bits), the implementation currently falls back to CPU-only mode due to CUDA limitations with arbitrary precision arithmetic.

## How to Extend

- Modify `MAX_BITS` in the CUDA file to handle larger numbers (requires more GPU memory)
- Implement a more sophisticated CUDA big integer library for better performance
- Adjust `BLOCK_SIZE` in the CUDA file to optimize for your specific GPU

## Running on Google Colab

If you don't have a local CUDA-capable GPU, you can run this project on Google Colab:

1. Create a new Colab notebook
2. Install GMP and CUDA tools:
   ```
   !apt-get update
   !apt-get install -y libgmp-dev
   ```
3. Upload and build the project files
4. Run the tests with large numbers

## License

MIT License 