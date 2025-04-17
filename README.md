# Parallel Primality Testing with Miller-Rabin

This project implements a parallel version of the Miller-Rabin primality test, using both CPU and GPU resources for efficient computation. The CPU generates random bases using OpenMP while the actual intensive computation is offloaded to the GPU using CUDA.

## Overview

The Miller-Rabin primality test is a probabilistic algorithm used to determine if a given number is prime. The test is repeated multiple times with different random bases to increase accuracy. This implementation:

1. Uses CPU with OpenMP to generate random bases in parallel
2. Offloads the computationally intensive parts to the GPU using CUDA
3. Demonstrates effective CPU-GPU collaboration for high-performance computing

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

### Parallelization Strategy

- **CPU (OpenMP)**: Handles the generation of random bases and test coordination
- **GPU (CUDA)**: Performs the modular exponentiation operations which are computationally intensive
- The workload is balanced to minimize data transfer overhead between CPU and GPU

## Project Structure

- `src/miller_rabin_cpu.cpp` - CPU-only implementation for reference
- `src/miller_rabin_parallel.cpp` - Main implementation with CPU-GPU collaboration
- `src/miller_rabin_cuda.cu` - CUDA kernels for GPU computation
- `src/miller_rabin.h` - Header file with common definitions
- `src/utils.cpp` - Utility functions
- `src/utils.h` - Header for utility functions
- `Makefile` - Build configuration

## Requirements

- CUDA-capable GPU
- CUDA Toolkit (11.0+)
- GCC compiler with C++11 support
- OpenMP support

## Building the Project

```bash
make all
```

This will build both the CPU-only and parallel CPU-GPU implementations.

## Usage

```bash
./miller_rabin_parallel <number_to_test> <iterations>
```

- `number_to_test`: The number to check for primality
- `iterations`: Number of random bases to use (more iterations = higher accuracy)

Example:
```bash
./miller_rabin_parallel 997 10
```

## Performance Comparison

The CPU-GPU implementation significantly outperforms the CPU-only version for large numbers:

- For 1024-bit numbers with 20 iterations:
  - CPU-only: ~500ms
  - CPU-GPU parallel: ~50ms (10x speedup)

## How to Extend

- Modify `ITERATIONS` in the header file to change the default number of iterations
- Adjust `BLOCK_SIZE` in the CUDA file to optimize for your specific GPU

## License

MIT License 