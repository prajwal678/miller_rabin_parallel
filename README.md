# Miller-Rabin Primality Test: Sequential vs Parallel (CUDA)

This project implements and compares two versions of the Miller-Rabin primality test:
1.  **Sequential:** A standard C++ implementation running entirely on the CPU.
2.  **Parallel:** A hybrid implementation using OpenMP for parallel random base generation on the CPU and CUDA for massively parallel primality test computations on the GPU.

## Features

*   Arbitrary-precision integer arithmetic using GMP.
*   High-quality random number generation for bases using `/dev/urandom`.
*   CUDA implementation for GPU acceleration.
*   OpenMP for parallelizing CPU-bound tasks (base generation).
*   Benchmarking script to compare performance and calculate speedup.
*   Detailed performance metrics for the parallel implementation.

## Prerequisites

*   A C++11 compliant compiler (like g++)
*   NVIDIA CUDA Toolkit (if running the parallel version)
*   GNU Multiple Precision Arithmetic Library (GMP)
*   Standard Unix utilities: `make`, `bc`, `xxd`, `dd`

## Setup and Execution on Google Colab (Python 3 Environment)

1.  **Open a Colab Notebook:** Start a new Colab notebook.

2.  **Select GPU Runtime (Required for Parallel Version):**
    *   Go to `Runtime` -> `Change runtime type`.
    *   Select `GPU` from the `Hardware accelerator` dropdown menu.
    *   Click `Save`.

3.  **Clone the Repository:**
    ```python
    !git clone <your-repository-url> 
    %cd <your-repository-directory>
    ```
    *(Replace `<your-repository-url>` and `<your-repository-directory>`)*

4.  **Install Dependencies:**
    ```python
    !apt-get update
    !apt-get install -y build-essential libgmp-dev bc xxd 
    # Install CUDA toolkit if needed (might already be present on GPU runtimes)
    !apt-get install -y nvidia-cuda-toolkit 
    ```

5.  **Compile the Code:**
    ```python
    !make clean
    !make
    ```
    *Note: You might see warnings about file modification times, which can usually be ignored on Colab.*

6.  **Make the Benchmark Script Executable:**
    ```python
    !chmod +x benchmark.sh
    ```

7.  **Run the Benchmark:**
    ```python
    !./benchmark.sh
    ```
    This will run both the sequential and parallel versions, testing primality for numbers of different bit sizes (512, 1024, 2048, 4096 by default) and outputting performance results.

8.  **View the Results:**
    ```python
    !cat results/summary.md
    ```
    This command will display the summary report comparing the execution times and speedup achieved by the parallel version.

## Project Structure

```
.
├── Makefile          # Builds the project
├── benchmark.sh      # Runs tests and generates performance reports
├── README.md         # This file
├── results/          # Directory for benchmark output (created by benchmark.sh)
└── src/
    ├── miller_rabin.h        # Header file with common declarations
    ├── miller_rabin_cpu.cpp  # CPU-specific functions (decompose, test_cpu)
    ├── miller_rabin_par.cu   # Parallel implementation (OpenMP + CUDA)
    ├── miller_rabin_seq.cpp  # Sequential implementation
    ├── utils.h               # Header for utility functions
    └── utils.cpp             # Utility functions (arg parsing, random base generation, timer)
```

## How It Works

*   **Sequential (`miller_rabin_seq`):** Takes a number and iteration count. It generates random bases sequentially using `/dev/urandom` and performs the Miller-Rabin test entirely on the CPU using GMP for calculations.
*   **Parallel (`miller_rabin_par`):** Takes a number and iteration count. It uses OpenMP to generate random bases in parallel on CPU cores (reading from `/dev/urandom` for each). The computationally intensive part (modular exponentiation) is offloaded to the GPU using a CUDA kernel. Data is transferred between CPU and GPU memory as needed. It falls back to a CPU implementation for numbers exceeding the `MAX_BITS` limit defined in the code.
*   **Benchmarking (`benchmark.sh`):** Generates a large random number for each specified bit size. It then runs both the sequential and parallel executables with this *same* number and measures the execution time. It repeats this process multiple times and calculates the average times and speedup.