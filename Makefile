CC = g++
NVCC = nvcc
CFLAGS = -std=c++11 -O3 -fopenmp
NVCCFLAGS = -std=c++11 -O3
LDFLAGS = -lcuda -lcudart

# Directories
SRC_DIR = src
BIN_DIR = bin

# Target executables
CPU_TARGET = $(BIN_DIR)/miller_rabin_cpu
PARALLEL_TARGET = $(BIN_DIR)/miller_rabin_parallel

# Source files
COMMON_SRC = $(SRC_DIR)/utils.cpp
CPU_SRC = $(SRC_DIR)/miller_rabin_cpu.cpp
CUDA_SRC = $(SRC_DIR)/miller_rabin_cuda.cu
PARALLEL_SRC = $(SRC_DIR)/miller_rabin_parallel.cpp

# Default target
all: dirs $(CPU_TARGET) $(PARALLEL_TARGET)

# Create bin directory if it doesn't exist
dirs:
	mkdir -p $(BIN_DIR)

# CPU-only version
$(CPU_TARGET): $(CPU_SRC) $(COMMON_SRC)
	$(CC) $(CFLAGS) -o $@ $^

# GPU version with CUDA
$(PARALLEL_TARGET): $(PARALLEL_SRC) $(CUDA_SRC) $(COMMON_SRC)
	$(NVCC) $(NVCCFLAGS) -o $@ $^ $(LDFLAGS)

# Clean
clean:
	rm -rf $(BIN_DIR)

.PHONY: all dirs clean 