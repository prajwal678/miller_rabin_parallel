CC = g++
NVCC = nvcc
CFLAGS = -std=c++11 -O3 -fopenmp
NVCCFLAGS = -std=c++11 -O3 -Xcompiler -fopenmp
LDFLAGS = -lgmp
CUDA_LDFLAGS = -lcuda -lcudart -lgmp

SRC_DIR = src
BIN_DIR = bin

SEQ_TARGET = $(BIN_DIR)/miller_rabin_seq
PAR_TARGET = $(BIN_DIR)/miller_rabin_par

COMMON_SRC = $(SRC_DIR)/utils.cpp
SEQ_SRC = $(SRC_DIR)/miller_rabin_seq.cpp $(SRC_DIR)/miller_rabin_cpu.cpp
PAR_SRC = $(SRC_DIR)/miller_rabin_par.cu

all: dirs $(SEQ_TARGET) $(PAR_TARGET)

dirs:
	mkdir -p $(BIN_DIR)

$(SEQ_TARGET): $(SEQ_SRC) $(COMMON_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(PAR_TARGET): $(PAR_SRC) $(COMMON_SRC) $(SRC_DIR)/miller_rabin_cpu.cpp
	$(NVCC) $(NVCCFLAGS) -o $@ $^ $(CUDA_LDFLAGS)

clean:
	rm -rf $(BIN_DIR)

.PHONY: all dirs clean 