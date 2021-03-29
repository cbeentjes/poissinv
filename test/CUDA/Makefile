# Source directory
SRC_DIR := ../../src/CUDA

INC	:= -I$(CUDA_HOME)/include -I$(SRC_DIR)
LIB	:= -L$(CUDA_HOME)/lib64 -lcudart -lm

NVCCFLAGS	:= -O3 -arch=sm_35 --ptxas-options=-v --use_fast_math

# Flag to turn off the correction step via CDF evaluation
# NVCCFLAGS := $(NVCCFLAGS) -D NO_CDF

# Flag to specify whether to run poissinv_timing with just DOUBLE or SINGLE routines
# NVCCFLAGS := $(NVCCFLAGS) -D DOUBLE
# NVCCFLAGS := $(NVCCFLAGS) -D SINGLE

# Flags to test accuracy by increments on a linear scale, rather than the standard
# increments on a log scale
# NVCCFLAGS := $(NVCCFLAGS) -D PLUS
# NVCCFLAGS := $(NVCCFLAGS) -D COUNT_LAMBDA=500000

# Dependencies
DEPS := $(SRC_DIR)/poissinv_cuda.h

all:	poissinv_timing poissinv_check

poissinv_timing:	poissinv_timing.cu $(DEPS) Makefile
		nvcc poissinv_timing.cu -o poissinv_timing $(INC) $(NVCCFLAGS) $(LIB)

poissinv_check:		poissinv_check.cu poissinv_quad.o $(DEPS) Makefile
		nvcc poissinv_check.cu poissinv_quad.o -o poissinv_check $(INC) $(NVCCFLAGS) $(LIB) -lquadmath

poissinv_quad.o:	../Serial/poissinv_quad.c
		gcc -c ../Serial/poissinv_quad.c

clean:
		rm -f poissinv_timing poissinv_check poissinv_quad.o
