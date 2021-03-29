# Source directory
SRC_DIR := ../../src/Serial

INC := -I$(SRC_DIR)
LIB := -lm

CC	:= gcc-10
# CC  := clang
# CC	:= icc

CFLAGS := -O3 -Winline -Wall
# CFLAGS := $(CFLAGS) -march=native

# Flag to turn off the correction step via CDF evaluation
# CFLAGS := $(CFLAGS) -D NO_CDF

# Flag to specify whether to run poissinv_timing with just DOUBLE or SINGLE routines
# CFLAGS := $(CFLAGS) -D DOUBLE
# CFLAGS := $(CFLAGS) -D SINGLE

# Flag to turn on the checks on the poisscinv function
CFLAGS := $(CFLAGS) -D COMPLEMENTARY

# Flag to turn on the checks on the vector versions
CFLAGS := $(CFLAGS) -D VECTOR

# Flag to turn off the checks on the scalar versions
# CFLAGS := $(CFLAGS) -D SKIP_SCALAR

# Flag to speed up the poissinv_check code by using shortcuts
CFLAGS := $(CFLAGS) -D FAST_CHECK
ifeq ($(CC), icc)
	# When using FAST_CHECK test need to tell Intel to not flush to zero
	CFLAGS := $(CFLAGS) -no-ftz
endif

# Flags to test accuracy by increments on a linear scale, rather than the standard
# increments on a log scale
# CFLAGS := $(CFLAGS) -D PLUS
# CFLAGS := $(CFLAGS) -D COUNT_LAMBDA=100000

# Dependencies
DEPS := $(SRC_DIR)/poissinv.h \
		$(SRC_DIR)/poissinvf.h

all:	poissinv_timing poissinv_check

poissinv_timing:	poissinv_timing.c $(DEPS) Makefile
	${CC} $(CFLAGS) poissinv_timing.c -o poissinv_timing $(INC) $(LIB)

poissinv_check:		poissinv_check.c poissinv_quad.o $(DEPS) Makefile
	${CC} $(CFLAGS) poissinv_check.c poissinv_quad.o -o poissinv_check $(INC) $(LIB) -lquadmath

poissinv_quad.o:	poissinv_quad.c Makefile
# CLANG does not find libquadmath so compile the poissinv_quad with gcc instead
ifeq ($(CC) , clang)
	gcc $(CFLAGS) -c poissinv_quad.c
else
	${CC} $(CFLAGS) -c poissinv_quad.c
endif

# Timing routine to check performance across different compilers
timing_compilers:   poissinv_timing.c $(DEPS) Makefile
	gcc-10 	$(CFLAGS) poissinv_timing.c -o poissinv_timing_gcc   $(INC) $(LIB)
	clang 	$(CFLAGS) poissinv_timing.c -o poissinv_timing_clang $(INC) $(LIB)
	icc 	$(CFLAGS) poissinv_timing.c -o poissinv_timing_icc   $(INC) $(LIB)

clean:
	rm -f poissinv_timing poissinv_check poissinv_quad.o poissinv_timing_*
