# FCST / Fast FCST
#   make            C library, tests and C benchmark (no dependencies)
#   make test       run the tests against a naive suffix tree
#   make bench_cst  comparison benchmark (needs SDSL-lite and CST-CN; see scripts/install_deps.sh)

CC       ?= gcc
CXX      ?= g++
CFLAGS   ?= -O3 -march=native -Wall -Wextra -std=c11
CXXFLAGS ?= -O3 -march=native -DNDEBUG -std=c++14 -Wno-deprecated-declarations -Wno-template-body

DEPS  ?= $(CURDIR)/deps
SDSL  ?= $(DEPS)/sdsl
CSTCN ?= $(DEPS)/CST-CN

BUILD = build
COBJS = $(BUILD)/bits.o $(BUILD)/wm.o $(BUILD)/fcst.o

all: $(BUILD)/libfcst.a $(BUILD)/test_fcst $(BUILD)/bench_c

$(BUILD):
	mkdir -p $@

$(BUILD)/%.o: src/c/%.c src/c/*.h | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/libfcst.a: $(COBJS)
	ar rcs $@ $^

$(BUILD)/test_fcst: tests/test_fcst.c $(BUILD)/libfcst.a
	$(CC) $(CFLAGS) -Isrc/c -o $@ $< $(BUILD)/libfcst.a

$(BUILD)/bench_c: bench/bench_c.c $(BUILD)/libfcst.a
	$(CC) $(CFLAGS) -D_POSIX_C_SOURCE=199309L -Isrc/c -o $@ $< $(BUILD)/libfcst.a

test: $(BUILD)/test_fcst
	./$(BUILD)/test_fcst

SDSL_FLAGS = -I$(SDSL)/include -L$(SDSL)/lib -lsdsl -ldivsufsort -ldivsufsort64

bench_cst: $(BUILD)/bench_cst $(BUILD)/csasize

$(BUILD)/bench_cst: bench/bench_cst.cpp src/sdsl/fcst_sdsl.hpp $(BUILD)/libfcst.a
	$(CXX) $(CXXFLAGS) -Isrc/c -Isrc/sdsl -I$(CSTCN)/include -o $@ $< $(BUILD)/libfcst.a $(SDSL_FLAGS)

$(BUILD)/csasize: bench/csasize.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(SDSL_FLAGS)

clean:
	rm -rf $(BUILD)

.PHONY: all test bench_cst clean
