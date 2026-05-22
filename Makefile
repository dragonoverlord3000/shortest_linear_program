CXX := g++
AR := ar

USE_MATHOPT ?= 0

ORTOOLS_CXXFLAGS :=
ORTOOLS_LDFLAGS :=
ORTOOLS_LDLIBS :=

ifeq ($(USE_MATHOPT),1)
ORTOOLS_DIR ?= $(HOME)/opt/or-tools_x86_64_Ubuntu-24.04_cpp_v9.12.4544

ORTOOLS_CXXFLAGS := -isystem $(ORTOOLS_DIR)/include -DOR_PROTO_DLL= -DSLP_WITH_MATHOPT
ORTOOLS_LDFLAGS := -Wl,-rpath,$(ORTOOLS_DIR)/lib
ORTOOLS_LDLIBS := $(sort $(wildcard $(ORTOOLS_DIR)/lib/*.so))
endif

CXXFLAGS := -std=c++23 -Ofast -march=native -Iinclude -Ithird_party -Wall -Wextra $(ORTOOLS_CXXFLAGS)

BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
LIB_DIR := $(BUILD_DIR)/lib
EXAMPLE_BIN_DIR := $(BUILD_DIR)/examples
APP_BIN_DIR := $(BUILD_DIR)/bin
BENCH_DIR := $(BUILD_DIR)/benchmarks

LIB := $(LIB_DIR)/libslp.a

SRC := src/slp/algorithm.cpp \
	src/slp/potential/gf2/core.cpp \
	src/slp/potential/gf2/greedy.cpp \
	src/slp/potential/gf2/backtrack.cpp \
	src/slp/potential/ternary/greedy.cpp \
	src/slp/boyar_peralta/gf2/core.cpp \
	src/slp/boyar_peralta/gf2/BP.cpp \
	src/slp/boyar_peralta/gf2/RNBP.cpp \
	src/slp/boyar_peralta/gf2/Ax.cpp \
	src/slp/paar/gf2/core.cpp \
	src/slp/paar/gf2/greedy.cpp \
	src/slp/mip/mathopt.cpp \
	src/slp/utils/utils.cpp \
	src/slp/preprocess/preprocess.cpp \
	src/slp/postprocess/postprocess.cpp \
	src/slp/framework/framework.cpp

OBJ := $(patsubst src/%.cpp,$(OBJ_DIR)/%.o,$(SRC))

EXAMPLE_SRC := $(wildcard examples/*.cpp)
EXAMPLES := $(patsubst examples/%.cpp,$(EXAMPLE_BIN_DIR)/%,$(EXAMPLE_SRC))
APP_SRC := $(wildcard apps/*.cpp)
APPS := $(patsubst apps/%.cpp,$(APP_BIN_DIR)/%,$(APP_SRC))
BENCH_BIN := $(BENCH_DIR)/runner
BENCH_SRC := benchmarks/runner.cpp \
             benchmarks/io.cpp \
             benchmarks/3x3_matmul/bench.cpp \
             benchmarks/3x3_matmul/matrix.cpp \
             benchmarks/crypt/bench.cpp \
             benchmarks/struct_matmul/bench.cpp
BENCH_3x3MATMUL_SCHEMES := benchmarks/3x3_matmul/schemes-tab
BENCH_3x3MATMUL_ARCHIVE := $(BUILD_DIR)/schemes-tab.tgz

all: $(LIB) $(APPS) $(EXAMPLES)

$(BUILD_DIR) $(OBJ_DIR) $(LIB_DIR) $(EXAMPLE_BIN_DIR) $(APP_BIN_DIR) $(BENCH_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: src/%.cpp | $(OBJ_DIR)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LIB): $(OBJ) | $(LIB_DIR)
	$(AR) rcs $@ $^

$(EXAMPLE_BIN_DIR)/%: examples/%.cpp $(LIB) | $(EXAMPLE_BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $< -L$(LIB_DIR) -lslp $(ORTOOLS_LDFLAGS) $(ORTOOLS_LDLIBS)

$(APP_BIN_DIR)/%: apps/%.cpp $(LIB) | $(APP_BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $< -L$(LIB_DIR) -lslp $(ORTOOLS_LDFLAGS) $(ORTOOLS_LDLIBS)

apps: $(APPS)

examples: $(EXAMPLES)

$(BENCH_BIN): $(BENCH_SRC) $(LIB) | $(BENCH_DIR)
	$(CXX) $(CXXFLAGS) \
		'-DSLP_COMPILER="$(CXX)"' \
		'-DSLP_CXXFLAGS="$(CXXFLAGS)"' \
		-o $@ $(BENCH_SRC) -L$(LIB_DIR) -lslp $(ORTOOLS_LDFLAGS) $(ORTOOLS_LDLIBS)

$(BENCH_3x3MATMUL_SCHEMES): | $(BUILD_DIR)
	rm -rf $@
	mkdir -p $@
	wget -O $(BENCH_3x3MATMUL_ARCHIVE) "http://www.algebra.uni-linz.ac.at/research/matrix-multiplication/schemes-tab.tgz"
	tar -xzf $(BENCH_3x3MATMUL_ARCHIVE) -C benchmarks/3x3_matmul/schemes-tab
	test -d $@

download-3x3: $(BENCH_3x3MATMUL_SCHEMES)
download-crypt:
	python3 benchmarks/crypt_download_helper.py
download-struct:
	python3 benchmarks/struct_matmul_download_helper.py
download-bernoulli:
	python3 benchmarks/bernoulli_random_download_helper.py

bench-full: $(BENCH_BIN) $(BENCH_3x3MATMUL_SCHEMES)
	./$(BENCH_BIN) $(BENCH_ARGS)

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean apps examples bench-full download-3x3 download-crypt download-struct download-bernoulli
