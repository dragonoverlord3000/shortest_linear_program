CXX := g++
AR := ar
CXXFLAGS := -std=c++23 -O3 -Iinclude -Ithird_party -Wall -Wextra

BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
LIB_DIR := $(BUILD_DIR)/lib
BIN_DIR := $(BUILD_DIR)/examples
BENCH_DIR := $(BUILD_DIR)/benchmarks

LIB := $(LIB_DIR)/libslp.a

SRC := src/slp/algorithm.cpp \
	src/slp/potential/gf2/core.cpp \
	src/slp/potential/gf2/greedy.cpp \
	src/slp/potential/gf2/backtrack.cpp \
	src/slp/boyar_peralta/gf2/core.cpp \
	src/slp/boyar_peralta/gf2/BP.cpp \
	src/slp/boyar_peralta/gf2/RNBP.cpp \
	src/slp/boyar_peralta/gf2/Ax.cpp \
	src/slp/paar/gf2/core.cpp \
	src/slp/paar/gf2/greedy.cpp \
	src/slp/utils/utils.cpp \
	src/slp/preprocess/preprocess.cpp \
	src/slp/postprocess/postprocess.cpp \
	src/slp/framework/framework.cpp

OBJ := $(patsubst src/%.cpp,$(OBJ_DIR)/%.o,$(SRC))

EXAMPLE_SRC := $(wildcard examples/*.cpp)
EXAMPLES := $(patsubst examples/%.cpp,$(BIN_DIR)/%,$(EXAMPLE_SRC))
BENCH_BIN := $(BENCH_DIR)/runner
BENCH_SRC := benchmarks/runner.cpp \
             benchmarks/3x3_matmul/bench.cpp \
             benchmarks/3x3_matmul/io.cpp \
             benchmarks/3x3_matmul/matrix.cpp
BENCH_3x3MATMUL_SCHEMES := benchmarks/3x3_matmul/schemes-tab
BENCH_3x3MATMUL_ARCHIVE := $(BUILD_DIR)/schemes-tab.tgz

all: $(LIB) $(EXAMPLES)

$(BUILD_DIR) $(OBJ_DIR) $(LIB_DIR) $(BIN_DIR) $(BENCH_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: src/%.cpp | $(OBJ_DIR)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LIB): $(OBJ) | $(LIB_DIR)
	$(AR) rcs $@ $^

$(BIN_DIR)/%: examples/%.cpp $(LIB) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $< -L$(LIB_DIR) -lslp

$(BENCH_BIN): $(BENCH_SRC) $(LIB) | $(BENCH_DIR)
	$(CXX) $(CXXFLAGS) \
		'-DSLP_COMPILER="$(CXX)"' \
		'-DSLP_CXXFLAGS="$(CXXFLAGS)"' \
		-o $@ $(BENCH_SRC) -L$(LIB_DIR) -lslp

$(BENCH_3x3MATMUL_SCHEMES): | $(BUILD_DIR)
	rm -rf $@
	mkdir -p $@
	wget -O $(BENCH_3x3MATMUL_ARCHIVE) "http://www.algebra.uni-linz.ac.at/research/matrix-multiplication/schemes-tab.tgz"
	tar -xzf $(BENCH_3x3MATMUL_ARCHIVE) -C benchmarks/3x3_matmul/schemes-tab
	test -d $@

download-3x3: $(BENCH_3x3MATMUL_SCHEMES)

bench-full: $(BENCH_BIN) $(BENCH_3x3MATMUL_SCHEMES)
	./$(BENCH_BIN) $(BENCH_ARGS)

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean bench-full download-3x3
