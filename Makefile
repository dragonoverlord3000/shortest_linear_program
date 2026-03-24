CXX := g++
AR := ar
CXXFLAGS := -std=c++23 -O3 -march=native -Iinclude -Ithird_party -Wall -Wextra

BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
LIB_DIR := $(BUILD_DIR)/lib
BIN_DIR := $(BUILD_DIR)/examples

LIB := $(LIB_DIR)/libslp.a

SRC := \
	src/slp/algorithm.cpp \
	src/slp/potential/gf2/core.cpp \
	src/slp/potential/gf2/greedy.cpp \
	src/slp/potential/gf2/backtrack.cpp

OBJ := $(patsubst src/%.cpp,$(OBJ_DIR)/%.o,$(SRC))

EXAMPLES := $(BIN_DIR)/use_gf2

all: $(LIB) $(EXAMPLES)

$(OBJ_DIR) $(LIB_DIR) $(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: src/%.cpp | $(OBJ_DIR)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LIB): $(OBJ) | $(LIB_DIR)
	$(AR) rcs $@ $^

$(BIN_DIR)/use_gf2: examples/use_gf2.cpp $(LIB) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $< -L$(LIB_DIR) -lslp

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean
