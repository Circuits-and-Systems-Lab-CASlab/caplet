CXX     = g++
MPICXX  = mpic++
LIB = -lmpi -lpthread /lib64/libopenblas.so.0
FLAG    = -g -Wall -I$(INCLUDE)  -fpermissive
DEF_OPENMP = -fopenmp -DCAPLET_OPENMP
DEF_MPI = -DCAPLET_MPI

SRC        = src
OBJ_MPI    = obj/mpi
OBJ_OPENMP = obj/openmp
BIN        = bin
INCLUDE    = include
EXAMPLE    = example

CAPLET_MPI_OBJ = \
	$(OBJ_MPI)/caplet.o \
	$(OBJ_MPI)/caplet_elem.o \
	$(OBJ_MPI)/caplet_int.o \
	$(OBJ_MPI)/caplet_widgets.o \
	$(OBJ_MPI)/main.o \

CAPLET_OPENMP_OBJ = \
	$(OBJ_OPENMP)/caplet.o \
	$(OBJ_OPENMP)/caplet_elem.o \
	$(OBJ_OPENMP)/caplet_int.o \
	$(OBJ_OPENMP)/caplet_widgets.o \
	$(OBJ_OPENMP)/main.o \

# Default values (optional, overridden when passed via CLI)
ENABLE_DEBUG_PRINTS ?= 0
ENABLE_P_MATRIX_PRINT ?= 0

# Condition 1: Check general debug prints
ifeq ($(strip $(ENABLE_DEBUG_PRINTS)), 1)
	FLAG += -DENABLE_DEBUG_PRINTS
endif

# Condition 2: Check P-matrix prints
ifeq ($(strip $(ENABLE_P_MATRIX_PRINT)), 1)
	FLAG += -DENABLE_P_MATRIX_PRINT
endif

all: capletMPI capletOpenMP 

mpi: capletMPI 

openmp: capletOpenMP 

capletMPI: $(CAPLET_MPI_OBJ)
	$(MPICXX) $(FLAG) $(DEF_MPI) -o $(BIN)/$@ $^ $(LIB)

$(OBJ_MPI)/%.o: $(SRC)/%.cpp
	$(MPICXX) $(FLAG) $(DEF_MPI) -c $< -o $@

capletOpenMP: $(CAPLET_OPENMP_OBJ)
	$(MPICXX) $(FLAG) $(DEF_OPENMP) -o $(BIN)/$@ $^ $(LIB) 

$(OBJ_OPENMP)/%.o: $(SRC)/%.cpp
	$(MPICXX) $(FLAG) $(DEF_OPENMP) -c $< -o $@

.phony: clean
clean:
	rm -rf $(OBJ_MPI)/* $(OBJ_OPENMP)/* 
