# Esercizio 10.1 - MPI Genetic Algorithm
CXX = mpicxx
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra
TARGET = main.exe
SRC = main.cpp random.cpp
NP ?= 4
MODE ?= migration
INTERVAL ?= 50
GEOMETRY ?= both

.PHONY: all run clean
all: $(TARGET)

$(TARGET): $(SRC) random.h
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

run: all
	mpiexec -np $(NP) ./$(TARGET) --mode $(MODE) --interval $(INTERVAL) --geometry $(GEOMETRY)

clean:
	rm -f $(TARGET)
