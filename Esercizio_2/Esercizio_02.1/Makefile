# Compiler
CXX = g++

# Compilation options
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Executable
TARGET = main

# Source files
SRC = main.cpp random.cpp

# Header files
HEADERS = $(wildcard *.h *.hpp)


# Default target
all: $(TARGET)


# Compile program
$(TARGET): $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)


# Compile and run
run: $(TARGET)
	./$(TARGET) $(ARGS)


# Debug version
debug:
	$(CXX) -std=c++17 -Wall -Wextra -g -O0 $(SRC) -o $(TARGET)


# Remove compiled executable
clean:
	rm -f $(TARGET)


# Remove executable and object files
cleanall: clean
	rm -f *.o


.PHONY: all run debug clean cleanall