CXX = g++

CXXFLAGS = -O3 -std=c++17 -Wall -Wextra
LDLIBS = -larmadillo

TARGET = simulator.exe

OBJ = system.o NSL_SIMULATOR.o particle.o random.o


# Default target
all: $(TARGET)


# Link executable
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET) $(LDLIBS)


# Compile source files
NSL_SIMULATOR.o: NSL_SIMULATOR.cpp system.h particle.h
	$(CXX) $(CXXFLAGS) -c NSL_SIMULATOR.cpp -o NSL_SIMULATOR.o

system.o: system.cpp system.h random.h particle.h
	$(CXX) $(CXXFLAGS) -c system.cpp -o system.o

particle.o: particle.cpp particle.h
	$(CXX) $(CXXFLAGS) -c particle.cpp -o particle.o

random.o: random.cpp random.h
	$(CXX) $(CXXFLAGS) -c random.cpp -o random.o


# Compile and run
run: $(TARGET)
	./$(TARGET)


# Remove object files and executable
clean:
	rm -f $(OBJ) $(TARGET)


# Remove output files
remove:
	rm -f ../OUTPUT/*.*
	rm -f ../OUTPUT/CONFIG/*.*


# Remove compiled files and output files
cleanall: clean remove


.PHONY: all run clean remove cleanall
