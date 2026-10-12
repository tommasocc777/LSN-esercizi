CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra
LDLIBS = -larmadillo

TARGET = simulator.exe
OBJ = system.o particle.o random.o NSL_SIMULATOR.o

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET) $(LDLIBS)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all run clean
