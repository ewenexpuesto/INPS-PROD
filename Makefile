CC := g++
CXXFLAGS := -Wall -O3 -std=c++14 -Iarmadillo-code-15.0.x/armadillo-code-15.0.x/include \
	-DARMA_DONT_USE_WRAPPER -DARMA_DONT_USE_BLAS -DARMA_DONT_USE_LAPACK
LDFLAGS :=

TARGET := main
SOURCES := src/Main.cpp src/Basis.cpp src/Hermit.cpp src/Poly.cpp src/Solver.cpp
OBJECTS := $(SOURCES:.cpp=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJECTS)
