CC := g++
CXXFLAGS := -std=c++20 -O2 -DARMA_DONT_USE_WRAPPER -DARMA_DONT_USE_BLAS -DARMA_DONT_USE_LAPACK

TARGET := main
SOURCES := Main.cpp Basis.cpp Hermit.cpp Poly.cpp Rho.cpp
INCLUDES := -Iarmadillo-code-15.0.x/armadillo-code-15.0.x/include

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CC) $(CXXFLAGS) $(SOURCES) $(INCLUDES) -o $(TARGET)

clean:
	rm -f $(TARGET)
