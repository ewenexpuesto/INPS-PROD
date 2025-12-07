CXX := g++

TARGET := main
SOURCES := Main.cpp Basis.cpp Hermit.cpp Poly.cpp
INCLUDES := -Iarmadillo-code-15.0.x/armadillo-code-15.0.x/include

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(SOURCES) $(INCLUDES) -o $(TARGET)

clean:
	rm -f $(TARGET)
