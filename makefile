CXX = g++
CXXLIBS = -lgsl -lgslcblas -lm
CXXOPTS = -march=native -mtune=native -O3
OBJECTS = main.o gauss.o power.o harppi.o

build: $(OBJECTS)
	$(CXX) $(CXXLIBS) $(CXXOPTS) $^ -o gausspk
	
%.o: %.cpp
	$(CXX) $(CXXOPTS) -c $< -o $@
	
clean:
	rm *.o gausspk
