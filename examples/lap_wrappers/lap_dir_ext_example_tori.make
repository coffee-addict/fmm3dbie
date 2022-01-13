EXEC = int2-ext
HOST = gcc-openmp

LIBS = -lfmm3d -lfmm3dbie 
ifeq ($(HOST),gcc)
    FC=gfortran -L${LDF} 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -std=legacy 
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy 
endif

-include make.inc

SURF=../../src/surface_routs

.PHONY: all clean 

OBJECTS = lap_dir_ext_example_tori.o

%.o : %.f90  
	$(FC) -c $(FFLAGS) $< -o $@

all: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) -L$(LD_LIBRARY_PATH) $(LIBS) 
	./$(EXEC)  

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)

