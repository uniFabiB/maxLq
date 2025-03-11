compiler = mpif90
progName = prog

#OBJ = global_variables.o initialize.o data_ops.o fftwfunction.o function_ops.o optimization.o maxdLpdtKappaTest.o
OBJ = global_variables.o initialize.o data_ops.o fftwfunction.o function_ops.o optimization.o maxdLpdt.o

BASIC_LIB  = -lm -lmpi

DEB_OPTS = -g -O0
#DEB_OPTS = -g -O0 -fcheck=all -Wall -fbacktrace
# -g enable debugging
# -O0 optimization off (for faster debugging)
EDI_OPTS = -ffree-line-length-512 
# -ffree-line-length-512 otherwise line split error
# -fallow-argument-mismatch turns argument missmatch in mpi to warnings (not needed anymore since use mpi instead of include mpif.h)
OPTIONS = $(DEB_OPTS) $(EDI_OPTS)

FFTW_DIR      = -I$(EBROOTFFTWMPI)/include
FFTW_LIB      = -lfftw3_mpi -lfftw3 -lm		# Load FFTW3 Library
NETCDF_DIR    = -I$(EBROOTNETCDFMINFORTRAN)/include
NETCDF_LIB    = -lnetcdf -lnetcdff

$(progName): $(OBJ)
	$(compiler) $(OPTIONS) $(OBJ) $(BASIC_LIB) $(NETCDF_LIB) $(FFTW_LIB) -o $@

%.o: %.f90
	$(compiler) -c $(OPTIONS) $(NETCDF_DIR) $(FFTW_DIR) $<

clean: 
	rm -f *.o *.mod *.nc *.dat $(progName)
	rm -f -r output/*.dat output/*.log output/*.nc output/constraintDirs
