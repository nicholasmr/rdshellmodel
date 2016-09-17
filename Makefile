#----------------------
# GENERAL SETUP
#----------------------

OUTDIR=../../debug/
MULTITHREADED=0

# Integration status
RESUME=0
DO_AGGREGATIONS=0
RESUME_AGGREGATIONS=0

# The model
MODEL=12
NSH=148
LAMBDA_m=1.1
LAMBDA_g=$(LAMBDA_m)
LAMBDA_G=$(LAMBDA_m)
DEMARC_MODEL=0
DEMARC_VALUE=0
FIXED_P=1
FIXED_Q=2

MODELCONFIG=$(MODEL) $(LAMBDA_m) $(LAMBDA_g) $(LAMBDA_G) $(DEMARC_MODEL) $(DEMARC_VALUE) $(FIXED_P) $(FIXED_Q) $(NSH)

# Forcing, dissipation, etc.
MACROS__MAIN__FORCE=-DFTYPE=-2 -DFSH=48 -DDFSH=1 -DFMAGNITUDE=1 
MACROS__MAIN__DISSP=-DDTYPE=3 -DVISC_SMALLSCALE=1.0D-8 -DVISC_LARGESCALE=1.0D+3 -DLSKDEP=-4
MACROS__MAIN__NUMERICS=-DNT=100 -DNTI=1000 -DDT=1.0D-7 -DKZERO=1
MACROS__MAIN__OTHER=-DRESUME=$(RESUME) -DMULTITHREADED=$(MULTITHREADED) -DOUTDIR=$(OUTDIR) -DDO_AGGREGATIONS=$(DO_AGGREGATIONS) -DRESUME_AGGREGATIONS=$(RESUME_AGGREGATIONS)

MACROS__MAIN=$(MACROS__MAIN__FORCE) $(MACROS__MAIN__DISSP) $(MACROS__MAIN__NUMERICS) $(MACROS__MAIN__OTHER)

#----------------------
# COMPILER SETUP: GFORTRAN
#----------------------

#COMPILER=gfortran 
#OPTS=-funroll-loops -fdefault-real-8 -cpp -O2 -ftree-vectorize -ftree-vectorizer-verbose=2 -m64 -march=native -msse2 
##OPTS=-fimplicit-none -funroll-loops -ffree-line-length-none -finit-local-zero -finline-functions -fdefault-real-8 -cpp -O2
##OPTM=-fbacktrace -Wsurprising -Wunderflow -Walign-commons -fcheck=all #-fcheck=bounds

#----------------------
# COMPILER SETUP: IFORT  (ALWAYS USE ifort IF POSSIBLE. IT IS THE ONLY COMPILER WHICH CAN VECTORISE DOUBLE COMPLEX ARRAYS CORRECTLY.). NO REAL EFFECT FOR -ftz of -O3, SO DON'T USE.
#----------------------

# CIC ICEMONSTER: /opt/intel/composer_xe_2011_sp1.9.293/bin/intel64/ifort ---- MAKE SURE THE LD PATH CONTAINS: export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/intel/compiler/lib/intel64"
COMPILER=/opt/intel/composer_xe_2011_sp1.9.293/bin/intel64/ifort
OPTS=-real-size 64 -fpp -O2  -free -xHost -funroll-loops -align all    -vec-report3 -diag-disable 8291 -diag-disable 8290   -mcmodel=medium -i-dynamic   -debug all
ifeq ($(MULTITHREADED),1)
	OPTM=-openmp
endif

#----------------------
# LIBS
#----------------------
# INSTALLING NETCDF LIB:
#----------------------
# - Download both netcdf-fortran-4.4.2 and netcdf-4.3.3. 
# - Compile non-fortran lib by using ifort: see https://software.intel.com/en-us/articles/performance-tools-for-software-developers-building-netcdf-with-the-intel-compilers
#   ---> MAKE SURE THAT ALL COMPILER FLAG FIELDS HAVE: -mcmodel=medium -i-dynamic
#   ---> run: export NCDIR=/cicdisk/home/rathmann/libs/netcdflibnew
#   ---> run: source /opt/intel/composer_xe_2011_sp1.9.293/bin/intel64/../compilervars.sh intel64
#   ---> run: ./configure --disable-netcdf-4 --prefix=${NCDIR}
#   ---> run: make && make install
# - Compile fortran lib by using ifort (same link as above, but in addition see https://www.unidata.ucar.edu/software/netcdf/docs/building_netcdf_fortran.html)
#   ---> run: export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}
#   ---> run: CPPFLAGS=-I${NCDIR}/include LDFLAGS=-L${NCDIR}/lib ./configure --prefix=${NCDIR}
#   ---> run: make && make install
#----------------------

#LIBDIRS=-I../../libs/netcdflib/include -L../../libs/netcdflib/lib # OLD version, does not support mcmodel=medium needed for large dumps.
LIBDIRS=-I../../libs/netcdflibnew/include -L../../libs/netcdflibnew/lib
LIBS=-lnetcdff -lnetcdf

#----------------------
# MAKE PROGRAM
#----------------------
# export NCDIR=/cicdisk/home/rathmann/libs/netcdflibnew
# export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}  # MUST BE SET FOR SHARED LIBS TO BE RECOGNIZED ON COMPILE-TIME. 
# NOTE THAT FOR SLURM SETUPS THE FRONTEND ENV. VARS. ARE PASSED ON DIRECTLY TO NODES.
#----------------------

F90MAIN=main
F90MODELS=shellmodels
F90MODELSPLAIN=$(F90MODELS)_template
F90CPP=$(F90MAIN)_cpp

shellmodel: $(F90MAIN).f90
	mkdir -p $(OUTDIR)
	perl setup.pl $(MODELCONFIG) $(F90MODELSPLAIN).f90 $(F90MODELS).f90
	gfortran -E -cpp $(MACROS__MAIN) $(F90MAIN).f90 | sed 's/LC/\&\n/g' | sed "s/dump\\.nc/dump$(EXT).nc/g" | sed "s/resume\\.nc/resume$(EXT).nc/g" | sed '/^#/d' | sed '/^!/d' > $(F90CPP).f90
	$(COMPILER) $(OPTS) $(OPTM) $(LIBDIRS) -o $(OUTDIR)/shellmodel$(EXT) $(F90CPP).f90 $(LIBS)
	cp $(F90CPP).f90 $(OUTDIR)/shellmodel$(EXT).f90
	cp $(F90CPP).f90 $(F90MAIN)_latest.f90
	rm $(F90CPP).f90 $(F90MODELS).f90

