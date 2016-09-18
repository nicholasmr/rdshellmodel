#If you have questions or are experiencing any problems, please contact <rathmann@nbi.ku.dk>. 

#----------------------
# GENERAL SETUP
#----------------------

OUTDIR=./
MULTITHREADED=0         # Use OpenMP multithreading?

# Integration status
#-------------------
RESUME=0                # Resume from previously dumped .nc file? (file must be re-named from "dump.nc" to "resume.nc").
DO_AGGREGATIONS=0       # Aggregate over all time steps correlator values (corr_*), u^+_n and u^-_n (u*_abs), Pi_^E_n (Eflux) and structure functions (structfuncs) (see main.f90 for details, also note some of these may manually be deactivated in shellmodels_template.f90).
RESUME_AGGREGATIONS=0   # If resuming from resume.nc, should the already saved aggregated statistics be neglected or not?   

# The model
#----------
MODEL=15                # Models 11-14 refer to the individual submodels 1-4, whereas model 15 refers to the coupled submodels 1-4 (Rathmann and Ditlevsen, 2016).  
NSH=148                 # Number of shells
LAMBDA_m=1.1            # Lambda value used for epsilon_{p,q}^{s',s''} and xi_{p,q}^{s',s''}.
LAMBDA_g=$(LAMBDA_m)    # Lambda value used for g_{p,q}^{s',s''}
LAMBDA_G=$(LAMBDA_m)    # Lambda value used for G_{p,q}^{s',s''}
FIXED_P=1               # Comma separated list (no spaces) of p-values to include (if more than one is specified a multi triad shape system is simulated). A value of zero (0) implies simulating all values allowed by the above specified Lambda. 
FIXED_Q=2               # ... similarly as for FIXED_P, but for q.
DEMARC_MODEL=0          # Demarcate submodel weights (g) towards of away from this submodel (see setup.pl for details). Zero (0) implies no demarcation.
DEMARC_VALUE=0          # Demarcation value.

# Forcing
#--------
FTYPE=-2        # Forcing type - see shellmodels_template.f90 for details.
FSH=48          # Forcing shell
DFSH=1          # Number of shells to force above FSH (1=>only shell FSH)
FMAGNITUDE=1    # Forcing magnitude

# Dissipation
#------------
DTYPE=3                 # Dissipation type - see shellmodels_template.f90 for details.
VISC_SMALLSCALE=1.0D-8  # Small-scale viscocities
VISC_LARGESCALE=1.0D+3  # Large-scale viscocities
LSKDEP=-4               # Large-scale k-dependence for the dissipation (the small-scale one is fixed at -2 as the Navier-Stokes equation prescribes, i.e. nu*k^2*u_n^{+,-})

# Numerics
#---------
NT=100          # Number of outer-loop integration steps for which velocity profile is saved (see main.f90 for details).
NTI=1000        # Number of inner-loop integration steps 
DT=1.0D-7       # dt
KZERO=1         # k_0 (smallest wave-number resolved)

#---------------------------
MODELCONFIG=$(MODEL) $(LAMBDA_m) $(LAMBDA_g) $(LAMBDA_G) $(DEMARC_MODEL) $(DEMARC_VALUE) $(FIXED_P) $(FIXED_Q) $(NSH)
MACROS__MAIN__FORCE=-DFTYPE=$(FTYPE) -DFSH=$(FSH) -DDFSH=$(DFSH) -DFMAGNITUDE=$(FMAGNITUDE)  
MACROS__MAIN__DISSP=-DDTYPE=$(DTYPE) -DVISC_SMALLSCALE=$(VISC_SMALLSCALE) -DVISC_LARGESCALE=$(VISC_LARGESCALE) -DLSKDEP=$(LSKDEP)
MACROS__MAIN__NUMERICS=-DNT=$(NT) -DNTI=$(NTI) -DDT=$(DT) -DKZERO=$(KZERO)
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

