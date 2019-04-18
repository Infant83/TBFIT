# $^ : dependency list
# $@ : target

#----- Change options and library path according to your system ------------#
#-----------------------------------
# Compiler options and bin path    |
#---------------------------------------------------------------------------|
#OPTIONS= -fpp -DMPI -mcmodel=large # for curion2
################# Possible options ##########################################
#  -DSPGLIB      : printout spacegroup information in the initial stages
#                : if you want to use this option, please link SPGLIB 
#                  library path properly in the "Dependencies" section below
#  -DMPI         : MPI paralallism activation 
#                  If -DSCALAPACK option is activated: k-point + eigenvalue 
#                  parallism will be imployed,
#                  otherwise, only the k-point parallization will be performed.   
#  -DF08         : Fortran 2008 language is accepted or not
#  -DMKL_SPARSE  : use MKL_SPARSE library for the sparse matrix routines
#                  This option will save your memory 'enormously'.
#                  Before activate this option, make sure that the file
#                  mkl_spblas.f90 is exist in $(MKLPATH)/include/mkl_spblas.f90
#                  If this option is activated, you can use EWINDOW tag in
#                  your input file. See the manual for the details.
#  -DSCALAPACK   : use ScaLAPACK library for the eigenvalue parallism 
#                  !!! WARN !!! do not use in the current version: it is upon
#                               developing stage now.
#############################################################################
#OPTIONS= -fpp -DF08 -DSPGLIB -DMKL_SPARSE #-DSCALAPACK
OPTIONS= -fpp -DMPI -DF08 -DSPGLIB -DMKL_SPARSE #-DSCALAPACK
F90    = mpif90 $(OPTIONS)
FFLAG  = -O2 -heap-arrays -nogen-interfaces

#OPTIONS= -cpp -DMPI -DF08 -DSPGLIB #-DMKL_SPARSE -DSCALAPACK
#F90    = mpif90-openmpi-mp $(OPTIONS)
#FFLAG  = -O2 -ffree-line-length-512 -fmax-stack-var-size=32768
BIN    = ~/code/bin
#---------------------------------------------------------------------------|

#-----------------------------------
# Dependencies: LAPACK, SPGLIB     |
#---------------------------------------------------------------------------|
#SPGLIB    = -L/Users/Infant/code/lib/ -lsymspg
SPGLIB    = -L/Users/Infant/tbfit_fortran/LIB/spglib-master -lsymspg
MKLPATH   = $(MKLROOT)
LAPACK    = -L$(MKLPATH)/lib/ \
            -lmkl_intel_lp64 -lmkl_sequential \
            -lmkl_core -liomp5
BLAS      = 
INCLUDE   = -I$(MKLPATH)/include
#SCALAPACK = /Users/Infant/tbfit_fortran/LIB/scalapack-2.0.2/libscalapack.a
SCALAPACK = /Users/Infant/tbfit_fortran/LIB/scala_home/libscalapack.a
#---------------------------------------------------------------------------|


######################### Do not modify below ###############################
#-----------------------------------
# Objects                          |
#---------------------------------------------------------------------------|
MKL_SP =$(findstring -DMKL_SPARSE,$(OPTIONS))
MKL_SPARSE = mkl_spblas.o
ifeq ($(MKL_SP),-DMKL_SPARSE)
  SP_MOD = mkl_spblas.o
else 
  SP_MOD = 
endif
SCALAPACK_USE=$(findstring -DSCALAPACK,$(OPTIONS))

MPI_MOD= blacs_basics.o mpi_basics.o mpi_setup.o 
TEST   = test.o
MODULE = $(MPI_MOD) memory.o time.o version.o $(SP_MOD) \
		 parameters.o read_incar.o orbital_wavefunction.o \
		 kronecker_prod.o phase_factor.o do_math.o \
         sorting.o berry_phase.o sparse_tool.o pikaia_module.o geodesiclm.o
READER = parse.o read_input.o read_param.o read_poscar.o read_kpoint.o \
		 read_energy.o set_weight.o get_site_number.o find_nn.o
WRITER = plot_eigen_state.o plot_stm_image.o set_ribbon_geom.o print_energy.o \
		 print_wcc.o print_zak_phase.o print_berry_curvature.o replot_dos.o
GET    = get_tij.o get_eig.o get_dos.o get_soc.o get_param_class.o \
		 get_cc_param.o get_berry_curvature.o get_wcc.o get_zak_phase.o \
         get_z2_invariant.o get_parity.o get_hamk_sparse.o get_effective_ham.o \
         e_onsite.o
SYMM   = spglib_interface.o get_symmetry.o 
FITTING_LIB= get_fit.o minpack_sub.o lmdif.o genetic_alorithm.o

ifeq ($(SCALAPACK_USE), -DSCALAPACK)
  SCALAPACK_LIB= $(SCALAPACK)
  SCALAPACK_OBJ= #scalapack_initialize.o
else
  SCALAPACK_LIB=
  SCALAPACK_OBJ= 
endif

SPG    =$(findstring -DSPGLIB,$(OPTIONS))

ifeq ($(SPG),-DSPGLIB)
  SPGLIB_=  $(SPGLIB)
  OBJECTS=  $(MODULE) tbfit.o tbpack.o $(READER) $(WRITER) $(GET) \
                      $(FITTING_LIB) $(SCALAPACK_OBJ) $(TEST) $(SYMM)
else
  SPGLIB_= 
  OBJECTS=  $(MODULE) tbfit.o tbpack.o $(READER) $(WRITER) $(GET) \
                      $(FITTING_LIB) $(SCALAPACK_OBJ) $(TEST)
endif


#---------------------------------------------------------------------------|

#-----------------------------------
# Suffix rules                     |
#-----------------------------------
ifeq ($(MKL_SP),-DMKL_SPARSE)
.SUFFIXES: $(MKL_SPARSE)
$(MKL_SPARSE): $(MKLPATH)/include/mkl_spblas.f90
	$(F90) $(FFLAG) -c $<
endif
.SUFFIXES: .f .f90 
%.o: %.f90
	$(F90) $(FFLAG) -c $<


#-----------------------------------
# Targets                          |
#-----------------------------------
#$(BIN)/tbfit: $(OBJECTS) 
tbfit: $(OBJECTS) 
	$(F90) -o $@ $^ $(BLAS) $(LAPACK) $(SCALAPACK_LIB) $(SPGLIB_) $(INCLUDE)
	cp tbfit $(BIN)

poscar2bs: poscar2bs.o
	$(F90) -o $@ $^ 
	cp pc2xyz $(BIN)

all: $(OBJECTS)
	$(F90) -o tbfit $^ $(BLAS) $(LAPACK) $(SCALAPACK_LIB) $(SPGLIB_) $(INCLUDE)
	$(F90) -o poscar2bs poscar2bs.f90
	cp tbfit $(BIN)
	cp pc2xyz $(BIN)

clean:
	rm $(BIN)/tbfit *.o *.mod

clean_poscar2xyz:
	rm $(BIN)/pc2xyz poscar2bs.o
