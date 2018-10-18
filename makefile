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
#  -DF08         : Fortran 2008 language is accepted or not
#  -DMKL_SPARSE  : use MKL_SPARSE library for the sparse matrix routines
#                  This option will save your memory 'enormously'.
#                  Before activate this option, make sure that the file
#                  mkl_spblas.f90 is exist in $(MKLPATH)/include/mkl_spblas.f90
#                  If this option is activated, you can use EWINDOW tag in
#                  your input file. See the manual for the details.
#############################################################################
OPTIONS= -fpp -DMPI -DF08 -DSPGLIB -DMKL_SPARSE
F90    = mpif90 $(OPTIONS)
FFLAG  = -O2 -heap-arrays -nogen-interfaces
BIN    = ~/code/bin
#---------------------------------------------------------------------------|

#-----------------------------------
# Dependencies: LAPACK, SPGLIB     |
#---------------------------------------------------------------------------|
SPGLIB = -L/Users/Infant/code/lib/ -lsymspg
MKLPATH= $(MKLROOT)
LAPACK = -L$(MKLPATH)/lib/ \
         -lmkl_intel_lp64 -lmkl_sequential \
         -lmkl_core -liomp5
BLAS   = 
INCLUDE= -I$(MKLPATH)/include
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

MPI_MOD= mpi_basics.o mpi_setup.o
TEST   = test.o
MODULE = $(MPI_MOD) time.o $(SP_MOD) \
		 parameters.o read_incar.o orbital_wavefunction.o \
		 kronecker_prod.o phase_factor.o do_math.o \
         sorting.o berry_phase.o sparse_tool.o
READER = parse.o read_input.o read_param.o read_poscar.o read_kpoint.o \
		 read_energy.o set_weight.o get_site_number.o find_nn.o 
WRITER = plot_eigen_state.o plot_stm_image.o set_ribbon_geom.o print_energy.o \
		 print_wcc.o print_zak_phase.o print_berry_curvature.o
GET    = get_tij.o get_eig.o get_dos.o get_soc.o get_param_class.o \
		 get_cc_param.o get_berry_curvature.o get_wcc.o get_zak_phase.o \
         get_z2_invariant.o get_parity.o get_hamk_sparse.o
SYMM   = spglib_interface.o get_symmetry.o 
MINPACK_LIB= get_fit.o minpack_sub.o lmdif.o

SPG    =$(findstring -DSPGLIB,$(OPTIONS))

ifeq ($(SPG),-DSPGLIB)
  SPGLIB_=  $(SPGLIB)
  OBJECTS=  $(MODULE) tbfit.o tbpack.o $(READER) $(WRITER) $(GET) \
                      e_onsite.o $(MINPACK_LIB) $(TEST) $(SYMM)
else
  SPGLIB_= 
  OBJECTS=  $(MODULE) tbfit.o tbpack.o $(READER) $(WRITER) $(GET) \
                      e_onsite.o $(MINPACK_LIB) $(TEST)
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
	$(F90) -o $@ $^ $(BLAS) $(LAPACK) $(SPGLIB_) $(INCLUDE)
	cp tbfit $(BIN)

clean:
	rm $(BIN)/tbfit *.o *.mod

