# $^ : dependency list
# $@ : target

#----- Change options and library path according to your system ------------#
#-----------------------------------
# Compiler options and bin path    |
#---------------------------------------------------------------------------|
#OPTIONS= -fpp -DMPI -mcmodel=large # for curion2
################# Possible options ##########################################
#        -DSPGLIB : printout spacegroup information in the initial stages
#                 : if you want to use this option, please link SPGLIB 
#                   library path properly in the "Dependencies" section below
#        -DMPI    : MPI paralallism activation
#############################################################################

OPTIONS= -fpp -DMPI -DSPGLIB
F90    = mpif90 $(OPTIONS)
FFLAG  = -O3 -heap-arrays -nogen-interfaces
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
#---------------------------------------------------------------------------|


######################### Do not modify below ###############################
#-----------------------------------
# Objects                          |
#---------------------------------------------------------------------------|
MPI_MOD= mpi_basics.o mpi_setup.o
TEST   = test.o
MODULE = $(MPI_MOD) \
		 parameters.o read_incar.o orbital_wavefunction.o \
		 kronecker_prod.o phase_factor.o do_math.o \
         sorting.o berry_phase.o spglib_interface.o 
READER = read_input.o read_param.o read_poscar.o read_kpoint.o \
		 read_energy.o set_weight.o get_site_number.o find_nn.o 
WRITER = plot_eigen_state.o plot_stm_image.o set_ribbon_geom.o print_energy.o \
		 print_wcc.o print_zak_phase.o print_berry_curvature.o
GET    = get_tij.o get_eig.o get_dos.o get_soc.o get_param_class.o \
		 get_cc_param.o get_berry_curvature.o get_wcc.o get_zak_phase.o \
         get_z2_invariant.o get_parity.o get_symmetry.o
MINPACK_LIB= get_fit.o minpack_sub.o lmdif.o

OBJECTS=$(MODULE) tbfit.o tbpack.o $(READER) $(WRITER) $(GET) \
        e_onsite.o $(MINPACK_LIB) $(TEST)
#---------------------------------------------------------------------------|

#CC     = gcc
#AR     = ar
#CFLAG  = -g -O3
#-----------------------------------
# Suffix rules                     |
#-----------------------------------
.SUFFIXES: .f .f90 #.mod
%.o: %.f
	$(F90) $(FFLAG) -c $<
%.o: %.f90
	$(F90) $(FFLAG) -c $<

#-----------------------------------
# Targets                          |
#-----------------------------------
$(BIN)/tbfit: $(OBJECTS) 
	$(F90) -o $@ $^ $(LAPACK) $(SPGLIB)

#mv $@ $(BIN)

clean:
	rm $(BIN)/tbfit *.o *.mod

