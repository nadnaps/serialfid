#     ===========================================
#     Makefile for serialfid
#     ===========================================


#     -----------------------------------
#     Compiler 

#FC=h5pfc -O3 -cpp -I/usr/local/hdf5/include -I/usr/include #optimised option
FC=h5pfc -O0 -g -fbacktrace -fbounds-check -cpp -fdefault-real-8 -fdefault-double-8 -I/usr/local/hdf5/include -I/usr/include #debug option



#     -----------------------------------
#     Library 

LINKS = -L/usr/local/fftw3/lib -lfftw3 -llapack -lz -lhdf5_fortran -lhdf5


#     -----------------------------------
#     Make program

PROGRAM = serfid

OBJECTS = param.o aux_routines.o main.o ReadInputFile.o Logs.o InitRoutines.o MemRoutines.o \
		  HdfAuxRoutines.o CreateGrid.o CompCoeff.o InitPressSolv.o CheckRoutines.o QuitRoutines.o \
		  TimeMarcher.o ExplicitTermsVY.o ExplicitTermsVZ.o ImplicitAndUpdateVY.o ImplicitAndUpdateVZ.o \
		  PerTridiagSolve.o SolveTridY.o SolveTridXYZ.o SolveTridZ.o \
		  CalcLocalDivergence.o SolvePressureCorrection.o CorrectVelocity.o CorrectPressure.o \

MODULES = param.o 


#     -----------------------------------
#     Linking

$(PROGRAM): $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) $(LINKS) -o $@

#     ----------------------------------
#     Dependencies

param.o: param.F90
	$(FC) -c param.F90

%.o: %.F90 $(MODULES)
	$(FC) -c $<

#     ---------------------------------
#     Clean up

clean :
	rm *.o 
	rm *.mod
	rm serfid
