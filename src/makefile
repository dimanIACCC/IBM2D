#
# Makefile IBM2D programm for Linux
#
# Two variants of making:
# 1. for instaled Boost library with devs packages (version > 1.60.0)
# 2. for manually compiled Boost library (version > 1.60.0)

# Compiler
cc = g++

# Object files
obj = stdafx.o String.o GeomVec.o Exact_solutions.o Boundary_initial_conditions.o Matrix.o Parameters.o helmholtz.o PredictVel.o BiCGStab.o SolidBody.o \
      Calculate_press.o CalculateForce.o Calculate_u_p.o Output.o BodyOfProgram.o \
      ImmersedBoundaryMeth.o

MKLROOT = /opt/intel/mkl/10.2.5.035

BOOSTROOT = /mnt/storage/home/dskuranakov/IBM2D/boost_1_64_0_linux

ibm2d: $(obj)
	$(cc) -std=c++11 -fopenmp -O2 -o ibm2d  $(obj)  -Wl,--no-as-needed -L${MKLROOT}/lib/em64t  -lmkl_solver_ilp64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm  -L${BOOSTROOT}/stage/lib -lboost_system -lboost_filesystem
#	mv ibm2d ../bin

%.o : %.cpp
	$(cc) -std=c++11 -fopenmp -c -DMKL_ILP64 -m64  -I${BOOSTROOT}/  -I${MKLROOT}/include/  $<

clean:
	rm -f *.o