#
# Makefile diflow programm for Linux
#
# Two variants of making:
# 1. for instaled Boost library with devs packages (version > 1.60.0)
# 2. for manually compiled Boost library (version > 1.60.0)

# Compiler
cc = g++

# Object files
obj = stdafx.o String.o GeomVec.o Matrix.o BiCGStab.o Parameters.o SolidBody.o PredictVel.o \
      Calculate_press.o CalculateForce.o Calculate_u_p.o Output.o BodyOfProgram.o \
      Testing.o ImmersedBoundaryMeth.o

ibm: $(obj)
#	$(cc) -O2 -o ibm $(obj)
	$(cc) -O2 -o ibm $(obj) -L/usr/local/lib/ -lboost_system -lboost_filesystem
	mv ibm ./../bin

%.o : %.cpp
	$(cc) -std=c++11 -c $<

clean:
	rm -f *.o