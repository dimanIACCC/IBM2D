#
# Makefile 
#

# Compiler
cc = g++

# Object files
obj = stdafx.o BiCGStab.o GeomVec.o String.o Parameters.o SolidBody.o PredictVel.o \
      Calculate_press.o CalculateForce.o Output.o ImmersedBoundaryMeth.o

ibm: $(obj)
	$(cc) -O2 -o ibm $(obj)
	mv ibm ./../bin

%.o : %.cpp
	$(cc) -std=c++11 -c $<

clean:
	rm -f *.o