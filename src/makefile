#
# Makefile 
#

# Compiler
cc = g++

# Object files
obj = stdafx.o BiCGStab.o SolidBody.o Calculate_A.o Calculate_B.o CalculateForce.o Calculate_press.o \
      Output.o ImmersedBoundaryMeth.o

ibm: $(obj)
	$(cc) -O2 -o ibm $(obj)
	mv ibm ./../bin

%.o : %.cpp
	$(cc) -std=c++11 -c $<

clean:
	rm -f *.o