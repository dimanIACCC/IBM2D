The project is devoted to numerical simulation of disperse flows in the channels of various shapes.



============================================================================================
How to compile:
The project required the Boost library (version 1.65 recommended) https://www.boost.org/ 
============================================================================================
How to start:
1)Create the input file that may consist of :
Re 					= 	Reynolds number (integer value)
L 					= 	length of the channel (integer value)
H 					= 	height of the channel (integer value)
N1					= 	number of nodes on X-axis (integer value)
N2 					=	number of nodes on Y-axis (integer value)
d_t 				= 	time step (double value)
Nn 					= 	number of nodes on each solid body (integer value)
rho 				= 	density (double value)
r					= 	radius of solid body (double value)
output_step 		=	frequency of writing output files (integer value)
N_max 				=	maximum of time iterations (integer value)
N_Zeidel			=	maximum of iterations in Zeidel method (integer value)
Zeidel_eps 			=	level of discrepancy in Zeidel method (double value)
InelasticCollision 	=	 (boolean)
k_dist 				=	 1.1
AddSolids_N 		=	how many solids will be added every [AddSolids_interval]
AddSolids_start 	=	number of time-step when the solids addition will be started
AddSolids_interval 	=	time interval of solid addition
BC 					=	border conditions: {u_infinity, u_inflow, periodic}

2)Create Solids.txt(if necessary)

3)Start:
There are several key`s to start program:
 -dir=[path] : read input files such as {input.txt, Solids.txt} from [path] and write all results into [path]
 -h=[path]   : run the program using a hibernation file, please check up the existence of {hibernation.txt} into [path]
				so the program continues the calculation
				
here some examples of correct launch-lines:
ibm.exe -dir=MyResultFolder
ibm.exe -h=D:\MyResultFolder			
==================================================================================================
While program working it creates: 
folder "Solids" 	-	contains files for each solid body (0.plt, 1.plt,...)
"hibernation.txt" 	-	file for continuation calculation
"log.txt" 			- 	file keep information of launch and calculation (such as time of start, date, etc.)
"step%N%.plt" 		-	output file of N-iteration cosist of pressure, velocity, force in calculation area   