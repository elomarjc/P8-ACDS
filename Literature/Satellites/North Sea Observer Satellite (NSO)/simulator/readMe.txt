%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NSO Simulator  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


-For automatic setup of path, start matlab in the path of startup.m, which sets up the paths.
-To compile the mex-files just run startup('hep').
-A prerequisite is a preinstallment of Yalmip and SeDuMi.

The simulator source files are sorted in designated directories, which name reflects the purpose
for which they are used. The different directories are 
______________________________________________________________________________________________________
|	Directory Name 	|	Content Description									|
|-----------------------------------------------------------------------------------------------------|
|	root  		|	Files used to ensure simulation on various versions of Matlab, 		|
|				|	the magnetic field model data and startup.m      				|
|	lib			|	Contains a simulink library with the components of the simulator.		|
|     lib/albedo_toolbox| 	Contains a source model for the albedo toolbox used in			|
|      			|	the modeling developed by M.Sc.E.E., Ph.D., Dan Bhanderi, AAU.		|
|     lib/pictures  	|	Contains jpeg-files used in different simulink block masks. 		|
|     lib/utils"  	|	Contains mex- and m-files used in the mathematical modeling of		|
|     			|	the satellite sorted into subdirectories.						|
|     misc 			|	Miscellaneous script files for, e.g , compiling mex-files			|
|     			|	printing simulink models or illustrating vectors.				|
|     design 		|	Files used to create different controllers.					|
|     models 		|	Contains various simulink model files used as templates for test cases. |
|     test 			|	Contains simulink model files with specific controllers			|
|				|	used to evaluate the individual controllers.					|
|_______________________|_____________________________________________________________________________|

Generally a description of the m-files used in the simulator can also be found in the
individual files.