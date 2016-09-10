************************************************************************************************************
*************************         		 PHASM-DB		         ***************************
********************  POLARIZATION FROM HANLE EFFECT FOR STELLAR MAGNETISM - DATABASE  *********************

**********************          	    VIEU Thibault 2016		              **********************
**********************  Instituto de Astrofisica de Canarias / Université Paris-Sud   **********************
**********************               Contact : thibault.vieu@u-psud.fr                **********************
************************************************************************************************************

------------------------------------------------    INFOS    -----------------------------------------------

- PHASM-DB computes a database which contains the Hanle integrated linear polarization signals emitted from a stellar atmosphere in the presence of a global dipolar magnetic field.
- If you want to change the stellar parameters, you need to modify the c++ code, which is not hard since all the global physical constants are gathered at the very beginning.


-----------------------------------------  BEFORE RUNNING THE .EXE  ----------------------------------------

- Make sure that the file "rotation_axis.txt" is with the .exe
- Make sure that the four files containing the tables of the exponential integrals are in the folder "exp_int" which is in the same folder as the .exe


-------------------------------------------------  STEPS  --------------------------------------------------

- You can easily personalize the steps and the limits of the database at the beginning of the main() function in main.cpp
- Default :
    double stepB = 20;
    double steppsi = 90;
    double Bmin=0,Bmax=20,psimin=0,psimax=90;


------------------------------------------------  RESULT  --------------------------------------------------

- The result of PHASM-DB is a folder 'BDD' which contains five folders 'a', 'i', 'p', 'q', 'u', which contain files. Those files are the integrated signatures along the rotation.
- The files are named X_B_PSI.txt, where X is the Stokes parameter (depends on the folder), B the magnetic field strength and PSI the misalignment with respect to the rotation axis.

- You can plot a given signature from its result file.
- You can use the entire database to do an inversion with PHASM-IP.


---------------------------------------------  RECOMPILATION  ----------------------------------------------

- All the code is written in a single main.c++ file in order to make the recompilation easier.
- When recompiling the code, you may have troubles with static arrays, especially on Mac. Try dynamic ones.


---------------------------------------------  ANOTHER STAR  -----------------------------------------------

- You can quite easily modify the code to get the polarization signals from a different line or/and star (Sun by default).
- Make sure this line/star is suitable for our approximations.
- Modify the physical constants in the first part of the c++ file.
- Modify the arrays in the function "hydro".


************************************************************************************************************
------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------