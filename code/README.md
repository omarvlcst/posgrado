This code is thought for use in Linux (tested and successful in Ubuntu 22, gcc and gfortran version 9). 

The Makefile makefile includes all the instructions for gcc compiling and creating a Python module called: *"Trayectoria_Py_P"* that will pass the Fortran subroutines into functions reusable inside a Python framework and environment.

An additional directory called "subroutines" includes all the Fortran subroutines that are used in the code, and a respective Makefile for generating *"Trayectoria.out"* as an executable program in Fortran.

<h3>For converting the wrapped f2py module into Python, simply type down in a Linux terminal the command: </h3>

~$ chmod 777 *.for

~$ make Trayectoria_Py_P
