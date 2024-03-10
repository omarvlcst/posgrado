<h3> Dependencies and Requirements for running the code w/o trouble: </h3>

- gfortran in a version of at least 8 (Working in a future version with the Fortran compiler of IBM, coming soon)
- cmake
- Don't overrun in too many points for the integrator. In the file called "*size.for*", customize the array of points for the variables *xd* and *thetad*, so that they don't exceed the amount of memory left in your computer.
- This code is thought for being compiled and executed in a Linux terminal (tested successfully in Ubuntu 22, gcc and gfortran version 9). 

The Makefile makefile includes all the instructions for GFORTRAN compiling and creating a Python module called: *"Trayectoria_Py_P"* that will pass the Fortran subroutines into functions reusable inside a Python framework and environment. An IPython notebook (Jupyter notebook) is included with the routines being cast in each cell where needed, and the module being called at the beginning. This notebook can be executed in a friendly environment, like Jupyter Notebook or a similar code editor able to read such a notebook like for example, Visual Studio Code.

An additional directory called "subroutines" includes all the Fortran subroutines that are used in the code, and a respective Makefile for generating *"Trayectoria.out"* as an executable program in Fortran.

<h3>For compiling  all the .for files in Fortran and end up converting the wrapped f2py module into Python as an executable, simply type down in a Linux terminal the command: </h3>

~$ chmod 777 *.for

~$ make Trayectoria_Py_P

<h3> After editing out any piece of the code and trying to recompile again with the newest version of the files, first clean the memory by typing: </h3>

~$ make clean

And then, the instruction "*make Trayectoria_Py_P*" will recompile and execute the f2py module succesfully in its latest version.
