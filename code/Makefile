######################################################################################
CC = gfortran

#OPTIONS = -O3 -ffinite-math-only -fno-trapping-math -fPIC -flto -march=native \
		  -finit-local-zero -fno-strict-aliasing -fopenmp -fopenmp-simd

OPTIONS = -O3 -fPIC -finit-local-zero -mcmodel=medium

#OPTIONS_PY = --opt="-O3 -ffinite-math-only -fno-trapping-math -fPIC -flto -march=native \
					-finit-local-zero -fno-strict-aliasing -fopenmp -fopenmp-simd"

OPTIONS_PY = --opt="-O3 -fPIC -finit-local-zero -mcmodel=medium"

######################################################################################
## With OpenMP on Mac:

#OPTIONS_PAR_PY = --opt="-fopenmp -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000  -lgomp -lpthread "  -lgomp -lpthread
#OPTIONS_PAR    =        -fopenmp -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000  -lgomp -lpthread

#OPTIONS_PAR_PY =
#OPTIONS_PAR =

######################################################################################
## With OpenMP on Linux:
#OPTIONS_PAR = -O3 -w -finit-local-zero -fopenmp -lgomp -lpthread
#OPTIONS_PAR = -fopenmp
OPTIONS_PAR_PY = --opt="-fopenmp -lgomp -lpthread"
#OPTIONS_PAR_PY = --opt="-fopenmp"
## Run these two command in the terminal (or put them in the ALIAS):
## export OMP_STACKSIZE=1024M
## ulimit -s 1024000
## OMP_STACKSIZE works for slave threads
## ulimit works for the master thread


MODULES = #-J./Modules

Trayectoria_Py: main_2805.for odeint.o rkqs.o rkck.o ultimo.o polint.o midsqu.o OMP_Void.o
	f2py --fcompiler=gfortran $(OPTIONS_PY) -c -m Trayectoria_Py main_2805.for odeint.o rkqs.o rkck.o \
	OMP_Void.o ultimo.o polint.o midsqu.o 

Trayectoria_Py_P: main_2805.for odeint.o rkqs.o rkck.o ultimo.o polint.o midsqu.o 
	f2py --fcompiler=gfortran $(OPTIONS_PY) $(OPTIONS_PAR_PY) -c -m Trayectoria_Py_P main_2805.for odeint.o rkqs.o rkck.o \
	ultimo.o polint.o midsqu.o 
#	-lgomp -lpthread
#	-L/opt/local/lib -lgomp -lpthread

qromotest.out: mainqromo.for ultimo.o polint.o midsqu.o
	gfortran $(OPTIONS) $(OPTIONS_PAR) -o qromotest.out mainqromo.for \
	ultimo.o polint.o midsqu.o

Trayectoria.out: main_2805.for odeint.o rkqs.o rkck.o ultimo.o polint.o midsqu.o
	gfortran $(OPTIONS) $(OPTIONS_PAR) -o Trayectoria.out main_2805.for odeint.o rkqs.o rkck.o \
	ultimo.o polint.o midsqu.o

#Trayectoria.out: main_ultpts.for
#	gfortran $(OPTIONS) $(OPTIONS_PAR) -o Trayectoria.out main_ultpts.for

#Compute.o: Compute.for
#	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c Compute.o

OMP_Void.o: OMP_Void.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c OMP_Void.for

odeint0.o: odeint0.for
	$(CC) $(OPTIONS) -c odeint0.for

odeint.o: odeint.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c odeint.for

rkqs.o: rkqs.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c rkqs.for

rkqs0.o: rkqs0.for
	$(CC) $(OPTIONS) -c rkqs0.for

rkck.o: rkck.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c rkck.for

rkck0.o: rkck0.for
	$(CC) $(OPTIONS) -c rkck0.for

ultimo.o: ultimo.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c ultimo.for

ultimo0.o: ultimo0.for
	$(CC) $(OPTIONS) -c ultimo0.for

polint.o: polint.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c polint.for

polint0.o: polint0.for
	$(CC) $(OPTIONS) -c polint0.for

midsqu.o: midsqu.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c midsqu.for

midsqu0.o: midsqu0.for
	$(CC) $(OPTIONS) -c midsqu0.for

clean: 
	rm -f *o *out *so
