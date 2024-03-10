######################################################################################
CC = gfortran

OPTIONS = -O3 -fPIC -finit-local-zero -mcmodel=medium

OPTIONS_PY = --opt="-O3 -fPIC -finit-local-zero -mcmodel=medium"

######################################################################################

MODULES = #-J./Modules

Trayectoria_Py_P: main.for odeint.o rkqs.o rkck.o ultimo.o polint.o midsqu.o 
	f2py --fcompiler=gfortran $(OPTIONS_PY) $(OPTIONS_PAR_PY) -c -m Trayectoria_Py_P main.for odeint.o rkqs.o rkck.o \
	ultimo.o polint.o midsqu.o 
#	-lgomp -lpthread
#	-L/opt/local/lib -lgomp -lpthread

odeint.o: odeint.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c odeint.for

rkqs.o: rkqs.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c rkqs.for

rkck.o: rkck.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c rkck.for

ultimo.o: ultimo.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c ultimo.for

polint.o: polint.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c polint.for

midsqu.o: midsqu.for
	$(CC) $(OPTIONS) $(OPTIONS_PAR) -c midsqu.for

clean: 
	rm -f *o *out *so
