emu.out: emu.c main.c pade.c pk_to_xi.c
	gcc -o emu.out -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm -g main.c emu.c pade.c pk_to_xi.c

DeltaSigma.o: DeltaSigma.c gamma_t.h pade.c  emu.c likelihood.c
	gcc -c -g DeltaSigma.c -I/usr/local/include

DeltaSigma.out: DeltaSigma.o likelihood.c
	gcc -o DeltaSigma.out pade.c  emu.c -I/usr/local/include DeltaSigma.o -L/usr/local/lib -lgsl -lgslcblas -lm -g


gamma_t.o: gamma_t.c gamma_t.h pade.c  emu.c 
	gcc -c -g gamma_t.c -I/usr/local/include

gamma_t.out: gamma_t.o
	gcc -o gamma_t.out pade.c  emu.c -I/usr/local/include gamma_t.o -L/usr/local/lib -lgsl -lgslcblas -lm -g

