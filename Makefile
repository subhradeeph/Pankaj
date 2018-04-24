CC = gcc

CFLAGS1 = -std=c99 -g

all: SpinDe

SpinDe:	1D_FFT_Combined_CODE.o get_input.o init_conf.o evolve.o free_energy_density.o out_conf.o interfacial_energy.o
	$(CC) $(CFLAGS1) 1D_FFT_Combined_CODE.o get_input.o init_conf.o evolve.o free_energy_density.o out_conf.o interfacial_energy.o -lfftw3 -lm -o SpinDe

1D_FFT_Combined_CODE.o:	1D_FFT_Combined_CODE.c
			$(CC) $(CFLAGS1) -c 1D_FFT_Combined_CODE.c

%.o:	%.c
	$(CC) $(CFLAGS1) -c $<	

#init_conf.o:	init_conf.c
#	$(CC) $(CFLAGS1) -c init_conf.c

#evolve.o:	evolve.c
#	$(CC) $(CFLAGS1) -c evolve.c

#out_conf.o:	out_conf.c
#	$(CC) $(CFLAGS1) -c out_conf.c

#interfacial_energy.o:	interfacial_energy.c
#	$(CC) $(CFLAGS1) -c interfacial_energy.c

clean:	
	rm -rf *o SpinDe
	rm -rf conf.* 
