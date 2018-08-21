TARGET = LbmCavity3D
USER = $(shell whoami)
DATE = `date +"%Y-%m-%d_%H:%M"`
CC = sw5cc
LD = mpicc 

CFLAGS =  -O3 -host -I/usr/sw-mpp/mpi2/include/ -lm -msimd -OPT:ieee_arith=1 -fno-math-errno
SCFLAGS = -O3 -slave -msimd -OPT:ieee_arith=1 -lm_slave -fno-math-errno
LDFLAGS = -O3 -msimd -OPT:ieee_arith=1 -lm_slave -fno-math-errno

COBJ = LbmCavity3D.o Collide.o Parallel.o  Stream.o LbmMaster.o HchTimer.o
SOBJ = LbmSlave.o

LIB = lib/liblbm.a

VPATH = ./header/:./master/:./slave/

$(TARGET): $(COBJ) $(SOBJ)
	$(LD) $(COBJ) $(SOBJ) $(LIB) -o $(TARGET) 
#	rm $(COBJ) $(SOBJ)
	
$(COBJ): %.o: %.c
	$(CC) $(CFLAGS) -c $<

$(SOBJ): %.o: %.c
	$(CC) $(SCFLAGS) -c $<

run:
	make
	mkdir ./code-$(DATE)
	cp -r ./header ./code-$(DATE)
	cp -r ./master ./code-$(DATE)
	cp -r ./slave ./code-$(DATE)
	cp -r makefile ./code-$(DATE)
	bsub -I -b -q q_sw_cpc_2 -cgsp 64 -n 16 -np 4  -share_size 6500 -host_stack 500 -J test -o run-$(DATE).log ./LbmCavity3D $(USER)

#-------------------------------------*
.PHONY : clean clear
clean:
	-rm -rf $(COBJ) $(SOBJ)
	
