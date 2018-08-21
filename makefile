TARGET = LbmCavity3D
USER = $(shell whoami)
CC = sw5cc
LD = mpicc 

CFLAGS =  -O3 -host -I/usr/sw-mpp/mpi2/include/ -lm

OBJ = LbmCavity3D.o Parallel.o MasterController.o MasterCollide.o MasterStream.o SlaveController.o SlaveStream.o SlaveCollide.o HchTimer.o

LIB = lib/liblbm.a

$(TARGET): $(OBJ)
	$(LD) $(OBJ) $(LIB) -o $(TARGET) 
	# rm $(OBJ)

SlaveController.o: SlaveController.c
	sw5cc -c -slave $^

%.o:%.c
	$(CC) $(CFLAGS) -c $<

run:$(TARGET)
	bsub -I -b -q q_sw_cpc_2 -cgsp 64 -n 16 -np 4  -share_size 6500 -host_stack 500 -J test ./LbmCavity3D $(USER)

#-------------------------------------*
.PHONY : clean clear
clean:
	-rm -rf $(TARGET) $(OBJ) 
	
