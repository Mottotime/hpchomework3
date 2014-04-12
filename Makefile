
hw3:hw3.f90
	mpif90 -O2 -o run hw3.f90
clean:
	rm run
run:
	mpirun -np 4 -machine mpdhost ./run
