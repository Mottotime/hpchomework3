hw3:hw3.f90
	mpif90 -o run hw3.f90
run:hw3
	mpirun -np 4 ./run