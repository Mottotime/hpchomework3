all:hw3_1d hw3_2d_s hw3_2d_a hw3_2d_s_sr

hw3_1d:hw3_1d.f90
	mpif90 -o run_1d hw3_1d.f90

hw3_2d_s:hw3_2d_sync.f90
	mpif90 -o run_2d_s hw3_2d_sync.f90

hw3_2d_s_sr:hw3_2d_sync_sr.f90
	mpif90 -o run_2d_s hw3_2d_sync_sr.f90

hw3_2d_a:hw3_2d_async.f90
	mpif90 -o run_2d_a hw3_2d_async.f90

run:
	mpirun -np 4 -machine mpdhost ./run
