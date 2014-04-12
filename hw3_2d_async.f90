program main
implicit none

include 'mpif.h'
include 'utils.f90'

integer nb,tsize,steps,bsize ! num of blocks, total size, steps, block size. Size is num of columns or rows.
integer myid,nproc,rc,ierr 
real,allocatable:: a(:,:),b(:,:),temp1(:),temp2(:)

integer bi,bj ! 2d index of blocks
integer left,right,up,down

real*8 time1,time2,time,gtime

call start_mpi(myid,nproc,tsize,nb,steps,ierr,bsize)

allocate(a(bsize+2,bsize+2),stat=ierr)
allocate(b(bsize+2,bsize+2),stat=ierr)
allocate(temp1(bsize))
allocate(temp2(bsize))
if (ierr/=0) then
    print *,"Unable to allocate."
    call flush(6)
    stop 
end if

bi=myid/nb
bj=MOD(myid,nb)

call block_neighbour(bi,bj,myid,left,right,up,down,nb)
call block_value(bsize,a,bi,bj,nb)

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
time1=MPI_WTIME()
call nonblock_jacobi(bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,ierr)
time2=MPI_WTIME()
time=time2-time1
call MPI_REDUCE(time,gtime,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)

if(myid==0)then
    print *,gtime,'is the running time.'
end if

print *,myid,'Exit!'
call MPI_Finalize(rc)

end program !main

subroutine nonblock_jacobi(bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,ierr)
    implicit none
    include 'mpif.h'

    integer:: bi,bj,left,right,up,down,bsize,steps,nb,ierr
    real:: a(bsize+2,bsize+2),b(bsize+2,bsize+2),temp1(bsize),temp2(bsize)

    integer begin_col,end_col,begin_row,end_row
    integer n,i,j
    integer isend,irecv
    integer status(MPI_STATUS_SIZE),status1(MPI_STATUS_SIZE),status2(MPI_STATUS_SIZE)
    integer::tag1=3,tag2=4,tag3=5,tag4=6

    begin_col=2
    end_col=bsize+1
    begin_row=2
    end_row=bsize+1

    if(bj==0)then
    begin_col=3
    end if
    if(bj==nb-1)then
    end_col=bsize
    endif
    if(bi==0)then
    begin_row=3
    end if
    if(bi==nb-1)then
    end_row=bsize
    end if

    do n=0,steps
        !send data to right
        if(bj>0)then
        call MPI_IRECV(a(2,1),bsize,MPI_REAL,left,tag1,&
                MPI_COMM_WORLD,irecv,ierr)
        call MPI_WAIT(irecv,status,ierr)
        end if
        if(bj<nb-1)then
        call MPI_ISEND(a(2,bsize+1),bsize,MPI_REAL,right,tag1,&
                MPI_COMM_WORLD,isend,ierr)
        call MPI_WAIT(isend,status,ierr)
        end if
        !send data to left
        if(bj<nb-1)then
        call MPI_RECV(a(2,bsize+2),bsize,MPI_REAL,right,tag2,&
                MPI_COMM_WORLD,status,ierr)
        end if
        if(bj>0)then
        call MPI_SEND(a(2,2),bsize,MPI_REAL,left,tag2,&
                MPI_COMM_WORLD,ierr)
        end if
        do i=1,bsize
            temp1(i)=a(2,i+1)
        end do

        !send data up
        if(bi<nb-1)then
        call MPI_RECV(temp2,bsize,MPI_REAL,down,tag3,&
                MPI_COMM_WORLD,status,ierr)
        end if
        if(bi>0)then
        call MPI_SEND(temp1,bsize,MPI_REAL,up,tag3,&
                MPI_COMM_WORLD,ierr)
        end if
        do i=1,bsize
            a(bsize+2,i+1)=temp2(i)
        end do

        !send data down
        do i=1,bsize
            temp1(i)=a(bsize+1,i+1)
        end do
        if(bi>0)then
        call MPI_RECV(temp2,bsize,MPI_REAL,up,tag4,&
                MPI_COMM_WORLD,status,ierr)
        end if
        if(bi<nb-1)then
        call MPI_SEND(temp1,bsize,MPI_REAL,down,tag4,&
                MPI_COMM_WORLD,ierr)
        end if
        do i=1,bsize
            a(1,i+1)=temp2(i)
        end do

        !calculation
        do j=begin_col,end_col
            do i=begin_row,end_row
                b(i,j)=(a(i,j+1)+a(i,j-1)+a(i+1,j)+a(i-1,j))*0.25
            end do
        end do

        do j=begin_col,end_col
            do i=begin_row,end_row
                a(i,j)=b(i,j)
            end do
        end do

    end do
end ! subroutine nonblock_jacobi
