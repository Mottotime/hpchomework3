subroutine start_mpi(myid,nproc,tsize,nb,steps,ierr,bsize)
    implicit none
    include 'mpif.h'

    integer:: myid,nproc,tsize,nb,steps,ierr,bsize

    integer namelen
    character(len=MPI_MAX_PROCESSOR_NAME) myname
    namelist /para/ tsize,nb,steps

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    call MPI_GET_PROCESSOR_NAME(myname,namelen,ierr)

    if (myid == 0)then
        open(10,file="namelist")
        read(10,nml=para)
        print*,tsize,nb,steps
    end if
    call MPI_BCAST(tsize,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(steps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if (nb>0) then
        bsize=tsize/nb
    else
        print *,'number of blocks < 0'
        stop
    end if

    if(nproc/=nb*nb)then
        print *,'Proc num /= num of block ^2'
    end if

    print *, "I am proc (",myid,"), named (",myname
    call flush(6)
    return
end !subroutine start_mpi

subroutine block_neighbour(bi,bj,myid,left,right,up,down,nb)
    implicit none
    include 'mpif.h'
    integer:: bi,bj,myid,left,right,up,down,nb

    if (bj>0)then
    left=myid-1
    else
        left=MPI_PROC_NULL
    end if
    if (bj<nb-1)then
        right=myid+1
    else
        right=MPI_PROC_NULL
    end if
    if (bi<nb-1)then
        down=myid+nb
    else
        down=MPI_PROC_NULL
    end if
    if (bi>0)then
        up=myid-nb
    else
        up=MPI_PROC_NULL
    end if

    !print *,myid,bi,bj,left,right,up,down
    !call flush(6)
end !subroutine block_neighbour

subroutine block_value(bsize,a,bi,bj,nb)
    implicit none

    integer::bsize,bi,bj,nb
    real::a(bsize+2,bsize+2)

    integer i,j
    do j=1,bsize+2
        do i=1,bsize+2
        a(i,j)=0.0
        end do
    end do

    if (bj==0)then
        do i=2,bsize+1
        a(i,2)=8.0
        end do
    end if

    if (bj==nb-1)then
        do i=2,bsize+1
        a(i,bsize+1)=8.0
        end do
    end if

    if(bi==0)then
        do i=2,bsize+1
        a(1,i)=8.0
        end do
    end if

    if(bi==nb-1)then
        do i=2,bsize+1
        a(bsize+1,i)=8.0
        end do
    end if
end !subroutine block_value