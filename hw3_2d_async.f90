program main
implicit none

include 'mpif.h'

integer nb,tsize,steps,bsize ! num of blocks, total size, steps, block size. Size is num of columns or rows.
integer myid,nproc,rc,ierr 
real,allocatable:: a(:,:),b(:,:),temp1(:),temp2(:)

integer bi,bj ! 2d index of blocks
integer left,right,up,down

real*8 time1,time2,time,gtime

call start_mpi(myid,nproc,tsize,nb,steps,ierr,bsize)
print *,myid,"1"
call flush(6)

allocate(a(bsize+2,bsize+2),stat=ierr)
allocate(b(bsize+2,bsize+2),stat=ierr)
allocate(temp1(bsize))
allocate(temp2(bsize))
if (ierr/=0) then
    print *,"Unable to allocate."
    call flush(6)
    stop 
end if

print *,myid,"2"
call flush(6)

bi=myid/nb
bj=MOD(myid,nb)

call block_neighbour(bi,bj,myid,left,right,up,down,nb)
print *,myid,"3"
call flush(6)

call block_value(bsize,a,bi,bj,nb)
print *,myid,"4"
call flush(6)

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
time1=MPI_WTIME()
print *,myid,"5"
call flush(6)

call nonblock_jacobi(bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,ierr)

print *,myid,"6"
call flush(6)
time2=MPI_WTIME()
time=time2-time1
call MPI_REDUCE(time,gtime,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
print *,myid,"8"
if(myid==0)then
    print *,gtime,'is the running time.'
end if

print *,myid,'Exit!'
call MPI_Finalize(rc)

end program !main

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

    print *,myid,"11"
    if (nb>0) then
        bsize=tsize/nb
    else
        print *,'number of blocks < 0'
        stop
    end if

    if(nproc/=nb*nb)then
        print *,'Proc num /= num of block ^2'
    end if

    print *, "I am proc (",myid,"), named (",myname,"), in",nproc,"processors."
    call flush(6)
    return
end !subroutine start_mpi

!subroutine ab_allocation(a,b,bsize,ierr)
!    implicit none
!
!    integer:: bsize,ierr
!    real,allocatable:: a(:,:),b(:,:),temp1(:),temp2(:)
!
!    allocate(a(bsize+2,bsize+2),stat=ierr)
!    if(ierr/=0)then
!        print *,'Unsuccessful allocation!'
!        return
!    end if
!    allocate(b(bsize+2,bsize+2),stat=ierr)
!    if(ierr/=0)then
!        print *,'Unsuccessful allocation!'
!        return
!    end if
!    allocate(temp1(bsize))
!    if(ierr/=0)then
!        print *,'Unsuccessful allocation!'
!        return
!    end if
!    allocate(temp2(bsize))
!    if(ierr/=0)then
!        print *,'Unsuccessful allocation!'
!        return
!    end if
!end !subroutine ab_allocation

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
