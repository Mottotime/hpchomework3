program main
implicit none

include 'mpif.h'

integer nb,tsize,steps,bsize ! num of blocks, total size, steps, block size. Size is num of columns or rows.
integer myid,nproc,rc,ierr 
real,allocatable:: a(:,:),b(:,:),temp1(:),temp2(:)

integer bi,bj ! 2d index of blocks
integer left,right,up,down

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

call time_jacobi(1,bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,myid,ierr)
call time_jacobi(2,bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,myid,ierr)
call time_jacobi(3,bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,myid,ierr)

print *,myid,'Exit!'
call MPI_Finalize(rc)

end program !main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine time_jacobi(comm_type,bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,myid,ierr)
    implicit none
    include 'mpif.h'

    integer:: comm_type,bi,bj,left,right,up,down,bsize,steps,nb,ierr,myid
    real:: a(bsize+2,bsize+2),b(bsize+2,bsize+2),temp1(bsize),temp2(bsize)
    real*8 time1,time2,time,gtime

    call block_value(bsize,a,bi,bj,nb)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    time1=MPI_WTIME()
    
    if(comm_type==1)then
        call block_jacobi(bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,ierr)
    elseif(comm_type==2)then
        call sr_jacobi(bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,ierr)
    else
        call nonblock_jacobi(bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,ierr)
    end if

    time2=MPI_WTIME()
    time=time2-time1
    call MPI_REDUCE(time,gtime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    if(myid==0)then
        print *,gtime,'is the running time of comm type = ',comm_type
   end if
   call flush(6)
end !subroutine time_jacobi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nonblock_jacobi(bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,ierr)
    implicit none
    include 'mpif.h'

    integer:: bi,bj,left,right,up,down,bsize,steps,nb,ierr
    real:: a(bsize+2,bsize+2),b(bsize+2,bsize+2),temp1(bsize),temp2(bsize)

    integer:: begin_col,end_col,begin_row,end_row
    integer n,i,j
    integer req(8)
    integer stat(MPI_STATUS_SIZE,8)
    !integer status(MPI_STATUS_SIZE),status1(MPI_STATUS_SIZE),status2(MPI_STATUS_SIZE)
    integer::tag1=3,tag2=4,tag3=5,tag4=6

    real:: temp3(bsize),temp4(bsize)

    call begin_end(begin_col,begin_row,end_col,end_row,bsize,nb,bi,bj)

    do n=0,steps
        !send data to right

	do i=1,8
	    req(i)=MPI_REQUEST_NULL
	end do

        call MPI_IRECV(a(2,1),bsize,MPI_REAL,left,tag1,&
                MPI_COMM_WORLD,req(1),ierr)
        !call MPI_WAIT(req(1),status,ierr)
        call MPI_ISEND(a(2,bsize+1),bsize,MPI_REAL,right,tag1,&
                MPI_COMM_WORLD,req(2),ierr)
        !call MPI_WAIT(req(2),status,ierr)
	!call MPI_WAITALL(2,req,status,ierr)

        !send data to left
        call MPI_IRECV(a(2,bsize+2),bsize,MPI_REAL,right,tag2,&
                MPI_COMM_WORLD,req(3),ierr)
        call MPI_ISEND(a(2,2),bsize,MPI_REAL,left,tag2,&
                MPI_COMM_WORLD,req(4),ierr)

        !send data up
        do i=1,bsize
            temp1(i)=a(2,i+1)
        end do
        call MPI_IRECV(temp2,bsize,MPI_REAL,down,tag3,&
                MPI_COMM_WORLD,req(5),ierr)
        call MPI_ISEND(temp1,bsize,MPI_REAL,up,tag3,&
                MPI_COMM_WORLD,req(6),ierr)
        do i=1,bsize
            a(bsize+2,i+1)=temp2(i)
        end do

!	call MPI_WAITALL(4,req(1),stat(1,1),ierr)

        !send data down
        do i=1,bsize
            temp3(i)=a(bsize+1,i+1)
        end do
        call MPI_IRECV(temp4,bsize,MPI_REAL,up,tag4,&
                MPI_COMM_WORLD,req(7),ierr)
        call MPI_ISEND(temp3,bsize,MPI_REAL,down,tag4,&
                MPI_COMM_WORLD,req(8),ierr)
        do i=1,bsize
            a(1,i+1)=temp4(i)
        end do

!	call MPI_WAITALL(4,req(5),stat(5,1),ierr)

        !calculation
        do j=begin_col+1,end_col-1
            do i=begin_row+1,end_row-1
                b(i,j)=(a(i,j+1)+a(i,j-1)+a(i+1,j)+a(i-1,j))*0.25
            end do
        end do

        do j=begin_col+2,end_col-2
            do i=begin_row+2,end_row-2
                a(i,j)=b(i,j)
            end do
        end do

	call MPI_WAITALL(8,req,stat,ierr)

	do i=begin_row,end_row
	    b(i,begin_col)=(a(i,begin_col-1)+a(i,begin_col+1)+a(i+1,begin_col)+a(i-1,begin_col))*0.25
	end do
	do i=begin_row,end_row
	    b(i,end_col)=(a(i,end_col-1)+a(i,end_col+1)+a(i+1,end_col)+a(i-1,end_col))*0.25
	end do
	do j=begin_col,end_col
	    b(begin_row,j)=(a(begin_row,j-1)+a(begin_row,j+1)+a(begin_row+1,j)+a(begin_row-1,j))*0.25
	end do
	do j=begin_col,end_col
	    b(end_row,j)=(a(end_row,j-1)+a(end_row,j+1)+a(end_row+1,j)+a(end_row-1,j))*0.25
	end do

	do i=begin_row,end_row
	    do j=begin_col,begin_col+1
		a(i,j)=b(i,j)
	    end do
	end do
	do i=begin_row,end_row
	    do j=end_col,end_col-1
		a(i,j)=b(i,j)
	    end do
	end do
	do j=begin_col+2,end_col-2
	    a(begin_row,j)=b(begin_row,j)
	    a(begin_row+1,j)=b(begin_row+1,j)
	    a(end_row-1,j)=b(end_row-1,j)
	    a(end_row,j)=b(end_row,j)
	end do
    end do
end ! subroutine nonblock_jacobi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine block_jacobi(bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,ierr)
    implicit none
    include 'mpif.h'

    integer:: bi,bj,left,right,up,down,bsize,steps,nb,ierr
    real:: a(bsize+2,bsize+2),b(bsize+2,bsize+2),temp1(bsize),temp2(bsize)

    integer begin_col,end_col,begin_row,end_row
    integer n,i,j
    integer status(MPI_STATUS_SIZE)
    integer::tag1=3,tag2=4,tag3=5,tag4=6

    call begin_end(begin_col,begin_row,end_col,end_row,bsize,nb,bi,bj)

    do n=0,steps
        !send data to right
        call MPI_RECV(a(2,1),bsize,MPI_REAL,left,tag1,&
                MPI_COMM_WORLD,status,ierr)
        call MPI_SEND(a(2,bsize+1),bsize,MPI_REAL,right,tag1,&
                MPI_COMM_WORLD,ierr)

        call MPI_RECV(a(2,bsize+2),bsize,MPI_REAL,right,tag2,&
                MPI_COMM_WORLD,status,ierr)
        call MPI_SEND(a(2,2),bsize,MPI_REAL,left,tag2,&
                MPI_COMM_WORLD,ierr)

        do i=1,bsize
            temp1(i)=a(2,i+1)
        end do
        call MPI_RECV(temp2,bsize,MPI_REAL,down,tag3,&
                MPI_COMM_WORLD,status,ierr)
        call MPI_SEND(temp1,bsize,MPI_REAL,up,tag3,&
                MPI_COMM_WORLD,ierr)
        do i=1,bsize
            a(bsize+2,i+1)=temp2(i)
        end do

        !send data down
        do i=1,bsize
            temp1(i)=a(bsize+1,i+1)
        end do
        call MPI_RECV(temp2,bsize,MPI_REAL,up,tag4,&
                MPI_COMM_WORLD,status,ierr)
        call MPI_SEND(temp1,bsize,MPI_REAL,down,tag4,&
                MPI_COMM_WORLD,ierr)
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
end ! subroutine block_jacobi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sr_jacobi(bi,bj,left,right,up,down,a,b,temp1,temp2,bsize,steps,nb,ierr)
    implicit none
    include 'mpif.h'

    integer:: bi,bj,left,right,up,down,bsize,steps,nb,ierr
    real:: a(bsize+2,bsize+2),b(bsize+2,bsize+2),temp1(bsize),temp2(bsize)

    integer begin_col,end_col,begin_row,end_row
    integer n,i,j
    integer status(MPI_STATUS_SIZE)
    integer::tag1=3,tag2=4,tag3=5,tag4=6

    call begin_end(begin_col,begin_row,end_col,end_row,bsize,nb,bi,bj)

    do n=0,steps

        call MPI_SENDRECV(a(2,bsize+1),bsize,MPI_REAL,right,tag1,&
               a(2,1),bsize,MPI_REAL,left,tag1,&
               MPI_COMM_WORLD,status,ierr)

        call MPI_SENDRECV(a(2,2),bsize,MPI_REAL,left,tag2,&
               a(2,bsize+2),bsize,MPI_REAL,right,tag2,&
               MPI_COMM_WORLD,status,ierr)

        do i=1,bsize
            temp1(i)=a(2,i+1)
        end do
        call MPI_SENDRECV(temp1,bsize,MPI_REAL,up,tag3,&
                temp2,bsize,MPI_REAL,down,tag3,&
                MPI_COMM_WORLD,status,ierr)
        do i=1,bsize
        a(bsize+2,i+1)=temp2(i)
        end do

        do i=1,bsize
            temp1(i)=a(bsize+1,i+1)
        end do
        call MPI_SENDRECV(temp1,bsize,MPI_REAL,down,tag4,&
               temp2,bsize,MPI_REAL,up,tag4,&
               MPI_COMM_WORLD,status,ierr)
        do i=1,bsize
        a(1,i+1)=temp2(i)
        end do

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
end ! subroutine sr_jacobi

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

subroutine begin_end(begin_col,begin_row,end_col,end_row,bsize,nb,bi,bj)
    implicit none
    integer::begin_col,begin_row,end_col,end_row,bsize,nb,bi,bj
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
end !subroutine begin_end
