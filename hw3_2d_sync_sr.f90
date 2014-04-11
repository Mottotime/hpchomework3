program main
!use mpi
implicit none

include 'mpif.h'

!integer,parameter:: N=9
!integer,parameter:: tsize=1024
!integer,parameter:: bsize=tsize/N
!integer,parameter:: steps=10
integer bsize,nb,tsize,steps
namelist /para/ tsize,nb,steps
real time1,time2,time,gtime

integer n,myid,nproc,i,j,rc,namelen
integer bi,bj ! 2d index of blocks
character(len=MPI_MAX_PROCESSOR_NAME) myname
!real a(tsize,bsize+2),b(tsize,bsize+2)
real,allocatable:: a(:,:),b(:,:),temp1(:),temp2(:)

integer begin_col,end_col,begin_row,end_row,ierr
integer status(MPI_STATUS_SIZE)
integer left,right,tag1,tag2,up,down,tag3,tag4
integer isend1,isend2,irecv1,irecv2,istat1,istat2,istat3,istat4
tag1=3
tag2=4
tag3=5
tag4=6

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
call MPI_GET_PROCESSOR_NAME(myname,namelen,ierr)

if (myid == 0)then
    read(*,nml=para)
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

!print *, "I am proc (",myid,"), named (",myname,"), in",nproc,"processors."
!call flush(6)

allocate(a(bsize+2,bsize+2),stat=ierr)
allocate(b(bsize+2,bsize+2),stat=ierr)
if(ierr/=0)then
    print *,'Unsuccessful allocation!'
    stop
end if
allocate(temp1(bsize))
allocate(temp2(bsize))

bi=myid/nb
bj=MOD(myid,nb)


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

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
time1=MPI_WTIME()

do n=0,steps
!    call MPI_SENDRECV(a(2,bsize+1),bsize,MPI_REAL,right,tag1,&
!			a(2,1),bsize,MPI_REAL,left,tag1,&
!			MPI_COMM_WORLD,status,ierr)
    if(bj>0)then
	call MPI_RECV(a(2,1),bsize,MPI_REAL,left,tag1,&
		    MPI_COMM_WORLD,status,ierr)
    end if
    if(bj<nb-1)then
	call MPI_SEND(a(2,bsize+1),bsize,MPI_REAL,right,tag1,&
			MPI_COMM_WORLD,ierr)
    end if
!    call MPI_SENDRECV(a(2,2),bsize,MPI_REAL,left,tag2,&
!			a(2,bsize+2),bsize,MPI_REAL,right,tag2,&
!			MPI_COMM_WORLD,status,ierr)
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
!    call MPI_SENDRECV(temp1,bsize,MPI_REAL,up,tag3,&
!			temp2,bsize,MPI_REAL,down,tag3,&
!			MPI_COMM_WORLD,status,ierr)
    
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

    do i=1,bsize
	temp1(i)=a(bsize+1,i+1)
    end do
!    call MPI_SENDRECV(temp1,bsize,MPI_REAL,down,tag4,&
!			temp2,bsize,MPI_REAL,up,tag4,&
!			MPI_COMM_WORLD,status,ierr)
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
!


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
time2=MPI_WTIME()
time=time2-time1
!print *,myid,time
call MPI_REDUCE(time,gtime,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
if(myid==0)then
    print *,gtime,'is the running time.'
end if
!
!do i=2,tsize-1
!    print *,myid,(a(i,j),j=begin_col,end_col)
!end do
print *,myid,'Exit!'
call MPI_Finalize(rc)

end program !main
