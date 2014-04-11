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

integer n,myid,nproc,i,j,rc,namelen
integer bi,bj ! 2d index of blocks
character(len=MPI_MAX_PROCESSOR_NAME) myname
!real a(tsize,bsize+2),b(tsize,bsize+2)

integer begin_col,end_col,ierr
integer status(MPI_STATUS_SIZE)
integer left,right,tag1,tag2
tag1=3
tag2=4


call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
call MPI_GET_PROCESSOR_NAME(myname,namelen,ierr)

if (myid == 0)then
    read(*,nml=para)
    print*,tsize,nb,steps
    if (nb>0) then
	bsize=tsize/nb
    else
	print *,'number of blocks < 0'
	stop
    end if
end if
call MPI_BCAST(tsize,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(nb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(steps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

if(nproc/=nb*nb)then
    print *,'Proc num /= num of block ^2'
end if

print *, "I am proc (",myid,"), named (",myname,"), in",nproc,"processors."
call flush(6)

bi=myid/nb
bj=MOD(myid,nb)

print *,myid,bi,bj

if (myid>0)then
    left=myid-1
else
    left=MPI_PROC_NULL
end if
if (myid<3)then
    right=myid+1
else
    right=MPI_PROC_NULL
end if


!do j=1,bsize+2
!    do i=1,tsize
!	a(i,j)=0.0
!    end do
!end do
!
!if (myid==0)then
!    do i=1,tsize
!	a(i,2)=8.0
!    end do
!end if
!
!if (myid==3)then
!    do i=1,tsize
!	a(i,bsize+1)=8.0
!    end do
!end if
!
!do i=1,bsize+2
!    a(1,i)=8.0
!    a(tsize,i)=8.0
!end do
!
!do n=1,steps
!    call MPI_SENDRECV(a(1,bsize+1),tsize,MPI_REAL,right,tag1,&
!			a(1,1),tsize,MPI_REAL,left,tag1,&
!			MPI_COMM_WORLD,status,ierr)
!    call MPI_SENDRECV(a(1,2),tsize,MPI_REAL,left,tag2,&
!			a(1,bsize+2),tsize,MPI_REAL,right,tag2,&
!			MPI_COMM_WORLD,status,ierr)
!    begin_col=2
!    end_col=bsize+1
!
!    if(myid==0)then
!	begin_col=3
!    endif
!    if(myid==3)then
!	end_col=bsize
!    endif
!    do j=begin_col,end_col
!	do i=2,tsize-1
!	    b(i,j)=(a(i,j+1)+a(i,j-1)+a(i+1,j)+a(i-1,j))*0.25
!	end do
!    end do
!    do j=begin_col,end_col
!	do i=2,tsize-1
!	    a(i,j)=b(i,j)
!	end do
!    end do
!end do
!
!do i=2,tsize-1
!    print *,myid,(a(i,j),j=begin_col,end_col)
!end do
call MPI_Finalize(rc)

end program !main
