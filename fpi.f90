!---------------------------------------------c
!  Exercise: Pi                               c
!                                             c
!  Compute the value of PI using the integral c
!  pi = 4* int 1/(1+x*x)    x in [0-1]        c
!                                             c
!  The integral is approximated by a sum of   c
!  n interval.                                c
!                                             c
!  The approximation to the integral in each  c
!  interval is: (1/n)*4/(1+x*x).              c
!---------------------------------------------c
program pigreco
    use mpi
    implicit none

    integer(selected_int_kind(18)) :: i, istart, iend
    integer(selected_int_kind(18)), parameter :: intervals=1e7
    integer:: ierr, nproc, myrank
    real(kind(1.d0)) :: dx,sum,x, total_sum
    real(kind(1.d0)) :: f,pi
    real(kind(1.d0)), parameter :: PI25DT = acos(-1.d0)
    real(kind(1.d0)) :: time1, time2

!mpi stuff
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

    if(myrank == 0) then
       write(*,*) "MPI version with tastks = ", nproc
       write(*,*) "Number of intervals     = ", intervals
    endif
    if(mod(intervals,nproc) /= 0) then
       if(myrank == 0) then
          write(*,*) "The number of process must divide", intervals, "exactly."
       endif
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
       stop
    endif

    sum=0.d0
    dx=1.d0/intervals
    time1 = MPI_WTIME()
    istart = (intervals/nproc)*myrank + 1
    iend   = (intervals/nproc)*(myrank+1)
    sum = 0.0;
    total_sum = 0.0;
    do i=iend, istart, -1
        x=dx*(i-0.5d0)
        f=4.d0/(1.d0+x*x)
        sum=sum+f
    end do
    call MPI_Reduce(sum,total_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    time2 = MPI_WTIME()
    pi=dx*total_sum
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(myrank == 0) then
       PRINT '(a13,2x,f30.25)',' Computed PI =', pi
       PRINT '(a13,2x,f30.25)',' The True PI =', PI25DT
       PRINT '(a13,2x,f30.25)',' Error        ', PI25DT-pi
       PRINT *, ' '
       PRINT *, 'Elapsed time ', time2-time1 ,' s'
    endif

    call MPI_FINALIZE(ierr)

end program

