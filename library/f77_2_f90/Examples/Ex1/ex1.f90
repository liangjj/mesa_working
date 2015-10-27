  PROGRAM test
! Example to allocate an array, pass it to a subroutine
! where it appears as a real and integer variable, fill
! it with some values and print.
! An explicit interface in the main program is not needed
! because the arrays in the subroutine are dimensioned.
  IMPLICIT NONE
  INTEGER  :: inp=5, iout=6   
  INTEGER n
  INTEGER, DIMENSION(:), ALLOCATABLE ::igrid
  INTEGER, DIMENSION(:), ALLOCATABLE ::jgrid
  open (inp,file='inp',status='old')
  open (iout,file='out',status='unknown')  
  read(inp,*) n
  ALLOCATE(igrid(2*n),jgrid(2*n*n))
  call mytest(igrid,igrid,jgrid,jgrid,n,iout)
  DEALLOCATE(igrid,jgrid)
  stop
  END PROGRAM test
  SUBROUTINE mytest(gridi,igrid,gridj,jgrid,n,iout)
  REAL*8, DIMENSION(n) :: gridi
  REAL*8, DIMENSION(n,n) :: gridj
  INTEGER, DIMENSION(n) :: igrid
  INTEGER, DIMENSION(n,n) :: jgrid
  REAL*8 sumgr
  Do i=1, n
    gridi(i)=i
  END DO
  sumgr=sum(gridi,1)
  write(iout,*) sumgr
  END SUBROUTINE mytest
















