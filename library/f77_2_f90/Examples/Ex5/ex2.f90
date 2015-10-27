  PROGRAM test
  IMPLICIT NONE
  INTEGER  :: inp=5, iout=6   
  INTEGER n, i
  INTEGER, DIMENSION(:,:), ALLOCATABLE ::ig
  open (inp,file='inp',status='old')
  open (iout,file='out',status='unknown')  
  read(inp,*) n
  ALLOCATE(ig(n,1))
  ig(:,1)=0
  write(iout,*) ig(:,1)
  ALLOCATE(ig(n,2))
  ig(:,2)=3
  write(iout,*) ig(:,2)
  deallocate(ig)
  stop
  END PROGRAM test















