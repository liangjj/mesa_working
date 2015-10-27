  PROGRAM test
  IMPLICIT NONE
  INTEGER  :: inp=5, iout=6   
  INTEGER n, i
  INTEGER, DIMENSION(:), POINTER ::ig
INTERFACE
  SUBROUTINE ex2_test(ig,n,inp,iout)   
  INTEGER, DIMENSION(:), POINTER :: ig
  INTEGER n, iout, nloc
  END SUBROUTINE ex2_test
END INTERFACE
  open (inp,file='inp',status='old')
  open (iout,file='out',status='unknown')  
  call ex2_test(ig,n,inp,iout)
  write(iout,*) (ig(i),i=1,20)
  deallocate(ig)
  stop
  END PROGRAM test















