  MODULE ex1_mod
  IMPLICIT NONE
  INTEGER  :: inp=5, iout=6   
  INTEGER  :: n
  END MODULE ex1_mod


  PROGRAM ex1
  USE ex1_mod 
  IMPLICIT NONE
  open (inp,file='inp',status='old')
  open (iout,file='out',status='unknown')  
  read(inp,*) n
  call test
  END PROGRAM ex1

  SUBROUTINE TEST
  USE ex1_mod 
  IMPLICIT NONE
  INTEGER   :: i
  INTEGER, DIMENSION(n)  :: itest 
  do i=1,n
   itest(i) = i
  end do
  write(iout,*) itest
  END SUBROUTINE test














