  Subroutine Make_Array(a,l,u) 
  USE accuracy
  USE input_output
  IMPLICIT NONE
  INTEGER                   :: l
  INTEGER                   :: u
  INTEGER                   :: i
  REAL(idp), DIMENSION(l:u) :: a
  DO i = l, u
     a(i) = i
  END DO
  write(iout,*) a(l:u)
END SUBROUTINE Make_Array
