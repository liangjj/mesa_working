  SUBROUTINE init_vec
  USE lanczos_global
  IMPLICIT NONE 
  INTEGER                     :: i
  do i=1,n
     v(i,0) = 1.d0
  end do
END SUBROUTINE init_vec
