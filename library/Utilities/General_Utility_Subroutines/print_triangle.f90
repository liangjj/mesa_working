!deck print_triangle
  subroutine print_triangle(title,a,n,iout)
  IMPLICIT NONE
  CHARACTER*80                :: title
  REAL*8, DIMENSION(:)        :: a
  INTEGER                     :: n
  INTEGER                     :: i
  INTEGER                     :: j
  INTEGER                     :: count
  INTEGER                     :: iout
  write(iout,1) title
  count = 0
  DO i=1,n
     write(iout,2) i
     write(iout,3) (a(j), j=count+1, count + i )
     count = count + i
  END DO
1 Format(a80)
2 Format(/,5x,'Row = ',i4)
3 Format( (15x,5f10.5) )
END SUBROUTINE print_triangle
