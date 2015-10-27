  MODULE test_global
  IMPLICIT NONE
  INTEGER  :: inp=5, iout=6
  INTEGER  :: m   
  INTEGER, DIMENSION(100) :: nr, nc
  INTEGER                 :: offset 
TYPE pntr
  REAL*8, DIMENSION(:,:), POINTER     :: g
  REAL*8, DIMENSION(:), POINTER       :: h
END TYPE pntr
TYPE(pntr), DIMENSION(:), ALLOCATABLE ::pg
  END MODULE test_global









