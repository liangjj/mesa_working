!deck Test_Cleb
  PROGRAM Test_Cleb
!
  IMPLICIT NONE
  INTEGER                                  :: j1, m1, j2, m2, j3, m3
  REAL*8                                   :: C_3J, F_3J

!
! Read in and compute C3J
!
  DO 
   Read(5,*) j1, m1, j2, m2, j3, m3
   IF (j1<0) THEN
     EXIT
   ELSE
     C_3J = F_3J(j1,m1,j2,m2,j3,m3,.false.)
     write(6,*) j1, m1, j2, m2, j3, m3, C_3J
   END IF
  END DO
  stop
END PROGRAM Test_Cleb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
