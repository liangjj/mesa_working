!deck tstovl.f
  SUBROUTINE tstovl
  USE io
  USE dvd_global
  USE dvd_prnt
  IMPLICIT NONE
  INTEGER                            :: i, j
  REAL*8                             :: ovl, sdot
DO  i=1,size
    DO  j=1,i
        ovl = sdot(n_dvd,vec(1,i),1,vec(1,j),1)
        WRITE(iout,1) i, j, ovl
    END DO
  END DO
1 FORMAT(1X,'i = ',i3,1X,'j = ',i3,1X,'overlap = ',e15.8)
END SUBROUTINE tstovl
