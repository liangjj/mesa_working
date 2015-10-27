!***********************************************************************
                           MODULE test_overlap
                           INTERFACE tst_ovl
                    MODULE PROCEDURE tst_ovl_d, tst_ovl_z
                       END INTERFACE tst_ovl
!***********************************************************************
                           CONTAINS
!***********************************************************************
!***********************************************************************
!***********************************************************************
!deck tst_ovl_z.f
  SUBROUTINE tst_ovl_z(v_in,nvec)
  USE arnoldi_global_rt
  IMPLICIT NONE
  INTEGER                                :: n, nvec, i, j 
  COMPLEX*16, DIMENSION(n3d,nvec)        :: v_in
  COMPLEX*16                             :: ovl, cdotc
  DO  i=1,nvec
      DO  j=1,i
          ovl = cdotc(n3d,v_in(1,i),1,v_in(1,j),1)
          WRITE(iout,1) i, j, ovl
      END DO
  END DO
1 FORMAT(1X,'i = ',i3,1X,'j = ',i3,1X,'overlap = ',2E15.8)
END SUBROUTINE tst_ovl_z
!
!deck tst_ovl_d.f
  SUBROUTINE tst_ovl_d(v_in,nvec)
  USE arnoldi_global_it
  IMPLICIT NONE
  INTEGER                                :: n, nvec, i, j 
  REAL*8, DIMENSION(n3d,nvec)            :: v_in
  REAL*8                                 :: ovl, ddot
  DO  i=1,nvec
      DO  j=1,i
          ovl = ddot(n3d,v_in(1,i),1,v_in(1,j),1)
          WRITE(iout,1) i, j, ovl
      END DO
  END DO
1 FORMAT(1X,'i = ',i3,1X,'j = ',i3,1X,'overlap = ',E15.8)
END SUBROUTINE tst_ovl_d
END MODULE test_overlap 
