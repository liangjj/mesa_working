!***********************************************************************
                           MODULE test_overlap_module
                           USE arnoldi_global
!***********************************************************************
                           INTERFACE test_ovlp
                    MODULE PROCEDURE tst_ovl_d, tst_ovl_z
                       END INTERFACE test_ovlp
!***********************************************************************
                           CONTAINS
!***********************************************************************
!***********************************************************************
!***********************************************************************
!deck tst_ovl_d.f
  SUBROUTINE tst_ovl_d(v_in,nvec)
  IMPLICIT NONE
  INTEGER                                :: n, nvec, i, j 
  REAL*8, DIMENSION(:,:)                 :: v_in
  REAL*8                                 :: ovl, ddot
  DO  i=1,nvec
      DO  j=1,i
          ovl = ddot(n3d,v_in(1,i),1,v_in(1,j),1)
          WRITE(iout,1) i, j, ovl
      END DO
  END DO
1 FORMAT(1X,'i = ',i3,1X,'j = ',i3,1X,'overlap = ',E15.8)
END SUBROUTINE tst_ovl_d
!***********************************************************************
!***********************************************************************
!deck tst_ovl_z.f
  SUBROUTINE tst_ovl_z(v_in,nvec)
  IMPLICIT NONE
  INTEGER                                :: n, nvec, i, j 
  COMPLEX*16, DIMENSION(:,:)             :: v_in
  COMPLEX*16                             :: ovl, cdotc
  DO  i=1,nvec
      DO  j=1,i
          ovl = cdotc(n3d,v_in(1,i),1,v_in(1,j),1)
          WRITE(iout,1) i, j, ovl
      END DO
  END DO
1 FORMAT(1X,'i = ',i3,1X,'j = ',i3,1X,'overlap = ',2E15.8)
END SUBROUTINE tst_ovl_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE test_overlap_module
