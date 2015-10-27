!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE auto_correlation
!**begin prologue     auto_correlation module
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       auto_correlation
!
IMPLICIT NONE
  INTERFACE auto_correlation_function
    SUBROUTINE auto_corr_1d(corr_fn,vector_0,vector)
      USE arnoldi_global
      USE fd_global,         ONLY : del
      USE dvr_shared  
      IMPLICIT NONE
      REAL*8, DIMENSION(nphy(1))                     :: vector_0, vector
      REAL*8                                         :: ddot, val, first, last
      REAL*8                                         :: a
      COMPLEX*16                                     :: corr_fn, dcmplx, conjg
    END SUBROUTINE auto_corr_1d
    SUBROUTINE auto_corr_2d(corr_fn,vector_0,vector,g_2)
      USE arnoldi_global
      USE fd_global,         ONLY : del
      USE dvr_shared  
      IMPLICIT NONE
      INTEGER                                        :: i 
      REAL*8, DIMENSION(nphy(2),nphy(1))             :: vector_0, vector
      REAL*8, DIMENSION(nphy(1))                     :: g_2
      REAL*8                                         :: ddot, val, first, last
      REAL*8                                         :: a
      COMPLEX*16                                     :: corr_fn, dcmplx, conjg
    END SUBROUTINE auto_corr_2d
    SUBROUTINE auto_corr_3d(corr_fn,vector_0,vector,g_2,g_3)
      USE arnoldi_global
      USE fd_global,         ONLY : del
      USE dvr_shared  
      IMPLICIT NONE
      INTEGER                                         :: i, j 
      REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))      :: vector_0, vector
      REAL*8, DIMENSION(nphy(1))                      :: g_2
      REAL*8, DIMENSION(nphy(2),nphy(1))              :: g_3
      REAL*8                                          :: ddot, val, first, last
      REAL*8                                          :: a
      COMPLEX*16                                      :: corr_fn, dcmplx, conjg
    END SUBROUTINE auto_corr_3d
  END INTERFACE auto_correlation_function
END MODULE auto_correlation

