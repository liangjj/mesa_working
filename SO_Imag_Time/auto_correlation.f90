!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             MODULE auto_correlation
!**begin prologue     auto_correlation
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             INTERFACE auto_correlation_function
             MODULE PROCEDURE auto_corr_1d_d, auto_corr_2d_d,       &
                              auto_corr_3d_d 
                             END INTERFACE auto_correlation_function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!deck auto_corr_1d_d.f
!***begin prologue     auto_corr_1d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       auto_corr_1d_d
  SUBROUTINE auto_corr_1d_d(corr_fn,vector_0,vector)
  USE arnoldi_global
  USE fd_global,         ONLY : del
  USE dvr_shared  
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))             :: vector_0, vector
  REAL*8                                 :: ddot, val
  REAL*8                                 :: a
  COMPLEX*16                             :: corr_fn, dcmplx, conjg
!
  call iosys('read real "initial state" from bec',         &
              nphy(1),vector_0,0,' ')
  a = ddot(nphy(1),vector_0,1,vector,1)       
  corr_fn = a
  val = corr_fn * corr_fn
  write(iout,1) corr_fn, val
1 FORMAT(/,1x,'<psi_0|psi_t>                      = ',e20.12,1x,e20.12, &
         /,1x,'<psi_0|psi_t> Conjg(<psi_0|psi_t>) = ',e20.12 )
END SUBROUTINE auto_corr_1d_d
!deck auto_corr_2d_d.f
!***begin prologue     auto_corr_2d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       auto_corr_2d_d
  SUBROUTINE auto_corr_2d_d(corr_fn,vector_0,vector,g_2)
  USE arnoldi_global
  USE fd_global,         ONLY : del
  USE dvr_shared  
  IMPLICIT NONE
  INTEGER                                  :: i 
  REAL*8, DIMENSION(nphy(2),nphy(1))       :: vector_0, vector
  REAL*8, DIMENSION(nphy(1))               :: g_2
  REAL*8                                   :: ddot, val
  REAL*8                                   :: a
  COMPLEX*16                               :: corr_fn, dcmplx, conjg
!
  call iosys('read real "initial state" from bec',     &
              nphy(2)*nphy(1),vector_0,0,' ')
  DO i=1,nphy(1)
     g_2(i) = ddot(nphy(2),vector_0(1,i),1,vector(1,i),1)
  END DO
  a = 0.d0
  DO i=1,nphy(1)
     a = a + g_2(i)
  END DO
  corr_fn = a
  val = corr_fn * corr_fn
  write(iout,1) corr_fn, val
1 FORMAT(/,1x,'<psi_0|psi_t>                      = ',e20.12,1x,e20.12, /,1x, &
              '<psi_0|psi_t> Conjg(<psi_0|psi_t>) = ',e20.12 )
END SUBROUTINE auto_corr_2d_d
!deck auto_corr_3d_d.f
!***begin prologue     auto_corr_3d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       auto_corr_3d_d
  SUBROUTINE auto_corr_3d_d(corr_fn,vector_0,vector,g_2,g_3)
  USE arnoldi_global
  USE fd_global,         ONLY : del
  USE dvr_shared  
  IMPLICIT NONE
  INTEGER                                           :: i, j 
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))        :: vector_0, vector
  REAL*8, DIMENSION(nphy(1))                        :: g_2
  REAL*8, DIMENSION(nphy(2),nphy(1))                :: g_3
  REAL*8                                            :: ddot, val 
  REAL*8                                            :: a
  COMPLEX*16                                        :: corr_fn, dcmplx
  COMPLEX*16                                        :: conjg
!
  call iosys('read real "initial state" from bec',        &
              nphy(3)*nphy(2)*nphy(1),vector_0,0,' ')
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_3(j,i) = ddot(nphy(3),vector_0(1,j,i),1,vector(1,j,i),1)
     END DO
  END DO
  g_2(:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_2(i) = g_2(i) + g_3(j,i)
     END DO
  END DO
  a = 0.d0
  DO i=1,nphy(1)
     a = a + g_2(i)
  END DO
  corr_fn = a
  val = corr_fn * corr_fn
  write(iout,1) corr_fn, val
1 FORMAT(/,1x,'<psi_0|psi_t>                      = ',e20.12,1x,e20.12, /,1x, &
              '<psi_0|psi_t> Conjg(<psi_0|psi_t>) = ',e20.12 )
END SUBROUTINE auto_corr_3d_d
END MODULE auto_correlation

