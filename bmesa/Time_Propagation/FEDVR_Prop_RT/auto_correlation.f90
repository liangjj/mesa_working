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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE auto_correlation_function
         MODULE PROCEDURE auto_corr_1d_z, auto_corr_2d_z, auto_corr_3d_z
                 END INTERFACE auto_correlation_function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck auto_corr_1d.f
!***begin prologue     auto_corr_1d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       auto_corr_1d_z
  SUBROUTINE auto_corr_1d_z(corr_fn,vector_0,vector)
  USE arnoldi_global
  USE fd_global,         ONLY : del
  USE dvr_shared  
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: vector_0, vector
  REAL*8                                 :: ddot, val, x, xx
  REAL*8, DIMENSION(4)                   :: a
  COMPLEX*16                             :: corr_fn, dcmplx, conjg
  INTEGER                                :: i
!
  call iosys('read real "initial state" from bec',         &
              2*nphy(1),vector_0,0,' ')
  a(1) = ddot(nphy(1),vector_0(1,1),1,vector(1,1),1)       
  a(2) = ddot(nphy(1),vector_0(1,2),1,vector(1,2),1)       
  a(3) = ddot(nphy(1),vector_0(1,2),1,vector(1,1),1)       
  a(4) = ddot(nphy(1),vector_0(1,1),1,vector(1,2),1)       
  a(1) = a(1) + a(2)
  a(2) = a(3) - a(4)
  corr_fn = dcmplx(a(1),a(2))
  val = corr_fn * conjg(corr_fn)
  write(iout,1) corr_fn, val
  x=0.d0
  xx = 0.d0
  DO i=1,nphy(1)
     val = (vector(i,1) * vector(i,1) + vector(i,2) * vector(i,2) ) 
     x = x + val * grid(1)%pt(i)
     xx = xx + val * grid(1)%pt(i) * grid(1)%pt(i)
  END DO
  xx = sqrt(xx - x * x)
  write(iout,2) x, xx 
  write(iplot(6),3) t0, x, xx
1 FORMAT(/,1x,'<psi_0|psi_t>                      = ',e20.12,1x,e20.12,  &
        /,1x, '<psi_0|psi_t> Conjg(<psi_0|psi_t>) = ',e20.12 )
2 FORMAT(/,1x,'<psi_t|x|psi_t>                    = ',e20.12,            &
        /,1x,'sqrt( <psi_t|x*x|psi_t> - <psi_t|x|psi_t> * <psi_t|x|psi_t> ) = ', &
                                                      e20.12)
3 FORMAT(5x,e15.8,5x,e15.8,5x,e15.8)
END SUBROUTINE auto_corr_1d_z
!deck auto_corr_2d_z.f
!***begin prologue     auto_corr_2d_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       auto_corr_2d_z
  SUBROUTINE auto_corr_2d_z(corr_fn,vector_0,vector,g_2)
  USE arnoldi_global
  USE fd_global,         ONLY : del
  USE dvr_shared  
  IMPLICIT NONE
  INTEGER                                :: i, j 
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: vector_0, vector
  REAL*8, DIMENSION(nphy(1),4)           :: g_2
  REAL*8                                 :: ddot, val, x, xx, y, yy
  REAL*8, DIMENSION(4)                   :: a
  COMPLEX*16                             :: corr_fn, dcmplx, conjg
!
  call iosys('read real "initial state" from bec',     &
              2*nphy(2)*nphy(1),vector_0,0,' ')
  g_2(1:nphy(1),:) = 0.d0  
  DO i=1,nphy(1)
     g_2(i,1) = ddot(nphy(2),vector_0(1,i,1),1,vector(1,i,1),1)
     g_2(i,2) = ddot(nphy(2),vector_0(1,i,2),1,vector(1,i,2),1)
     g_2(i,3) = ddot(nphy(2),vector_0(1,i,2),1,vector(1,i,1),1)
     g_2(i,4) = ddot(nphy(2),vector_0(1,i,1),1,vector(1,i,2),1)
  END DO
  a(:) = 0.d0
  DO i=1,nphy(1)
     a(:) = a(:) + g_2(i,:)
  END DO
  a(1) = a(1) + a(2)
  a(2) = a(3) - a(4)
  corr_fn = dcmplx(a(1),a(2))
  val = corr_fn * conjg(corr_fn)
  write(iout,1) corr_fn, val
  x=0.d0
  y=0.d0
  xx = 0.d0
  yy = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        val = ( vector(j,i,1) * vector(j,i,1)               &
                              +                             &
                vector(j,i,2) * vector(j,i,2) ) 
        x = x + val * grid(1)%pt(i)
        y = y + val * grid(2)%pt(j)
        xx = xx + val * grid(1)%pt(i) * grid(1)%pt(i)
        yy = yy + val * grid(2)%pt(j) * grid(2)%pt(j)
     END DO
  END DO
  xx = sqrt(xx - x * x)
  yy = sqrt(yy - y * y)
  write(iout,2) x, xx, y, yy 
  write(iplot(6),3) t0, x, xx, y, yy
1 FORMAT(/,1x,'<psi_0|psi_t>                      = ',e20.12,1x,e20.12, /,1x,      &
              '<psi_0|psi_t> Conjg(<psi_0|psi_t>) = ',e20.12 )
2 FORMAT(/,1x,'<psi_t|x|psi_t>                    = ',e20.12,                      &
         /,1x,'sqrt( <psi_t|x*x|psi_t> - <psi_t|x|psi_t> * <psi_t|x|psi_t> ) = ',  &
                                                      e20.12,                      &
         /,1x,'<psi_t|y|psi_t>                    = ',e20.12,                      &
         /,1x,'sqrt( <psi_t|y*y|psi_t> - <psi_t|y|psi_t> * <psi_t|y|psi_t> ) = ',  &
                                                      e20.12)
3 FORMAT(5x,e15.8,5x,e15.8,5x,e15.8,5x,e15.8,5x,e15.8)
END SUBROUTINE auto_corr_2d_z
!deck auto_corr_3d_z.f
!***begin prologue     auto_corr_3d_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       auto_corr_3d_z
  SUBROUTINE auto_corr_3d_z(corr_fn,vector_0,vector,g_2,g_3)
  USE arnoldi_global
  USE fd_global,         ONLY : del
  USE dvr_shared  
  IMPLICIT NONE
  INTEGER                                        :: i, j, k 
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: vector_0, vector
  REAL*8, DIMENSION(nphy(1),4)                   :: g_2
  REAL*8, DIMENSION(nphy(2),nphy(1),4)           :: g_3
  REAL*8                                         :: ddot, val, x, xx, y, yy, z, zz
  REAL*8, DIMENSION(4)                           :: a
  COMPLEX*16                                     :: corr_fn, dcmplx, conjg
!
  call iosys('read real "initial state" from bec',        &
              2*nphy(3)*nphy(2)*nphy(1),vector_0,0,' ')
  g_3(:,:,:) = 0.d0  
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_3(j,i,1) = ddot(nphy(3),vector_0(1,j,i,1),1,vector(1,j,i,1),1)
        g_3(j,i,2) = ddot(nphy(3),vector_0(1,j,i,2),1,vector(1,j,i,2),1)
        g_3(j,i,3) = ddot(nphy(3),vector_0(1,j,i,2),1,vector(1,j,i,1),1)
        g_3(j,i,4) = ddot(nphy(3),vector_0(1,j,i,1),1,vector(1,j,i,2),1)
     END DO
  END DO
  g_2(:,:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_2(i,:) = g_2(i,:) + g_3(j,i,:)
     END DO
  END DO
  a(:) = 0.d0
  DO i=1,nphy(1)
     a(:) = a(:) + g_2(i,:)
  END DO
  a(1) = a(1) + a(2)
  a(2) = a(3) - a(4)
  corr_fn = dcmplx(a(1),a(2))
  val = corr_fn * conjg(corr_fn)
  write(iout,1) corr_fn, val
  x=0.d0
  y=0.d0
  z=0.d0
  xx = 0.d0
  yy = 0.d0
  zz=0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        DO k=1,nphy(3) 
           val = ( vector(k,j,i,1) * vector(k,j,i,1)               &
                                   +                               &
                   vector(k,j,i,2) * vector(k,j,i,2) )  
           x = x + val * grid(1)%pt(i)
           y = y + val * grid(2)%pt(j)
           z = z + val * grid(3)%pt(k)
           xx = xx + val * grid(1)%pt(i) * grid(1)%pt(i)
           yy = yy + val * grid(2)%pt(j) * grid(2)%pt(j)
           zz = zz + val * grid(3)%pt(k) * grid(3)%pt(k)
        END DO
     END DO
  END DO
  xx = sqrt(xx - x * x)
  yy = sqrt(yy - y * y)
  zz = sqrt(zz - z * z)
  write(iout,2) x, xx, y, yy, z, zz 
  write(iplot(6),3) t0, x, xx, y, yy, z, zz
1 FORMAT(/,1x,'<psi_0|psi_t>                      = ',e20.12,1x,e20.12, /,1x, &
              '<psi_0|psi_t> Conjg(<psi_0|psi_t>) = ',e20.12 )
2 FORMAT(/,1x,'<psi_t|x|psi_t>                    = ',e20.12,                      &
         /,1x,'sqrt( <psi_t|x*x|psi_t> - <psi_t|x|psi_t> * <psi_t|x|psi_t> ) = ',  &
                                                      e20.12,                      &
         /,1x,'<psi_t|y|psi_t>                    = ',e20.12,                      &
         /,1x,'sqrt( <psi_t|y*y|psi_t> - <psi_t|y|psi_t> * <psi_t|y|psi_t> ) = ',  &
                                                      e20.12,                      &
         /,1x,'<psi_t|z|psi_t>                    = ',e20.12,                      &
         /,1x,'sqrt( <psi_t|z*z|psi_t> - <psi_t|z|psi_t> * <psi_t|z|psi_t> ) = ',  &
                                                      e20.12)
3 FORMAT(5x,e15.8,5x,e15.8,5x,e15.8,5x,e15.8,5x,e15.8,5x,e15.8)
END SUBROUTINE auto_corr_3d_z
END MODULE auto_correlation

