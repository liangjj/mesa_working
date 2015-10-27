!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     MODULE auto_correlation
                     USE arnoldi_global
                     USE dvr_global
                     USE dvr_shared  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
         MODULE PROCEDURE auto_corr_d,                                &
                          auto_corr_z
                 END INTERFACE auto_correlation_function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck auto_corr_d
!***begin prologue     auto_corr_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       auto_corr_d
  SUBROUTINE auto_corr_d(corr_fn,vector_0,vector)
  IMPLICIT NONE
  REAL*8                                 :: corr_fn
  REAL*8, DIMENSION(n3d)                 :: vector_0, vector
!
  IF(spdim == 1 ) THEN
     CALL auto_corr_1d_d(corr_fn,vector_0,vector)
  ELSE IF(spdim == 2 ) THEN
     CALL auto_corr_2d_d(corr_fn,vector_0,vector)
  ELSE IF(spdim == 3 ) THEN
     CALL auto_corr_3d_d(corr_fn,vector_0,vector)
  END IF
END SUBROUTINE auto_corr_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck auto_corr_z
!***begin prologue     auto_corr_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       auto_corr_z
  SUBROUTINE auto_corr_z(corr_fn,vector_0,vector)
  IMPLICIT NONE
  COMPLEX*16                       :: corr_fn
  COMPLEX*16, DIMENSION(n3d)       :: vector_0, vector
!
  IF(spdim == 1 ) THEN
     CALL auto_corr_1d_z(corr_fn,vector_0,vector)
  ELSE IF(spdim == 2 ) THEN
     CALL auto_corr_2d_z(corr_fn,vector_0,vector)
  ELSE IF(spdim == 3 ) THEN
     CALL auto_corr_3d_z(corr_fn,vector_0,vector)
  END IF
END SUBROUTINE auto_corr_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck auto_corr_1d_d
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
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))             :: vector_0, vector
  REAL*8                                 :: corr_fn, val, x, xx
  REAL*8                                 :: ddot
  INTEGER                                :: i
!
  call iosys('read real "initial state" from bec',         &
              nphy(1),vector_0,0,' ')
  corr_fn = ddot(nphy(1),vector_0,1,vector,1)
  val = corr_fn * corr_fn
  write(iout,1) corr_fn, val
  x=0.d0
  xx = 0.d0
  DO i=1,nphy(1)
     val = vector(i) * vector(i) 
     x = x + val * grid(1)%pt(i)
     xx = xx + val * grid(1)%pt(i) * grid(1)%pt(i)
  END DO
  xx = sqrt(xx - x * x)
  write(iout,2) x, xx 
  write(iplot(6),3) t0, x, xx
1 FORMAT(/,1x,'<psi_0|psi_t>                      = ',e20.12,            &
        /,1x, '<psi_0|psi_t> * (<psi_0|psi_t>)    = ',e20.12 )
2 FORMAT(/,1x,'<psi_t|x|psi_t>                    = ',e20.12,            &
        /,1x,'sqrt( <psi_t|x*x|psi_t> - <psi_t|x|psi_t> * <psi_t|x|psi_t> ) = ', &
                                                      e20.12)
3 FORMAT(5x,e15.8,5x,e15.8,5x,e15.8)
END SUBROUTINE auto_corr_1d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck auto_corr_1d_z
!***begin prologue     auto_corr_1d_z
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
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(1))         :: vector_0, vector
  REAL*8                                 :: val, x, xx
  COMPLEX*16                             :: cdotc, corr_fn, conjg
  INTEGER                                :: i
!
  call iosys('read real "initial state" from bec',         &
              2*nphy(1),vector_0,0,' ')
  corr_fn = cdotc(nphy(1),vector_0,1,vector,1)
  val = corr_fn * conjg(corr_fn)
  write(iout,1) corr_fn, val
  x=0.d0
  xx = 0.d0
  DO i=1,nphy(1)
     val = vector(i) * conjg(vector(i)) 
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  SUBROUTINE auto_corr_2d_d(corr_fn,vector_0,vector)
  IMPLICIT NONE
  INTEGER                                :: i, j 
  REAL*8, DIMENSION(nphy(2),nphy(1))     :: vector_0, vector
  REAL*8                                 :: corr_fn
  REAL*8                                 :: ddot, val, x, xx, y, yy
!
  call iosys('read real "initial state" from bec',     &
              nphy(2)*nphy(1),vector_0,0,' ')
  corr_fn = ddot(nphy(2)*nphy(1),vector_0,1,vector,1)
  val = corr_fn * corr_fn
  write(iout,1) corr_fn, val
  x=0.d0
  y=0.d0
  xx = 0.d0
  yy = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        val = vector(j,i) * vector(j,i)  
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
1 FORMAT(/,1x,'<psi_0|psi_t>                      = ',e20.12,/,1x,                 &
              '<psi_0|psi_t> * (<psi_0|psi_t>)    = ',e20.12 )
2 FORMAT(/,1x,'<psi_t|x|psi_t>                    = ',e20.12,                      &
         /,1x,'sqrt( <psi_t|x*x|psi_t> - <psi_t|x|psi_t> * <psi_t|x|psi_t> ) = ',  &
                                                      e20.12,                      &
         /,1x,'<psi_t|y|psi_t>                    = ',e20.12,                      &
         /,1x,'sqrt( <psi_t|y*y|psi_t> - <psi_t|y|psi_t> * <psi_t|y|psi_t> ) = ',  &
                                                      e20.12)
3 FORMAT(5x,e15.8,5x,e15.8,5x,e15.8,5x,e15.8,5x,e15.8)
END SUBROUTINE auto_corr_2d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  SUBROUTINE auto_corr_2d_z(corr_fn,vector_0,vector)
  IMPLICIT NONE
  INTEGER                                :: i, j 
  COMPLEX*16, DIMENSION(nphy(2),nphy(1)) :: vector_0, vector
  REAL*8                                 :: val, x, xx, y, yy
  COMPLEX*16                             :: cdotc, corr_fn, conjg
!
  call iosys('read real "initial state" from bec',     &
              2*nphy(2)*nphy(1),vector_0,0,' ')
  corr_fn = cdotc(nphy(2)*nphy(1),vector_0,1,vector,1)
  val = corr_fn * conjg(corr_fn)
  write(iout,1) corr_fn, val
  x=0.d0
  y=0.d0
  xx = 0.d0
  yy = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        val = vector(j,i) * conjg ( vector(j,i) ) 
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  SUBROUTINE auto_corr_3d_d(corr_fn,vector_0,vector)
  IMPLICIT NONE
  INTEGER                                        :: i, j, k 
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))     :: vector_0, vector
  REAL*8                                         :: corr_fn
  REAL*8                                         :: ddot, val, x, xx, y, yy, z, zz
!
  call iosys('read real "initial state" from bec',        &
              nphy(3)*nphy(2)*nphy(1),vector_0,0,' ')
  corr_fn = ddot(nphy(3)*nphy(2)*nphy(1),vector_0,1,vector,1)
  val = corr_fn * corr_fn
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
           val = vector(k,j,i) * vector(k,j,i)                     
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
1 FORMAT(/,1x,'<psi_0|psi_t>                      = ',e20.12,/,1x,                 &
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
END SUBROUTINE auto_corr_3d_d
!***********************************************************************
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
  SUBROUTINE auto_corr_3d_z(corr_fn,vector_0,vector)
  IMPLICIT NONE
  INTEGER                                        :: i, j, k 
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1)) :: vector_0, vector
  REAL*8                                         :: val, x, xx, y, yy, z, zz
  COMPLEX*16                                     :: cdotc, corr_fn, conjg
!
  call iosys('read real "initial state" from bec',        &
              2*nphy(3)*nphy(2)*nphy(1),vector_0,0,' ')
  corr_fn = cdotc(nphy(3)*nphy(2)*nphy(1),vector_0,1,vector,1)
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
           val = vector(k,j,i) * conjg( vector(k,j,i) )  
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
!***********************************************************************
!***********************************************************************
END MODULE auto_correlation

