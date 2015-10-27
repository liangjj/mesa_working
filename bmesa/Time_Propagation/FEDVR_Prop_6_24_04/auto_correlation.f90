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
  IMPLICIT NONE
  CONTAINS
!deck auto_correlation_function
!***begin prolgue      auto_correlation_function
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate autocorrelation function.
!***
!***references
!***routines called
!***end prologue       auto_correlation_function
!
  SUBROUTINE auto_correlation_function
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
!
!
  call iosys('read real "initial state" from bec',2*n3d,v_scr,0,' ')
  IF(spdim==1) THEN
     CALL auto_1d(auto_corr,v_scr,psi)
  ELSE IF(spdim==2) THEN
     CALL auto_2d(auto_corr,v_scr,psi,fac)
  ELSE IF(spdim==3) THEN
     CALL auto_3d(auto_corr,v_scr,psi,fac,v_1)
  END IF
END SUBROUTINE auto_correlation_function
!deck auto_1d.f
!***begin prologue     auto_1d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       auto_1d
  SUBROUTINE auto_1d(corr_fn,vector_0,vector)
  USE arnoldi_global
  USE fd_global,         ONLY : del
  USE dvr_shared  
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: vector_0, vector
  REAL*8                                 :: ddot, val, first, last
  REAL*8, DIMENSION(4)                   :: a
  COMPLEX*16                             :: corr_fn, dcmplx, conjg
!
  a(1) = ddot(nphy(1),vector_0(1,1),1,vector(1,1),1)       
  a(2) = ddot(nphy(1),vector_0(1,2),1,vector(1,2),1)       
  a(3) = ddot(nphy(1),vector_0(1,2),1,vector(1,1),1)       
  a(4) = ddot(nphy(1),vector_0(1,1),1,vector(1,2),1)       
  IF (typke == 'fd' ) THEN
      first = vector_0(1,1) * vector(1,1)                  
      last  = vector_0(nphy(1),1) * vector(nphy(1),1)
      a(1) = a(1) - .5d0 * ( first + last )
      first = vector_0(1,2) * vector(1,2)
      last  = vector_0(nphy(1),2) * vector(nphy(1),2)
      a(2) = a(2) - .5d0 * ( first + last )
      first = vector_0(1,2) * vector(1,1)
      last  = vector_0(nphy(1),2) * vector(nphy(1),1)
      a(3) = a(3) - .5d0 * ( first + last )
      first = vector_0(1,2) * vector(1,2)
      last  = vector_0(nphy(1),2) * vector(nphy(1),2)
      a(4) = a(4) - .5d0 * ( first + last )
      a(:) = a(:) * del
  END IF
  a(1) = a(1) + a(2)
  a(2) = a(3) - a(4)
  corr_fn = dcmplx(a(1),a(2))
  val = corr_fn * conjg(corr_fn)
  write(iout,1) corr_fn, val
1 FORMAT(/,1x,'<psi_0|psi_t>                      = ',e20.12,1x,e20.12, /,1x, &
              '<psi_0|psi_t> Conjg(<psi_0|psi_t>) = ',e20.12 )
END SUBROUTINE auto_1d
!deck auto_2d.f
!***begin prologue     auto_2d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       auto_2d
  SUBROUTINE auto_2d(corr_fn,vector_0,vector,fac_2d)
  USE arnoldi_global
  USE fd_global,         ONLY : del
  USE dvr_shared  
  IMPLICIT NONE
  INTEGER                                :: i 
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: vector_0, vector
  REAL*8, DIMENSION(nphy(1),4)           :: fac_2d
  REAL*8                                 :: ddot, val, first, last
  REAL*8, DIMENSION(4)                   :: a
  COMPLEX*16                             :: corr_fn, dcmplx, conjg
!
  fac_2d(:,:) = 0.d0  
  DO i=1,nphy(1)
     fac_2d(i,1) = ddot(nphy(2),vector_0(1,i,1),1,vector(1,i,1),1)
     fac_2d(i,2) = ddot(nphy(2),vector_0(1,i,2),1,vector(1,i,2),1)
     fac_2d(i,3) = ddot(nphy(2),vector_0(1,i,2),1,vector(1,i,1),1)
     fac_2d(i,4) = ddot(nphy(2),vector_0(1,i,1),1,vector(1,i,2),1)
  END DO
  IF (typke == 'fd' ) THEN
      DO i=1,nphy(1)
         first = vector_0(1,i,1) * vector(1,i,1)                  
         last  = vector_0(nphy(2),i,1) * vector(nphy(2),i,1)
         fac_2d(i,1) = fac_2d(i,1) - .5d0 * ( first + last )
         first = vector_0(1,i,2) * vector(1,i,2)
         last  = vector_0(nphy(2),i,2) * vector(nphy(2),i,2)
         fac_2d(i,2) = fac_2d(i,2) - .5d0 * ( first + last )
         first = vector_0(1,i,2) * vector(1,i,1)
         last  = vector_0(nphy(2),i,2) * vector(nphy(2),i,1)
         fac_2d(i,3) = fac_2d(i,3) - .5d0 * ( first + last )
         first = vector_0(1,i,1) * vector(1,i,2)
         last  = vector_0(nphy(2),i,1) * vector(nphy(2),i,2)
         fac_2d(i,4) = fac_2d(i,4) - .5d0 * ( first + last )
      END DO
      fac_2d(:,:) = del * fac_2d(:,:)
  END IF
  a(:) = 0.d0
  DO i=1,nphy(1)
     a(:) = a(:) + fac_2d(i,:)
  END DO
  IF (typke == 'fd' ) THEN
      a(:) = a(:) - .5d0 * ( fac_2d(1,:) + fac_2d(nphy(1),:) ) 
      a(:) = del * a(:)
  END IF
  a(1) = a(1) + a(2)
  a(2) = a(3) - a(4)
  corr_fn = dcmplx(a(1),a(2))
  val = corr_fn * conjg(corr_fn)
  write(iout,1) corr_fn, val
1 FORMAT(/,1x,'<psi_0|psi_t>                      = ',e20.12,1x,e20.12, /,1x, &
              '<psi_0|psi_t> Conjg(<psi_0|psi_t>) = ',e20.12 )
END SUBROUTINE auto_2d
!deck auto_3d.f
!***begin prologue     auto_3d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       auto_3d
  SUBROUTINE auto_3d(corr_fn,vector_0,vector,fac_3d,v_3d)
  USE arnoldi_global
  USE fd_global,         ONLY : del
  USE dvr_shared  
  IMPLICIT NONE
  INTEGER                                        :: i, j 
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: vector_0, vector
  REAL*8, DIMENSION(nphy(2),nphy(1),4)           :: fac_3d
  REAL*8, DIMENSION(nphy(1),4)                   :: v_3d
  REAL*8                                         :: ddot, val, first, last
  REAL*8, DIMENSION(4)                           :: a
  COMPLEX*16                                     :: corr_fn, dcmplx, conjg
!
  fac_3d(:,:,:) = 0.d0  
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        fac_3d(j,i,1) = ddot(nphy(3),vector_0(1,j,i,1),1,vector(1,j,i,1),1)
        fac_3d(j,i,2) = ddot(nphy(3),vector_0(1,j,i,2),1,vector(1,j,i,2),1)
        fac_3d(j,i,3) = ddot(nphy(3),vector_0(1,j,i,2),1,vector(1,j,i,1),1)
        fac_3d(j,i,4) = ddot(nphy(3),vector_0(1,j,i,1),1,vector(1,j,i,2),1)
     END DO
  END DO
  IF (typke == 'fd' ) THEN
      DO i=1,nphy(1)
         DO j=1,nphy(2)
            first = vector_0(1,j,i,1) * vector(1,j,i,1)                  
            last  = vector_0(nphy(3),j,i,1) * vector(nphy(3),j,i,1)
            fac_3d(j,i,1) = fac_3d(j,i,1) - .5d0 * ( first + last )
            first = vector_0(1,j,i,2) * vector(1,j,i,2)
            last  = vector_0(nphy(3),j,i,2) * vector(nphy(3),j,i,2)
            fac_3d(j,i,2) = fac_3d(j,i,2) - .5d0 * ( first + last )
            first = vector_0(1,j,i,2) * vector(1,j,i,1)
            last  = vector_0(nphy(3),j,i,2) * vector(nphy(3),j,i,1)
            fac_3d(j,i,3) = fac_3d(j,i,3) - .5d0 * ( first + last )
            first = vector_0(1,j,i,1) * vector(1,j,i,2)
            last  = vector_0(nphy(3),j,i,1) * vector(nphy(3),j,i,2)
            fac_3d(j,i,4) = fac_3d(j,i,4) - .5d0 * ( first + last )
         END DO
      END DO
      fac_3d(:,:,:) = del * fac_3d(:,:,:)
  END IF
  v_3d(:,:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        v_3d(i,:) = v_3d(i,:) + fac_3d(j,i,:)
     END DO
  END DO
  IF (typke == 'fd' ) THEN
      DO i=1,nphy(1)
         v_3d(i,:) = v_3d(i,:) - .5d0 * ( fac_3d(1,i,:) + fac_3d(nphy(2),i,:) )
      END DO
      v_3d(:,:) = del * v_3d(:,:)
  END IF
  a(:) = 0.d0
  DO i=1,nphy(1)
     a(:) = a(:) + v_3d(i,:)
  END DO
  IF (typke == 'fd' ) THEN
      a(:) = a(:) - .5d0 * ( v_3d(1,:) + v_3d(nphy(1),:) ) 
      a(:) = del * a(:)
  END IF
  a(1) = a(1) + a(2)
  a(2) = a(3) - a(4)
  corr_fn = dcmplx(a(1),a(2))
  val = corr_fn * conjg(corr_fn)
  write(iout,1) corr_fn, val
1 FORMAT(/,1x,'<psi_0|psi_t>                      = ',e20.12,1x,e20.12, /,1x, &
              '<psi_0|psi_t> Conjg(<psi_0|psi_t>) = ',e20.12 )
END SUBROUTINE auto_3d
END MODULE auto_correlation

