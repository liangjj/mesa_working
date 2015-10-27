!deck auto_08_2d.f
!***begin prologue     auto_08_2d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       auto_08_2d
  SUBROUTINE auto_08_2d(corr_fn,vector_0,vector,fac_2d)
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
END SUBROUTINE auto_08_2d
