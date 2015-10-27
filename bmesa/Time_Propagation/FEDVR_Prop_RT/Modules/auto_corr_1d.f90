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
!***end prologue       auto_corr_1d
  SUBROUTINE auto_corr_1d(corr_fn,vector_0,vector)
  USE arnoldi_global
  USE fd_global,         ONLY : del
  USE dvr_shared  
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: vector_0, vector
  REAL*8                                 :: ddot, val, first, last
  REAL*8, DIMENSION(4)                   :: a
  COMPLEX*16                             :: corr_fn, dcmplx, conjg
!
  call iosys('read real "initial state" from bec',         &
              2*nphy(1),vector_0,0,' ')
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
END SUBROUTINE auto_corr_1d
