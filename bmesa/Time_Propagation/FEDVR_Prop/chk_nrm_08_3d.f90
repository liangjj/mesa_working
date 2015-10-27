!deck chk_nrm_08_3d.f
!***begin prologue     chk_nrm_08_3d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for three dimension wavefunction
!***
!***references
!***routines called    sdot
!***end prologue       chk_nrm_08_3d
  SUBROUTINE chk_nrm_08_3d(v,fac_3d,v_3d)
  USE dvrprop_global
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: v
  REAL*8, DIMENSION(nphy(2),nphy(1))             :: fac_3d
  REAL*8, DIMENSION(nphy(1))                     :: v_3d
  REAL*8                                         :: ddot
  REAL*8                                         :: nrm, first, last
  INTEGER                                        :: i, j, k
!
  nrm = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        fac_3d(j,i) = ddot(nphy(3),v(1,j,i,1),1,v(1,j,i,1),1)   &
                                          +                  &
                   ddot(nphy(3),v(1,j,i,2),1,v(1,j,i,2),1)
     END DO
  END DO
  IF( typke == 'fd' ) then
      DO i=1,nphy(1)
         DO j=1,nphy(2)
            first = v(1,j,i,1) * v(1,j,i,1)                 &
                               +                            &
                    v(1,j,i,2) * v(1,j,i,2)
            last  = v(nphy(3),j,i,1) * v(nphy(3),j,i,1)     &
                                     +                      &
                    v(nphy(3),j,i,2) * v(nphy(3),j,i,2)
            fac_3d(j,i) = fac_3d(j,i) - .5d0 * ( first + last)
         END DO
      END DO
      fac_3d(:,:) = del * fac_3d(:,:)
  END IF
  v_3d(:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        v_3d(i) = v_3d(i) + fac_3d(j,i)
     END DO
  END DO
  IF( typke == 'fd' ) then
      v_3d(:) = v_3d(:) - .5d0 *( fac_3d(1,:) + fac_3d(nphy(2),:) )
      v_3d(:) = del * v_3d(:)
  END IF
  DO i =1,nphy(1)
     nrm = nrm + v_3d(i)
  END DO
  IF( typke == 'fd' ) then
      nrm = nrm -.5d0 * ( v_3d(1) + v_3d(nphy(1)) )
      nrm = del * nrm
  END IF
  write(iout,1) nrm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE chk_nrm_08_3d
