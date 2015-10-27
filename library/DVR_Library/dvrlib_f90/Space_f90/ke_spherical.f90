!deck ke_spherical.f
!***begin prologue     ke_spherical
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements for
!***                   spherical on (-1,1).
!***
!***description
!***references
!***routines called
!***end prologue       ke_spherical
!\begin{eqnarray}
!\end{eqnarray}
  SUBROUTINE ke_spherical(ker,dfar,ptr,wtr,coord,nr,region,nreg,unit_weight)
  USE dvr_global,     ONLY   : parity, mass, iout, m_val
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nr
  REAL*8, DIMENSION(nr,nr)               :: ker, dfar
  REAL*8, DIMENSION(nr)                  :: ptr, wtr
  REAL*8                                 :: ptr_i
  REAL*8                                 :: ptr_j
  CHARACTER (LEN=*)                      :: coord
  INTEGER                                :: region, nreg
  REAL*8                                 :: scale
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  INTEGER                                :: i, j, k
  LOGICAL                                :: unit_weight
!
  Write(iout,1) unit_weight
  ker=0.d0
  IF ( unit_weight) THEN
       DO  i=1,nr
           ptr_i = 1.d0 / ptr(i)
           DO j=1,i
              ptr_j = 1.d0 / ptr(j)
              DO k=1,nr
                 ker(i,j) = ker(i,j) - wtr(k) * dfar(k,i) * dfar(k,j) * ptr(k) * ptr(k)
              END DO
              ker(i,j) = ptr_i * ker(i,j) * ptr_j
              ker(j,i) = ker(i,j)
           END DO
       END DO
  ELSE
!
       DO  i=1,nr
           DO j=1,i
              DO k=1,nr
                 ker(i,j) = ker(i,j) - wtr(k) * dfar(k,i) * dfar(k,j)
              END DO
              ker(j,i) = ker(i,j) 
           END DO
       END DO
  END IF
  scale= -.5D0/mass
  ker=scale*ker
  IF(prn(3)) THEN
     title='unnormalized kinetic energy matrix for '// 'region = '//itoc(region)
     CALL prntrm(title,ker,nr,nr,nr,nr,iout)
  END IF
1 Format(/,5x,'Calling Spherical KE with unit_weight = ',l1)
END SUBROUTINE ke_spherical
