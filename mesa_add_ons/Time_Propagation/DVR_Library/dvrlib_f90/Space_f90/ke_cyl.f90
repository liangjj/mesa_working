!deck ke_cyl.f
!***begin prologue     ke_cyl
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements for
!***                   cylindrical case.
!***
!***description
!***references
!***routines called
!***end prologue       ke_cyl
!\begin{eqnarray}
!\end{eqnarray}
  SUBROUTINE ke_cyl(ker,far,dfar,ddfar,ptr,wtr,coord,nr,region,nreg)
  USE dvr_global,     ONLY   : parity, mass, iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nr
  REAL*8, DIMENSION(nr,nr)               :: ker, far, dfar, ddfar
  REAL*8, DIMENSION(nr)                  :: ptr, wtr
  CHARACTER (LEN=*)                      :: coord
  INTEGER                                :: region, nreg
  REAL*8                                 :: scale
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  INTEGER                                :: i, j, k
!
  ker=0.d0
  DO i=1,nr
     DO j=1,i
        DO k=1,nr
           ker(i,j) = ker(i,j) - wtr(k) * dfar(k,i) *  dfar(k,j)
        END DO
        ker(j,i)=ker(i,j)
     END DO
  END DO
!  ker=0.d0
!  DO i=1,nr
!     ker(i,:) = ker(i,:) + wtr(i) * far(i,i) * ( ddfar(i,:) + dfar(i,:)/ptr(i) )
!  END DO
!  ker(nr,:) = ker(nr,:) - ptr(nr) * far(nr,nr) * dfar(nr,:)
!  ker(1,:)  = ker(1,:)  + ptr(1)  * far(1,1)   * dfar(1,:)
  scale= -.5D0/mass
  ker=scale*ker
  IF(prn(3)) THEN
     title='unnormalized kinetic energy matrix for '// 'region = '//itoc(region)
     CALL prntrm(title,ker,nr,nr,nr,nr,iout)
  END IF
END SUBROUTINE ke_cyl
