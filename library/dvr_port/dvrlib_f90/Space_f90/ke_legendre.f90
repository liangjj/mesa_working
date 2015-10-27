!deck ke_legendre.f
!***begin prologue     ke_legendre
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements for
!***                   legendre on (-1,1).
!***
!***description
!***references
!***routines called
!***end prologue       ke_legendre
!\begin{eqnarray}
!\end{eqnarray}
  SUBROUTINE ke_legendre(ker,far,dfbr,ddfbr,ptr,wtr,coord,nr,region)
  USE dvr_global,     ONLY   : parity, mass, output, m_val
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nr
  REAL*8, DIMENSION(nr,nr)               :: ker, far, dfbr, ddfbr
  REAL*8, DIMENSION(nr)                  :: ptr, wtr
  CHARACTER (LEN=*)                      :: coord
  INTEGER                                :: region
  REAL*8                                 :: scale
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  INTEGER                                :: i, j
!
  write(output,*) ptr
  ker=0.d0
  DO  i=1,nr
      ker(i,:) = ker(i,:) + wtr(i) * far(i,i) *                  &
                         ( ( 1.d0-ptr(i)*ptr(i) ) * ddfbr(i,:)   &
                              - 2.d0 * ptr(i) *dfbr(i,:) )
  END DO
  IF(m_val /= 0) then
     DO i=1,nr
        ker(i,i) = ker(i,i) - m_val*m_val/(1.d0 - ptr(i)*ptr(i))
     END DO
  END IF
  scale= -.5D0/mass
  ker=scale*ker
  IF(prn(3)) THEN
     title='unnormalized kinetic energy matrix for '// 'region = '//itoc(region)
     CALL prntrm(title,ker,nr,nr,nr,nr,output)
  END IF
END SUBROUTINE ke_legendre
