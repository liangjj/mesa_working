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
  SUBROUTINE ke_cyl(kin,f,df,pt,wt,n,region)
  USE dvr_global,     ONLY   : parity, mass, iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n,n)                 :: kin, f, df
  REAL*8, DIMENSION(n)                   :: pt, wt
  INTEGER                                :: region
  REAL*8                                 :: scale
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  INTEGER                                :: i, j, k
!
  kin=0.d0
  DO i=1,n
     DO j=1,i
        DO k=1,n
           kin(i,j) = kin(i,j) - wt(k) * df(k,i) *  df(k,j)
        END DO
     END DO
  END DO
  scale= -.5D0/mass
  kin=scale*kin
  IF(prn(3)) THEN
     title='unnormalized kinetic energy matrix for '// 'region = '//itoc(region)
     CALL prntrm(title,kin,n,n,n,n,iout)
  END IF
END SUBROUTINE ke_cyl
