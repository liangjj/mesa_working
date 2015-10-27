! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Generalized Kinetic Energy Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck kemat.f
!***begin prologue     kemat
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements with
!***                   singularities removed and derivative bloch operator
!***                   added.
!***
!***description
!***references
!***routines called
!***end prologue       kemat
!\begin{eqnarray}
!\end{eqnarray}
  SUBROUTINE kemat(kin,f,df,ddf,pt,wt,n,region)
  USE dvr_global,     ONLY   : parity, mass, iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n,n)                 :: kin, f, df, ddf
  REAL*8, DIMENSION(n)                   :: pt, wt
  INTEGER                                :: region
  REAL*8                                 :: scale
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  INTEGER                                :: i, j, k
!
!  Code using second derivative and Block operator
!
!  kin=0.d0
!  DO  i=1,n
!      kin(i,:) = kin(i,:) + f(i,i)*ddf(i,:)*wt(i)
!  END DO
!
!    Add Bloch contributions
!
!     kin(n,:) = kin(n,:) - f(n,n) * df(n,:)
!     kin(1,:) = kin(1,:) + f(1,1)*df(1,:)
!  scale= -.5D0/mass
!  kin=scale*kin
!  IF(prn(3)) THEN
!     title='unnormalized kinetic energy matrix for '// 'region = '//itoc(region)
!     CALL prntrm(title,kin,n,n,n,n,iout)
!  END IF
!
! Code using first derivatives with no Bloch operator
!
  kin=0.d0
  DO i=1,n
     DO j=1,i
        DO k=1,n 
           kin(i,j) = kin(i,j) - df(k,i) * wt(k) * df(k,j)
        END DO
        kin(j,i) = kin(i,j)
     END DO
  END DO
  scale= -.5D0/mass
  kin=scale*kin
  IF(prn(3)) THEN
     title='unnormalized kinetic energy matrix for '// 'region = '//itoc(region)
     CALL prntrm(title,kin,n,n,n,n,iout)
  END IF
END SUBROUTINE kemat
