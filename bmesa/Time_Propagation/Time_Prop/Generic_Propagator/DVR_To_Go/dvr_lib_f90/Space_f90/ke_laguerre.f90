! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Kinetic Energy Matrix Elements with Laguerre Weight Function}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck ke_laguerre.f
!***begin prologue     ke_laguerre
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           kinetic energy with Laguerre weight function
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements
!***                   for a laguerre weight function.
!***
!***description
!***references
!***routines called
!***end prologue       ke_laguerre
  SUBROUTINE ke_laguerre(kin,f,df,ddf,wt,n,region)
  USE dvr_global,     ONLY   : mass, iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: n    
  REAL*8, DIMENSION(n,n)                 :: kin, f, df, ddf
  REAL*8, DIMENSION(n)                   :: wt
  REAL*8                                 :: scale
  REAL*8                                 :: frac=.25d0
  INTEGER                                :: i, region 
  CHARACTER (LEN=3)                      :: itoc
  CHARACTER (LEN=80)                     :: title
  DO  i=1,n
      kin(i,:) = kin(i,:) + f(i,i) * wt(i) * ( ddf(i,:) - df(i,:) )
  END DO
  DO i=1,n
     kin(i,i) = kin(i,i) + f(i,i) * wt(i) * frac * f(i,i)
  END DO
  scale= -.5D0/mass
  kin=scale*kin
  IF(prn(3)) THEN
     title='unnormalized kinetic energy matrix for ' &
            // 'region = '//itoc(region)
     CALL prntrm(title,kin,n,n,n,n,iout)
  END IF
END SUBROUTINE ke_laguerre

