! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Generalized First Derivative Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck p_hermite.f
!***begin prologue     p_hermite
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate generalized first derivative matrix elements.
!***references
!***routines called
!***end prologue       p_hermite
!\begin{eqnarray}
!\end{eqnarray}
  SUBROUTINE p_hermite(p_momr,fa,fb,dfb,pt,wt,n,region)
  USE dvr_global,     ONLY   : parity, mass, output
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n,n)                 :: p_momr, fa, fb, dfb
  REAL*8, DIMENSION(n)                   :: pt, wt
  INTEGER                                :: region
  REAL*8                                 :: scale
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  INTEGER                                :: i
!
  p_momr=0.d0
  DO  i=1,n
      p_momr(i,:) = p_momr(i,:) + fa(i,i) * wt(i) * dfb(i,:)
  END DO
  DO i=1,n
     p_momr(i,i) = p_momr(i,i) -  fa(i,i) * wt(i) * pt(i) * fb(i,i)
  END DO
  scale= .5D0
  p_momr = scale*p_momr
  IF(prn(3)) THEN
     title='unnormalized first derivative matrix for '// 'region = '//itoc(region)
     CALL prntrm(title,p_momr,n,n,n,n,output)
  END IF
END SUBROUTINE p_hermite
