! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Weighted Scalar Product}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck scaprd.f
!***begin prologue     scaprd
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           scalar product for generalized polynomials
!***author             schneider, barry (nsf)
!***source
!***purpose
!***
!***
!***description
!***
!***
!***
!***references

!***routines called
!***end prologue       scaprd
  FUNCTION scaprd (va,vb,wt,rwt,n)
  IMPLICIT NONE
  INTEGER                                :: n, i
  REAL*8, DIMENSION(n)                   :: va, vb, wt, rwt
  REAL*8                                 :: scaprd
!\begin{equation}
!   \langle V_{a} \mid V_{b} \rangle = \sum_{i} V_{a}(i) wt(i) rwt(i) V_{b}(i)
!\end{equation}
! $wt$ is the weight function and $rwt$ the reference weight function used to
! compute the recursion coefficients.
  scaprd=0.d0
  DO  i=1,n
      scaprd=scaprd+va(i)*wt(i)*rwt(i)*vb(i)
  END DO
END FUNCTION scaprd
