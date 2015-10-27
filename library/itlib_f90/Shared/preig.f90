! \documentclass{article}
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-06  Time: 08:30:01
 
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{PREIG: Print Converged Eigenvalues}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck preig.f
!***begin prologue     preig
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           residual calculation
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       preig

SUBROUTINE preig(eig,num,begin,ipass)

REAL*8, INTENT(IN OUT)                   :: eig(*)
INTEGER, INTENT(IN)                      :: num
INTEGER, INTENT(IN)                      :: begin
INTEGER, INTENT(IN OUT)                  :: ipass
IMPLICIT INTEGER (a-z)


COMMON/io/inp, iout
WRITE(iout,1) ipass

actual=begin
DO  i=1,num
  actual=actual+1
  WRITE(iout,2) actual, eig(i)
END DO
RETURN
1    FORMAT(/,30X,'summary for converged roots on pass = ',i3,  &
    /,17X,'root',22X,'energy')
2    FORMAT(15X,i4,20X,f15.8)
END SUBROUTINE preig

