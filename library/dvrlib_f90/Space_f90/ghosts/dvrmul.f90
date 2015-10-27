! \documentclass{article}
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:41:59
 
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{DVR Matrix Vector Mutiply}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck dvrmul.f
!***begin prologue     dvrmul
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       dvrmul

SUBROUTINE dvrmul(vecin,vecout,diag,hbuf,ihbuf,n,nvc,nonz)


REAL*8, INTENT(IN)                       :: vecin(n,nvc)
REAL*8, INTENT(OUT)                      :: vecout(n,nvc)
REAL*8, INTENT(IN)                       :: diag(n)
REAL*8, INTENT(IN)                       :: hbuf(*)
INTEGER, INTENT(IN)                      :: ihbuf(2,*)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: nvc
INTEGER, INTENT(IN)                      :: nonz
IMPLICIT INTEGER (a-z)
REAL*8  hij


COMMON/io/inp, iout

CALL rzero(vecout,n*nvc)

!     do the diagonal

DO  i=1,nvc
  DO  j=1,n
    vecout(j,i) = vecout(j,i) + diag(j)*vecin(j,i)
  END DO
END DO

!     Now the off-diagonal in packed form.

DO  i=1,nonz
  ii=ihbuf(1,i)
  jj=ihbuf(2,i)
  hij=hbuf(i)
  DO  j=1,nvc
    vecout(ii,j) = vecout(ii,j) + hij*vecin(jj,j)
    vecout(jj,j) = vecout(jj,j) + hij*vecin(ii,j)
  END DO
END DO
RETURN
END SUBROUTINE dvrmul
















