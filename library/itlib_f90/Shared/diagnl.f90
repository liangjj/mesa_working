! \documentclass{article}
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-06  Time: 08:29:08
 
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{DIAGNL: Calculate DVR Diagonals}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck diagnl.f
!***begin prologue     diagnl
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           diagonal
!***author             schneider, barry (nsf)
!***source
!***description        calculate and store the full diagonal
!***                   for the 1,2 or 3 dimensional matrix.
!***                   then zero the diagonal elements of the
!***                   one-dimensional matrices to avoid overcounting.
!***references

!***routines called
!***end prologue       diagnl

SUBROUTINE diagnl(diag,hx,hy,hz,v,n,nx,ny,nz,dim)

REAL*8, INTENT(OUT)                      :: diag(n)
REAL*8, INTENT(IN OUT)                   :: hx(nx,nx)
REAL*8, INTENT(OUT)                      :: hy(ny,ny)
REAL*8, INTENT(OUT)                      :: hz(nz,nz)
REAL*8, INTENT(IN OUT)                   :: v(n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: nx
INTEGER, INTENT(IN)                      :: ny
INTEGER, INTENT(IN)                      :: nz
INTEGER, INTENT(IN)                      :: dim
IMPLICIT INTEGER (a-z)


COMMON/io/inp, iout

CALL copy(v,diag,n)
IF(dim == 1) THEN
  DO  i=1,nx
    diag(i) = diag(i) + hx(i,i)
    hx(i,i) = 0.d0
  END DO
ELSE IF(dim == 2) THEN
  count=0
  DO  i=1,nx
    DO  j=1,ny
      count = count + 1
      diag(count) = diag(count ) + hx(i,i) + hy(j,j)
    END DO
  END DO
  DO  i=1,nx
    hx(i,i) = 0.d0
  END DO
  DO  i=1,ny
    hy(i,i) = 0.d0
  END DO
ELSE IF(dim == 3) THEN
  count=0
  DO  i=1,nx
    DO  j=1,ny
      DO  k=1,nz
        count = count +1
        diag(count) = diag(count) + hx(i,i) + hy(j,j)  &
            + hz(k,k)
      END DO
    END DO
  END DO
  DO  i=1,nx
    hx(i,i) = 0.d0
  END DO
  DO  i=1,ny
    hy(i,i) = 0.d0
  END DO
  DO  i=1,nz
    hz(i,i) = 0.d0
  END DO
ELSE
  CALL lnkerr('error in dimension')
END IF
RETURN
END SUBROUTINE diagnl

