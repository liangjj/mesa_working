! \documentclass{article}
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-03-16  Time: 12:04:26
 
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{TSTOVL: Compute Overlap Matrices}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck tstovl.f

SUBROUTINE tstovl(vec,n,nvc)

REAL*8, INTENT(IN)                       :: vec(n,nvc)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: nvc
IMPLICIT INTEGER (a-z)
REAL*8  ovl, sdot

COMMON/io/inp, iout

DO  i=1,nvc
  DO  j=1,i
    ovl = sdot(n,vec(1,i),1,vec(1,j),1)
    WRITE(iout,1) i, j, ovl
  END DO
END DO
RETURN
1 FORMAT(1X,'i = ',i3,1X,'j = ',i3,1X,'overlap = ',e15.8)
END SUBROUTINE tstovl
