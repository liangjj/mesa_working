! \documentclass{article}
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-06  Time: 08:30:19
 
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{RDCIV: Read in Converged Davidson Vectors}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck rdciv.f
!***begin prologue     rdciv
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           read ci vectors
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       rdciv

SUBROUTINE rdciv(vec,code,n,nroot)

REAL*8, INTENT(IN OUT)                   :: vec(n,*)
CHARACTER (LEN=*), INTENT(IN OUT)        :: code
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: nroot
IMPLICIT INTEGER (a-z)

CHARACTER (LEN=4) :: itoc


COMMON/io/inp, iout

DO  i=1,nroot
  CALL iosys('read real "'//code//itoc(i)//'" from rwf', n,vec(1,i),0,' ')
!         write(iout,1) i, (vec(j,i),j=1,n)
END DO
RETURN
1    FORMAT(/,1X,'ci vector = ',i3,/,5E15.8)
END SUBROUTINE rdciv

