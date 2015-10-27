! \documentclass{article}
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-03-16  Time: 12:03:05
 
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{INVEC: Read in CI Vectors for Davidson Code}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck invec.f
!***begin prologue     invec
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           read in vectors
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       invec

SUBROUTINE invec(vec,TYPE,n,nvc,prnt)

REAL*8, INTENT(IN OUT)                   :: vec(n,nvc)
CHARACTER (LEN=*), INTENT(IN OUT)        :: TYPE
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: nvc
INTEGER, INTENT(IN)                      :: prnt
IMPLICIT INTEGER (a-z)

CHARACTER (LEN=4) :: itoc
CHARACTER (LEN=80) :: title


COMMON/io/inp, iout

DO  i=1,nvc
  CALL iosys('read real "'//TYPE//itoc(i)//'" from rwf', n,vec(1,i),0,' ')
END DO
IF(prnt) THEN
  title='input vectors'
  CALL prntrm(title,vec,n,nvc,n,nvc,iout)
END IF
RETURN
END SUBROUTINE invec

