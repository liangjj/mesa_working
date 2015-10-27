! \documentclass{article}
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-03-16  Time: 12:04:05
 
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{RDIAG: Driver for Diagonalization of Small Davidson Matrix}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck rdiag.f
!***begin prologue     rdiag
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           diagonalization
!***author             schneider, barry (nsf)
!***source
!***purpose            driver for real diagonalization.
!***
!***references

!***routines called
!**end prologue       rdiag

SUBROUTINE rdiag(mat,eig,vec,work,rep,iter,n,m,prnt)

REAL*8, INTENT(IN)                       :: mat(n,*)
REAL*8, INTENT(IN)                       :: eig(*)
REAL*8, INTENT(OUT)                      :: vec(n,*)
REAL*8, INTENT(OUT)                      :: work(*)
REAL*8, INTENT(IN)                       :: rep
INTEGER, INTENT(IN OUT)                  :: iter
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: m
LOGICAL, INTENT(IN)                      :: prnt
IMPLICIT INTEGER (a-z)
CHARACTER (LEN=80) :: title
CHARACTER (LEN=3) :: itoc



COMMON/io/inp, iout

CALL dsyev('v','l',m,mat,n,eig,work,5*m,info)
IF(info /= 0) THEN
  CALL lnkerr('error from direct diagonalization routine')
END IF
DO  i=1,m
  DO  j=1,m
    vec(j,i)=mat(j,i)
  END DO
END DO
IF(prnt) THEN
  DO  i=1,m
    work(i)=eig(i)+rep
  END DO
  title='eigenvalues of small matrix iteration = ' //itoc(iter)
  CALL prntfm(title,work,m,1,m,1,iout)
END IF
RETURN
END SUBROUTINE rdiag

