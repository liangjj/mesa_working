! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{ Coordinate Functions and Derivatives}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck cpoly.f
!***begin prologue     cpoly
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate coordinate functions.
!***
!***description
!***references
!***routines called
!***end prologue       cpoly

  SUBROUTINE cpoly(cp,dcp,ddcp,pt,arg,n,npt,parity,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                             :: n, npt
  REAL*8, DIMENSION(npt,0:n)          :: cp, dcp, ddcp
  REAL*8, DIMENSION(npt)              :: pt, arg
  CHARACTER (LEN=*)                   :: parity
  LOGICAL                             :: prn
  REAL*8                              :: fac, fac2, tmp
  CHARACTER (LEN=80)                  :: title
  IF(parity /= 'none') THEN
     arg = pt*pt
  ELSE
     arg = pt
  END IF
  CALL lgngr(cp,dcp,ddcp,arg,arg,npt,npt,.false.,'all')
  IF(prn) THEN
     title='coordinate function'
     CALL prntfm(title,cp(1,0),npt,n+1,npt,n+1,output)
     title='first derivative of coordinate function'
     CALL prntfm(title,dcp(1,0),npt,n+1,npt,n+1,output)
     title='second derivative of coordinate function'
     CALL prntfm(title,ddcp(1,0),npt,n+1,npt,n+1,output)
   END IF
END SUBROUTINE cpoly
