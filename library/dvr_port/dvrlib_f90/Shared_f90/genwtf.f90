! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Generalized Weight Functions for Orthogonal Functions}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck genwtf.f
!***begin prologue     genwtf
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenwtfctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate weight functions and derivatives on grid.
!***
!***description
!***references
!***routines called
!***end prologue       genwtf
  SUBROUTINE genwtf(pt,arg,wt,dwt,ddwt,alpha,beta,edge,npt,wtfn,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: npt 
  REAL*8, DIMENSION(npt)                 :: pt, arg, wt, dwt, ddwt
  REAL*8                                 :: alpha, beta
  REAL*8, DIMENSION(2)                   :: edge
  CHARACTER (LEN=*)                      :: wtfn
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
!     Get the weight functions and their derivatives
  WRITE(output,1) wtfn
  CALL genrwt(wt,dwt,ddwt,arg,wtfn,alpha,beta,.true., edge,npt)
  IF(prn) THEN
     title='weight function'
     CALL prntfm(title,wt,npt,1,npt,1,output)
     title='first derivative of weight function'
     CALL prntfm(title,dwt,npt,1,npt,1,output)
     title='second derivative of weight function'
     CALL prntfm(title,ddwt,npt,1,npt,1,output)
  END IF
1    FORMAT(/,1X,'type weight function = ',a32)
END SUBROUTINE genwtf
