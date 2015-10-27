! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Filv}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck filv.f
!***begin prologue     filv
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            fill up the global potential from the 
!***                   regional arrays.
!***                   
!***references

!***routines called
!***end prologue       filv

   SUBROUTINE filv(v,vi,ni,nfun,nglobal,start)
!
  USE input_output
!
  IMPLICIT NONE
  INTEGER                                   :: ni, nfun, nglobal, start
  REAL*8, DIMENSION(nglobal)                :: v
  REAL*8, DIMENSION(ni)                     :: vi
  INTEGER                                   :: cnti, i
!
  CHARACTER (LEN=80) :: title
  cnti=start
  DO  i=1,nfun-1
      v(i)=vi(cnti)
      cnti=cnti+1
  END DO
  v(nfun) = vi(cnti)
END SUBROUTINE filv
