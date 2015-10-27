! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Points}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck points.f
!**begin prologue     points
!**date written       951229   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            points
!**references
!**routines called
!**end prologue       points
  SUBROUTINE points(q_0,q,wt,n)
  USE fd_global,    ONLY    : input, output, del, edge
  USE fd_prnt
  IMPLICIT NONE
  INTEGER                                :: n  
  REAL*8                                 :: q_0
  REAL*8, DIMENSION(n)                   :: q, wt
  REAL*8, DIMENSION(10)                  :: tmp
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: i  
  q_0=edge(1)
  q(1)= q_0 + del
  DO  i=2,n
      q(i) = q(i-1) + del
  END DO

!     Trapezoidal rule

  CALL int_2(tmp,del)
  wt=2.d0*tmp(1)
  IF(prn_fd_log(1)) THEN
     title= 'points'
     CALL prntrm(title,q,n,1,n,1,output)
  END IF
  if(prn_fd_log(2)) then
     title='weights'
     call prntrm(title,wt,n,1,n,1,output)
  END IF
END SUBROUTINE points

