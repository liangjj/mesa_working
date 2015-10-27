! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Band3}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck band3
!**begin prologue     band3
!**date written       951229   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            3 point second derivative
!**
!**references
!**routines called
!**end prologue       band3

  SUBROUTINE band3(band,n)
  USE fd_global,   ONLY    : inp, iout, d, del
  USE fd_prnt
  IMPLICIT NONE
  INTEGER                                  :: n
  REAL*8, DIMENSION(2,n)                   :: band
  INTEGER                                  :: i
  CHARACTER (LEN=80)                       :: title
  d(1) = -2.d0/(del*del)
  d(2) = 1.d0/(del*del)
!     fill the diagonal and off-diagonal working bands.
  do i=1,2
     band(i,:) = d(i)
  end do
!     print

  IF(prn_fd_log(4)) THEN
     title='hamiltonian bands'
     CALL prntrm(title,band,2,n,2,n,iout)
  END IF
END SUBROUTINE band3




