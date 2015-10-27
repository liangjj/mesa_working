! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Band5}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck band5.f
!**begin prologue     band5
!**date written       951229   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            5 point second derivative
!**
!**references
!**routines called
!**end prologue       band5

  SUBROUTINE band5(band,n)
  USE fd_global,   ONLY    : inp, iout, d, del 
  USE fd_prnt
  IMPLICIT NONE
  INTEGER                                  :: n
  REAL*8, DIMENSION(3,n)                   :: band
  CHARACTER (LEN=80)                       :: title
  INTEGER                                  :: i 
  d(3) = -1.d0/(12.d0*del*del)
  d(2) = -16.d0*d(3)
  d(1) = 30.d0*d(3)
! The structure of the input for general band matrices is different
! than for the tridiagonal case.  They are stored as column vectors
! beginning at the diagonal element M(i,i) and proceeding down the column
! in the usual storage pattern.  Since there are only nrow non-zero
! elements for any column, that is all the storage required for the
! matrix.  If eigenvectors are wanted, provide storage in eigvec for them.
! Do all but the last two columns
  do i=1,3
     band(i,:)=d(i)
  end do
  IF(prn_fd_log(4)) THEN
     title='hamiltonian columns'
     CALL prntrm(title,band,3,n,3,n,iout)
  END IF
END SUBROUTINE band5

