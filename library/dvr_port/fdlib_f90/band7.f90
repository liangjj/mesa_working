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
!deck band7.f
!**begin prologue     band7
!**date written       951229   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            7 point second derivative
!**
!**references
!**routines called
!**end prologue       band7
  SUBROUTINE band7(band,n)
  USE fd_global,   ONLY    : input, output, d, del
  USE fd_prnt
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(4,n)                 :: band
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: i
  d(4) = 2.d0/(720.d0*del*del)
  d(3) = -54.d0*d(4)
  d(2) = 540.d0*d(4)
  d(1) = -980.d0*d(4)
  d(4) = 4.d0*d(4)
! The structure of the input for general band matrices is different
! than for the tridiagonal case.  They are stored as column vectors
! beginning at the diagonal element M(i,i) and proceeding down the column
! in the usual storage pattern.  Since there are only nrow non-zero
! elements for any column, that is all the storage required for the
! matrix.  If eigenvectors are wanted, provide storage in eigvec for them.
! Do all but the last three columns
  DO  i=1,4
      band(i,:) = d(i)
  END DO
  IF(prn_fd_log(4)) THEN
     title='hamiltonian columns'
     CALL prntrm(title,band,4,n,4,n,output)
  END IF
END SUBROUTINE band7

