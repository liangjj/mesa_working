! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{LANCZOS_MULT: Multiplication of Hamiltonian on a Vector}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck lanczos_mult
!**begin prologue     lanczos_mult
!**date written       010828   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           matrix vector
!**author             schneider, barry (nsf)
!**source             tproplib
!**purpose            matrix-vector, multiply, dvr.
!**description        multiply a hamiltonian matrix onto a vector.
!**end prologue       lanczos_mult
  SUBROUTINE lanczos_mult(v_in,v_out)
  USE lanczos_global
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n)             :: v_in, v_out
  call cebc(v_out,matrix,v_in,n,n,1)
  title='hamiltonian on vectors'
  CALL prntcm(title,v_out,n,1,n,1,iout)
END SUBROUTINE lanczos_mult



