! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Recursion Coefficient Interval Conversion}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck cnvtab.f
  SUBROUTINE cnvtab(a,b,endpts,n)
  USE input_output
  INTEGER                      :: n
  REAL*8, DIMENSION(n)         :: a, b
  REAL*8, DIMENSION(2)         :: endpts
  REAL*8                       :: f1, f2
  f1 = endpts(2)-endpts(1)
  f2 = endpts(1) + endpts(2)
  f2 = - f2/f1
  f1= 2.d0/f1
  a = ( a - f2 )/f1
  b = b/f1
END SUBROUTINE cnvtab
