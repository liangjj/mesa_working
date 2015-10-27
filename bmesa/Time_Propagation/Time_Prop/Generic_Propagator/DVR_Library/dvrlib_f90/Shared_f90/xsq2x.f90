! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Quadratic Back to Linear Variable}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck xsq2x.f
  SUBROUTINE xsq2x(x,END,n)
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: x
  REAL*8, DIMENSION(2)                   :: END
  x =  SQRT(x)
  END(1)=SQRT(END(1))
  END(2)=SQRT(END(2))
END SUBROUTINE xsq2x
