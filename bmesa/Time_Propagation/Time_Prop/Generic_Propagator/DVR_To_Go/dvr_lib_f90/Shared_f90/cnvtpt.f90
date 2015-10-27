! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Point and Weight Interval Conversion}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck cnvtpt.f
  SUBROUTINE cnvtpt(pt,wt,endpts,n)
  USE input_output
  IMPLICIT NONE
  INTEGER                             :: n
  REAL*8, DIMENSION(n)                :: pt, wt
  REAL*8, DIMENSION(2)                :: endpts
  REAL*8                              :: f1, f2
  f1 = ( endpts(2)-endpts(1) )*.5D0
  f2 = ( endpts(1) + endpts(2) )*.5D0
  pt =  f1*pt + f2
  wt = wt*f1
END SUBROUTINE cnvtpt
