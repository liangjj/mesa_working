c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Recursion Coefficient Interval Conversion}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck cnvtab.f
      subroutine cnvtab(a,b,endpts,n)
      implicit integer(a-z)
      real*8 a, b , endpts, f1, f2
      dimension  a(n), b(n), endpts(2)
      common/io/inp,iout
      f1 = endpts(2)-endpts(1)
      f2 = endpts(1) + endpts(2)
      f2 = - f2/f1
      f1= 2.d0/f1
      do 10 i=1,n
         a(i) = ( a(i) - f2 )/f1
         b(i) = b(i)/f1
 10   continue
      return
      end
