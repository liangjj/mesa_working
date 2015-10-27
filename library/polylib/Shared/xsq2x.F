c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Quadratic Back to Linear Variable}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck xsq2x.f
      subroutine xsq2x(x,end,n)
      implicit integer(a-z)
      real*8 x, end
      dimension  x(n), end(2)
      common/io/inp,iout
      do 10 i=1,n
         x(i) =  sqrt(x(i))
 10   continue   
      end(1)=sqrt(end(1))
      end(2)=sqrt(end(2))
      return
      end
