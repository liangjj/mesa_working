c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Points and Weights for A Quadratic VAriable}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck toquad.f
      subroutine toquad(pt,endpts,n)
      implicit integer(a-z)
      real*8 pt, endpts
      dimension  pt(n), endpts(2)
      common/io/inp,iout
      do 10 i=1,n
         pt(i) =  sqrt(pt(i))
 10   continue   
      endpts(1)=sqrt(endpts(1))
      endpts(2)=sqrt(endpts(2))
      return
      end
