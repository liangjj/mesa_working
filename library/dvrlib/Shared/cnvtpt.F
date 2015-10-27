c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Point and Weight Interval Conversion}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck cnvtpt.f
      subroutine cnvtpt(pt,wt,endpts,n)
      implicit integer(a-z)
      real*8 pt, wt , endpts, f1, f2
      dimension  pt(n), wt(n), endpts(2)
      common/io/inp,iout
      f1 = ( endpts(2)-endpts(1) )*.5d0
      f2 = ( endpts(1) + endpts(2) )*.5d0
      do 10 i=1,n
         pt(i) =  f1*pt(i) + f2
         wt(i) = wt(i)*f1
 10   continue
      return
      end
