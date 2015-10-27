c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Modification of Recursion Coefficients for Lobatto Quadrature}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck modab.f
      subroutine modab(a,b,kpts,endpts,n)
      implicit integer(a-z)
      real*8 a, b , endpts
      real*8 gbslve, gam, t1
      dimension  a(n), b(n), endpts(2)
      common/io/inp,iout
      if (kpts.eq.0) then
          return
      elseif (kpts.eq.1) then

c         only a(n) must be changed

          a(n) =gbslve(endpts(1), n, a, b)*b(n-1)**2 + endpts(1)
          return
      elseif(kpts.eq.2) then

c         a(n) and b(n-1) must be recomputed

          gam =gbslve(endpts(1), n, a, b)
          t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, a, b) 
     1                     - gam))
          b(n-1) =  sqrt(t1)
          a(n) = endpts(1) + gam*t1
          return
      endif
      return
      end
