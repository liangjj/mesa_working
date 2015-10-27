c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Points, Weights and Traansformation Matrix for 
c        Generalized Gauss Quadratures}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle    

*deck @(#)genq.f	1.1 9/9/91
      subroutine genq( a, b, w, z, fix, edge, muzero, n)
c
c        a, b     on input the recursion coefficients; will be destroyed
c
c     output parameters (both arrays of length n)
c
c        a        the desired points
c        w        the desired weights
c
      implicit integer (a-z)
      real*8 a, b, w, z, muzero, sumwt, edge
      logical fix
      dimension  a(n), b(n), z(n,n), w(n), edge(2), fix(2) 
      common/io/inp,iout

c     the method used is a ql-type method with origin shifting

      call copy(b,z,n)
      call copy(z,b(2),n-1)
      call rzero(w,n)
      call rzero(z,n*n)
      do 10 i=1,n
         z(i,i)=1.d0
 10   continue   
      call imtql2 (n, n, a, b, z, ierr)
      sumwt=0.d0
      do 20 i = 1, n
         w(i) = muzero * z(1,i) * z(1,i)
         sumwt=sumwt+w(i)
 20   continue   
      if(fix(1)) then
         a(1)=edge(1)
      endif
      if(fix(2)) then
         a(n)=edge(2)
      endif
      write(iout,1) sumwt
      return
 1    format(/,1x,'sum of the weights = ',e15.8)
      end
