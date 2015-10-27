*deck @(#)fndbrk.f	1.1 9/8/91
c***begin prologue     fndbrk
c***date written       910128   (yymmdd)
c***revision date               (yymmdd)
c***keywords           spline
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            calculate array of indices of largest breakpoint
c***                   to the left of values of x. used for efficiency
c***                   in bspline calculations.
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine fndbrk(x,break,ind,n,l)
      implicit integer (a-z)
      dimension x(n), break(l), ind(n)
      real *8 x, break
      common /io/ inp, iout
      do 10 i=1,n
         call interv(break,l,x(i),ind(i),ndum)
   10 continue     
      return
      end
