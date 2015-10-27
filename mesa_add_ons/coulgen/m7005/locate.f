*deck @(#)locate.f	1.1 9/8/91
c***begin prologue     locate
c***date written       920924   (yymmdd)
c***revision date               (yymmdd)
c***keywords           spline
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            calculate array of indices of index
c***                   to the left of values of x. used for efficiency
c***                   in spline calculations.
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine locate(xs,ns,x,n,ind)
      implicit integer (a-z)
      dimension xs(ns), x(n), ind(n)
      real *8 xs, x
      logical log1, log2
      common /io/ inp, iout
      do 100 pt=1,n
         jl=0
         ju=ns+1
   10    if (ju-jl.gt.1) then
             jm=(ju+jl)/2
             log1=xs(ns).gt.xs(1)
             log2=x(pt).gt.xs(jm)
             if (log1.eqv.log2) then
                 jl=jm
             else
                 ju=jm
             endif
             go to 10
         endif 
         ind(pt)=jl
  100 continue
      return
      end
