*deck @(#)splint.f	1.1 9/8/91
c***begin prologue     splint
c***date written       920924   (yymmdd)
c***revision date               (yymmdd)
c***keywords           spline
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            spline interpolation
c***
c***description        spline interpolate a function and its
c***                   derivative from the spline coefficients
c***                   and index array. 
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine splint(xs,ys,x,y,yp,c,ind,ns,n)
      implicit integer (a-z)
      dimension xs(ns), ys(ns), x(n), y(n), yp(n), c(n), ind(n)
      real *8 xs, ys, x, y, yp, c
      real*8 del, idel, a, b
      real*8 zero, one, three, sixth
      common /io/ inp, iout
      data zero, one, three, sixth / 0.d0, 1.d0, 3.d0,
     1                                .16666666666666666d0 /
      do 100 pt=1,n
         lo=ind(pt)
         hi=lo+1
         del=xs(hi)-xs(lo)
         idel=one/del
         if (del.eq.zero) then
             call lnkerr('error in spline interpolation')
         endif
         a=(xs(hi)-x(pt))*idel
         b=(x(pt)-xs(lo))*idel
         y(pt)=a*ys(lo)+b*ys(hi)+( (a*a*a-a)*c(lo) + (b*b*b-b)*c(hi) )
     1                            *del*del*sixth
         yp(pt)=(ys(hi)-ys(lo))*idel
     1                         - ( three*a*a -one )*c(lo)*del*sixth 
     2                         + ( three*b*b - one )*c(hi)*del*sixth
  100 continue
      return
      end
