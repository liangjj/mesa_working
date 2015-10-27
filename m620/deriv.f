*deck @(#)deriv.f	5.1  11/6/94
      function deriv(zl,site,vo)
      implicit none
c
      real*8 deriv
      real*8 zl(3),site(3),vo(4)
c
      real*8 x,y,z,dsqr,d1,d3,one
c
      parameter (one=1.0d+00)
c
      x=vo(1)-site(1)
      y=vo(2)-site(2)
      z=vo(3)-site(3)
      dsqr=x**2+y**2+z**2
c     --- reciprocals of d
      d1=one/sqrt(dsqr)
      d3=d1/dsqr
c     --- get contribution of site
      deriv=(zl(1)*x+zl(2)*y+zl(3)*z)*d3
c
c
      return
      end
