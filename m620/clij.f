*deck @(#)clij.f	5.1  11/6/94
      function clij(site,vo,l)
      implicit none
c
      real*8 clij
      integer l
      real*8 site(3),vo(4)
c
      real*8 x,y,z,dsqr,d1,d3,d5
      real*8 zero,one,half,three
c
      parameter (zero=0.0d+00,half=0.5d+00,one=1.0d+00,three=3.0d+00)
c
      x=vo(1)-site(1)
      y=vo(2)-site(2)
      z=vo(3)-site(3)
      dsqr=x**2+y**2+z**2
c     --- reciprocals of d
      d1=one/sqrt(dsqr)
      d3=d1/dsqr
      d5=d3/dsqr
c
      if(l.eq.1) then
         clij=d1
      else if(l.eq.2) then
         clij=x*d3
      else if(l.eq.3) then
         clij=y*d3
      else if(l.eq.4) then
         clij=z*d3
      else if(l.eq.5) then
         call lnkerr('m620: problem in clij. never vary qxx')
      else if(l.eq.6) then
         clij=half*(y*y-x*x)*d5
      else if(l.eq.7) then
         clij=half*(z*z-x*x)*d5
      else if(l.eq.8) then
         clij=three*x*y*d5
      else if(l.eq.9) then
         clij=three*x*z*d5
      else if(l.eq.10) then
         clij=three*y*z*d5
      endif
c
c
      return
      end
