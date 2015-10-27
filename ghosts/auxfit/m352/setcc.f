      subroutine setcc
c***begin prologue     setcc.f
c***date written       770101? (yymmdd)  
c***revision date      920427  (yymmdd)
c   27 april, 1992     rlm at lanl
c      
c***keywords           
c***author             l.e. mcmurchie and e.r. davidson (uw) 
c***source             @(#)setcc.f	2.1   4/29/92
c***purpose            
c***description
c     
c          set up a table of f(t*) expansion coefficients for a
c          7 term taylor series in t over the range 0<t<12. an
c          interval of 0.1 is used in t.
c    
c          the arrays are set up to handle m=0,16
c
c***references
c      l.e. mcmurchie and e.r. davidson, j.comp.phys., 26,218(1978).
c***routines called
c
c***end prologue       setcc.f
      implicit integer(a-z)
c
      real*8 cc,fq,fqi
      real*8 t,zm,fact,ym,yy,w
      real*8 one,two,f43,f45,f47,f2115
      real*8 pt1,pt05
      real*8 error
c
      parameter (one=1.0d0,two=2.0d0,f43=43.0d0,f45=45.0d0,f47=47.0d0)
      parameter (f2115=2115.0d0)
      parameter (pt1=0.1d0,pt05=0.05d0)
      parameter (error=1.d-13)
c
      common/over/cc(2760),fq(0:16),fqi(0:16)
c
      do 20 ii=1,120
         t=float(ii)*pt1 -pt05
         zm=exp(-t)
         fact=f47
         ym=two*t
         yy=ym/f2115
         w=one/f45
11       w=w+yy
         if(yy/w.ge.error) then
            fact=fact+two
            yy=yy*ym/fact
            go to 11
12       endif
         w=w*zm
         m=(ii-1)*23+1
         cc(m)=w
         fact=f43
         do 13 k=1,22
            cc(m+k)=(ym*cc(m+k-1)+zm)/fact
            fact=fact-two
13       continue
20    continue
c
c     ----- compute arrays fq and fqi used in formg -----
c           note that fq(m+1)=2m+1;fqi(m+1) is 1/(2m+1)
c
      fq(0)=one
      fqi(0)=one
      do 30 i=1,16
         fq(i)=fq(i-1)+two
         fqi(i)=one/fq(i)
  30  continue
c
c
      return
      end
