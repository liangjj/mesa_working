*deck %W%   %G%
      function ylm(l,m,theta,phi)
c***begin prologue     ylm
c***date written       910610  (yymmdd)
c***revision date      yymmdd
c
c***keywords           
c***author             martin, richard (lanl) 
c***source             %W%   %G%
c***purpose            compute spherical harmonic
c***description
c                                
c                      computes y (theta,phi).
c                                lm
c                      assumes theta and pi arrive in radians
c
c                      the phase convention is described in plgndr
c
c                                      m    *
c                      ylm(l,-m) = (-1) *ylm(l,m)
c
c***references
c                      w.h.press, b.p.flannery, s.a.teulosky,
c                      and w.t.vetterling, "numerical recipes", 
c                      cambridge university press, p.180,1987.
c
c***routines called
c                      plgndr
c
c***end prologue       ylm
c
c
      implicit none
c
c     ----- function returned -----
      complex*16 ylm
c     ----- arguments unchanged -----
      integer l,m
      real*8 theta,phi
c     ----- local variables -----
      real*8 realy,imy
      real*8 zero,one,two,four
      real*8 pi,norm,r,x,mphi
      logical called
c     ----- external functions -----
      real*8 plgndr
      real*8 sqrt
c
      external plgndr
c
      data zero/0.0d+00/, two/2.0d+00/
      data one/1.0d+00/,four/4.0d+00/
      data called/.false./
      save called,norm
c
c     ----- set up constants ----
      if(.not.called) then
         pi=4.0d+00*atan(1.0d+00)
         norm=one/sqrt(two*pi)
         called=.true.
      endif
c
c     ----- map theta to the interval (-1,1) -----
      x=cos(theta)
c
c     check for bad arguments
      if(abs(m).gt.l) then
         call lnkerr('bad arguments to ylm')
      end if
c
c
      r=norm*plgndr(l,m,x)
      mphi=float(m)*phi
      realy=r*cos(mphi)
      imy=r*sin(mphi)
      ylm=dcmplx(realy,imy)
c
c
      return
      end
