*deck @(#)sphere.f	5.3 4/17/95
      subroutine sphere(mxpts,p,w,l,npts)
c***begin prologue     sphere.f
c***date written       012095
c***revision date      4/17/95
c***keywords           spherical harmonics, numerical integration,
c                      angular quadrature 
c***author             r.l. martin, lanl
c***source             @(#)sphere.f	5.3 4/17/95
c***purpose            returns the points and weights associated with
c                      gaussian quadrature on the unit sphere.
c***description
c                      this routine returns the points and weights associated
c                      with integration on the unit sphere. in particular,
c                      quadratures of a given order l exactly integrate
c                      all surface harmonics of degree l or less.
c                      the jacobian [sin(theta)dtheta dphi] is implicitly 
c                      included, so the working formula for integrating a
c                      polynomial using this routine is simply Sum f(p(i))*w(i).
c
c                  mxpts  ---  leading dimension of p.
c                      p  ---  quadrature points (x/r,y/r,z/r), output.
c                      w  ---  quadrature weights, output.
c                      l  ---  order requested, input.
c                   npts  ---  number of points in the quadrature, input.
c                              the correspondence between order and npts is:
c                              l     npts
c                              ----------
c                              2        4    tetrahedral
c                              3        6    octahedral
c                              5       12    icosahedral formula
c                              5       16    abramowitz and stegun
c                              7       26    abramowitz and stegun
c                              9       38    lebedev
c                             11       50    lebedev
c                             13       74    lebedev
c                             15       86    lebedev
c                             17      110    lebedev
c                             19      146    lebedev
c                             23      194    lebedev
c                             29      302    lebedev
c                             35      434    lebedev/treutler,alrichs
c                             41      590    lebedev
c
c    
c
c***references
c                      v.i. lebedev, zh. vychisl. mat. mat. fiz, 15,48(1975)
c                                    ibid.,16,293(1976) -- note the erratum.
c                                    english translations in
c                                    U.S.S.R. Comput. Math. and Math. Phys.
c              
c                      see also
c                      a.d. mclaren, math. comput. 17, 361(1963).
c                      a.h. stroud, approximate calculation of multiple
c                                   integrals (prentice-hall,1971)
c***routines called
c                      none
c***end prologue       sphere.f
      implicit none
c
c     --- input variables ---
      integer mxpts,l,npts
c     --- input arrays (unmodified) ---
c     --- input arrays (modified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 p(mxpts,3),w(mxpts)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,j
      real*8 a1
      real*8 zero,half,one,two,three,four,five,six,eight,twelve,fifteen
      real*8 pi,atan,ctheta,stheta
      data zero/0.0d+00/, half/0.5d+00/, one/1.0d+00/ 
      data two/2.0d+00/, three/3.0d+00/, four/4.0d+00/ 
      data five/5.0d+00/
      data six/6.0d+00/, eight/8.0d+00/, twelve/12.0d+00/
      data fifteen/15.0d+00/
c
      save zero,half,one,two,three,four,six,eight,twelve,fifteen,pi
c
      integer inp,iout
      common/io/inp,iout
c
c     --- set the quadrature weights depending on order ---
      pi=four*atan(one)
      if(l.eq.2) then
c        tetrahedral formula
         a1=sqrt(one/three)
         p(1,1)=a1
         p(1,2)=a1
         p(1,3)=a1
c
         p(2,1)=-a1
         p(2,2)=-a1
         p(2,3)= a1
c
         p(3,1)=-a1
         p(3,2)= a1
         p(3,3)=-a1
c
         p(4,1)= a1
         p(4,2)=-a1
         p(4,3)=-a1
c
         do 10 i=1,4
            w(i)=four*pi/four
   10    continue
c        
      else if (l.eq.3) then
         call lebedev(mxpts,p,w,l,npts)
      else if(l.eq.5) then
         if(npts.eq.12) then
c           icoshaedral formula.
            ctheta=sqrt(three/fifteen)
            stheta=sqrt(twelve/fifteen)
            p(1,1)=zero
            p(1,2)=zero
            p(1,3)=one
c
            p(2,1)=stheta
            p(2,2)=zero
            p(2,3)=ctheta
c
            p(3,1)=cos(two*pi/five)*stheta
            p(3,2)=sin(two*pi/five)*stheta
            p(3,3)=ctheta
c
            p(4,1)=cos(four*pi/five)*stheta
            p(4,2)=sin(four*pi/five)*stheta
            p(4,3)=ctheta
c
            p(5,1)=cos(six*pi/five)*stheta
            p(5,2)=sin(six*pi/five)*stheta
            p(5,3)=ctheta
c
            p(6,1)=cos(eight*pi/five)*stheta
            p(6,2)=sin(eight*pi/five)*stheta
            p(6,3)=ctheta
c
c           generate remainder by inversion.
            do 20 i=1,6
               do 15 j=1,3
                  p(i+6,j)=-p(i,j)
   15          continue
               w(i)=(one/twelve)*four*pi
               w(i+6)=(one/twelve)*four*pi
   20       continue
         else
c           use the 18 point formula invariant under octahedral symmetry.
            call lebedev(mxpts,p,w,l,npts)
         endif
      else if(l.eq.7) then
         call lebedev(mxpts,p,w,l,npts)
      else if(l.eq.9) then
         call lebedev(mxpts,p,w,l,npts)
      else if(l.eq.11) then
         call lebedev(mxpts,p,w,l,npts)
      else if(l.eq.13) then
         call lebedev(mxpts,p,w,l,npts)
      else if(l.eq.15) then
         call lebedev(mxpts,p,w,l,npts)
      else if(l.eq.17) then
         call lebedev(mxpts,p,w,l,npts)
      else if(l.eq.19) then
         call lebedev(mxpts,p,w,l,npts)
      else if(l.eq.23) then
         call lebedev(mxpts,p,w,l,npts)
      else if(l.eq.29) then
         call lebedev(mxpts,p,w,l,npts)
      else if(l.eq.35) then
         call lebedev(mxpts,p,w,l,npts)
      else if(l.eq.41) then
         call lebedev(mxpts,p,w,l,npts)
      else
         call plnkerr(' order requested not yet available',3100)
      endif
c
c
      return
      end
