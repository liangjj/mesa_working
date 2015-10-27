*deck plgndr   %W%   %G%
      function plgndr(l,m,x)
c***begin prologue     plgndr
c***date written       910610  (yymmdd)
c***revision date      yymmdd
c
c***keywords           
c***author             press, flannery, teulosky, and vetterling
c***source             %W%   %G%
c***purpose            compute associated legendre polynomial 
c***description
c                                m
c                      computes p (x) by recursion.
c                                l
c                      the function is normalized over (-1.,1.)
c                      the phase is that of abramowitz and stegun,
c                      and gradshteyn and ryzhik.
c                      for either positive or negative m, 
c                      multiplication by exp(i*m*phi)/sqrt(two*pi)
c                      gives the spherical harmonics of condon and shortley,
c                      rose, messiah, and jackson.
c
c***references
c                      w.h.press, b.p.flannery, s.a.teulosky,
c                      and w.t.vetterling, "numerical recipes", 
c                      cambridge university press, p.180,1987.
c
c***routines called
c                      factrl
c
c***end prologue       plgndr
c
c
      implicit integer(a-z)
c
      real*8 plgndr
      real*8 x
c
      real*8 sqrt,abs
      real*8 pmm,sumx2,fact,pmmp1,pll
      real*8 one,two
      real*8 norm,factrl
      parameter (one=1.0d+00,two=2.0d+00)
c
      external factrl
c
c     ----- check for bad arguments -----
      if(     m.gt.l
     $   .or. abs(x).gt.one) then
         call lnkerr('bad arguments to plgndr')
      end if
c
c     generate p(m,m) to initiate the recursion in abs(m).
      mu=abs(m)
      pmm=one
      if(mu.gt.0) then
         sumx2=sqrt((one-x)*(one+x))
         fact=one
         do 10 i=1,mu
            pmm=-pmm*fact*sumx2
            fact=fact+two
   10    continue
      endif
c
c
c     return, or recur on p(m,m) to get p(l,m)
      if(l.eq.mu) then
         plgndr=pmm
      else
c        compute p(m+1,m)
         pmmp1=x*(2*mu+1)*pmm
         if(l.eq.mu+1) then
            plgndr=pmmp1
         else
c           compute p(m+2 ...,m) 
            do 20 ll=mu+2,l
               pll=(x*(2*ll-1)*pmmp1-(ll+mu-1)*pmm)/(ll-mu)
               pmm=pmmp1
               pmmp1=pll
   20       continue
            plgndr=pll
         endif
      endif
c
      norm=sqrt((2*l+1)*factrl(l-mu)/(two*factrl(l+mu)))
      plgndr=norm*plgndr
      if(m.lt.0.and.mod(m,2).ne.0) then
         plgndr=-plgndr
c        phase convention: plgndr=(-one)**m*norm*plgndr
      endif
c
c
      return
      end
