*deck @(#)vlgndr.f	5.1   11/6/94
      subroutine vlgndr(l,m,plm,x,n,scr)
c***begin prologue     vlgndr.f
c***date written       940304  (yymmdd)
c***revision date      11/6/94
c
c***keywords           legendre polynomials, spherical harmonics
c***author             martin, richard(lanl)
c***source             @(#)vlgndr.f	5.1   11/6/94
c***purpose            compute associated legendre polynomial 
c                      vectorized over the grid points.
c***description
c                                m
c                      computes p (x(1...n)) by recursion.
c                                l
c
c                      given a positive value for m, it returns all
c                      l values from m to l.
c                      the function is normalized over (-1.,1.)
c
c                      needs a scratch array of length n.
c
c                      --- phase conventions ---
c                      for negative m, an appropriate phase is
c                         p(l,m)=(-1**m)*p(l,abs(m))
c                      in which the latter is returned by this routine.
c                      this phase is that of abramowitz and stegun,
c                      and gradshteyn and ryzhik.
c                      with this convention for the legendre function,
c                      then multiplication by exp(i*m*phi)/sqrt(two*pi)
c                      gives the spherical harmonics of condon and shortley,
c                      rose, messiah, and jackson for either positive
c                      or negative m.
c
c
c***references
c                      w.h.press, b.p.flannery, s.a.teulosky,
c                      and w.t.vetterling, "numerical recipes", 
c                      cambridge university press, p.180,1987.
c
c***routines called
c                      factrl
c
c***end prologue       vlgndr.f
      implicit none
c     --- input variables -----
      integer l,m,n
c     --- input arrays (unmodified) ---
      real*8 x(n)
c     --- input arrays (scratch) ---
      real*8 scr(n)
c     --- output arrays ---
      real*8 plm(n,0:l)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
c
      integer i,ll
      integer inp,iout
      real*8 sqrt,abs
      real*8 one,two
      real*8 fact,norm,factrl
      parameter (one=1.0d+00,two=2.0d+00)
c
      common/io/inp,iout
c
      external factrl
c
c     --- check for bad arguments ---
      if(m.lt.0.or.l.lt.0.or.m.gt.l) then
         write(iout,*) 'vlgndr: l,m',l,m
         call lnkerr('bad argument to vlgndr')
      end if
      do 10 i=1,n
         if(abs(x(i)).gt.one) then
            write(iout,*) 'x must be bounded by -1<x<1',i,x(i)
            call lnkerr('bad argument to vlgndr')
         endif
   10 continue
c
c     --- generate p(m,m) to initiate the recursion in abs(m).
c         p(m,m)=(-1**m)*(2m-1!!)*(1-x*x)**m/2
      call vfill(plm(1,m),one,n)
      if(m.gt.0) then
         call vmul(scr,x,x,n)
         call sscal(n,-one,scr,1)
         call sadd(scr,scr,one,n)
         call vsqrt(scr,scr,n)
         fact=one
         do 15 i=1,m
            call vmul(plm(1,m),plm(1,m),scr,n)
c           call sscal(n,-fact,plm(1,m),1)
            call sscal(n,fact,plm(1,m),1)
            fact=fact+two
   15    continue
      endif
c 
c     --- return, or recur upwards on p(m,m) to get p(l,m)
c              p(m+1,m)=x(2m+1)p(m,m)
c         (l-m)p(l,m)  =x(2l-1)p(l-1,m) -(l+m-1)p(l-2,m)
      if(l.eq.m) then
      else
c        --- compute p(m+1,m)
         call vmul(plm(1,m+1),x,plm(1,m),n)
         fact=float(2*m+1)
         call sscal(n,fact,plm(1,m+1),1)
         if(l.eq.m+1) then
         else
c           compute p(m+2 ...,m) 
            do 20 ll=m+2,l
c              fact=float(2*ll-1)/float(ll-m)
               fact=float(2*ll-1)
               call vmul(plm(1,ll),x,plm(1,ll-1),n)
               call sscal(n,fact,plm(1,ll),1)
c              fact=-float(ll+m-1)/float(ll-m)
               fact=-float(ll+m-1)
               call saxpy(n,fact,plm(1,ll-2),1,plm(1,ll),1)
               fact=one/float(ll-m)
               call sscal(n,fact,plm(1,ll),1)
   20       continue
         endif
      endif
c
c     --- apply the normalization constant
c         norm = sqrt [(2l+1)(l-m)! / 2(l+m)!]
      do 30 ll=m,l
         norm=sqrt((2*ll+1)*factrl(ll-m)/(two*factrl(ll+m)))
         call sscal(n,norm,plm(1,ll),1)
   30 continue
c
c
      return
      end
