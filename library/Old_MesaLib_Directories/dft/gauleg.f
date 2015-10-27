*deck @(#)gauleg.f	5.1 11/6/94
      subroutine gauleg(p,w,lorder,npts)
c***begin prologue     gauleg.f
c***date written       930520  (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           gauss-legendre quadrature, integration, unit sphere
c***author             martin, richard (lanl)
c***source             @(#)gauleg.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       gauleg.f
      implicit none
c     --- input variables ---
      integer lorder,npts
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 p(3,npts),w(npts)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
c
      integer ltop
      integer ntheta,nphi
      integer theta,phi,i
      integer inp,iout
      parameter (ltop=101)
      real*8 ut(ltop),up(ltop)
      real*8 wt(ltop),wphi
      real*8 zero,half,one,two,four,pi
      logical debug
      parameter (zero=0.0d+00,half=0.5d+00,one=1.0d+00,two=2.0d+00)
      parameter (four=4.0d+00)
      parameter (debug=.true.)
c
      common/io/inp,iout
c
      external lgndrx
c
c     -----------------------------------------------------------------
      pi=4.0d+00*atan(1.0d+00)
c
c     --- in order to exactly integrate a spherical harmonic of
c         degree l, we need (l+1)/2 points in the theta quadrature,
c         and l+1 points in the phi quadrature.
c
c         also note that ntheta points in the gauss-legendre will 
c         exactly integrate any polynomial of degree 2*ntheta-1
c         in cos(theta).
      if(lorder+1.gt.ltop) then
         call lnkerr(' harmonic order too large for gauleg')
      endif
      if(mod(lorder+1,2).eq.0) then
         ntheta=(lorder+1)/2
      else
         ntheta=(lorder+1)/2 +1
      endif
c      nphi=lorder+1
      nphi=2*ntheta
c
c     --- were we given enough room ---
      if(ntheta*nphi.gt.npts) then
         call lnkerr(' need more room in gauleg')
      endif
c
c     --- calculate points and weights of gauss-legendre 
c           quadrature in theta
      do 5 theta=1,ntheta
         call lgndrx(ntheta,theta,wt(theta),ut(theta))
    5 continue
c
c     --- convert from the (0,1) interval of lgndrx to (-1,1). ---
      do 10 theta=1,ntheta
         wt(theta)=two*wt(theta)
         ut(theta)=two*ut(theta)-one
   10 continue
      if(debug) then
         write(iout,*)' gauss-legendre'
         write(iout,*) ' weights w=', (wt(i),i=1,ntheta)
         write(iout,*) ' points u=', (ut(i),i=1,ntheta)
      endif
c
c     --- calculate points and weights of simple quadrature on phi.
      wphi=two*pi/float(nphi)
      do 60 phi=1,nphi
         up(phi)=two*pi*(float(phi-1))/float(nphi)
   60 continue
      if(debug) then
         write(iout,*)' gauss weight wphi=', wphi
         write(iout,*) ' points u=', (up(phi),phi=1,nphi)
      endif
c
c     --- put them together ---
      i=0
      do 70 theta=1,ntheta
         do 65 phi=1,nphi
            i=i+1
            p(1,i)=cos(up(phi))*sin(acos(ut(theta)))
            p(2,i)=sin(up(phi))*sin(acos(ut(theta)))
            p(3,i)=ut(theta)
            w(i)=wt(theta)*wphi
   65    continue
   70 continue
c
      npts=ntheta*nphi
c
c
      return
      end
