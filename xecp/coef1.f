*deck %W% %G%
      subroutine coef1(ltop,ntheta,nphi,acoef1,len1,ptcf1,wt,ut,up,pl,
     $                 csn,sn,fre,fim)
      implicit integer(a-z)
c
c
c     calculate projection of cartesian powers on spherical harmonics.
c     r t pack, june 1992.
c
c
      complex*16 acoef1(len1)
      integer ptcf1(0:ltop,0:ltop,0:ltop)
      real*8 wt(ntheta),ut(ntheta),up(nphi)
      real*8 pl(ntheta,0:ltop,0:ltop)
      real*8 csn(nphi,0:ltop),sn(nphi,0:ltop)
      real*8 fre(ntheta,nphi),fim(ntheta,nphi)
c
c
      real*8 pi,wphi,x,y,ar,ai
      real*8 sin,cos,asin,acos,atan
      real*8 plgndr,gauleg
      real*8 zero,half,one,two,four
      logical debug
      parameter (zero=0.0d+00,half=0.5d+00,one=1.0d+00,two=2.0d+00)
      parameter (four=4.0d+00)
      parameter (debug=.false.)
c
      common/io/inp,iout
c
      external plgndr
      external gauleg
c
c     -----------------------------------------------------------------
      pi=four*atan(one)
c
c     ntheta points in the gauss-legendre will exactly integrate a
c     polynomial of degree 2*ntheta-1 in cos(theta). we expect polynomials
c     of order 2*ijmax.
c     the lower bound 4 arises from the lower limit
c     imposed in the routine lgndrx  which returns the points 
c     and weights.
      nt=max(ltop+1,4)
      if(nt.gt.ntheta) then
         call lnkerr(' ntheta too small for exact quadrature')
      endif
c
c     nphi points in the gauss-chebychev will take care of polynimials
c     of degree 2*nphi+1 in either sin(phi) or cos(phi).
      if(ltop+1.gt.nphi) then
         call lnkerr(' nphi too small for exact quadrature')
      endif
c
c     ----- calculate points and weights of gauss-legendre 
c           quadrature in theta
      do 5 theta=1,ntheta
         call lgndrx(ntheta,theta,wt(theta),ut(theta))
    5 continue
c     ----- convert from the (0,1) interval of lgndrx to (-1,1). -----
      do 10 theta=1,ntheta
         wt(theta)=two*wt(theta)
         ut(theta)=two*ut(theta)-one
   10 continue
c     ----- construct matrix of associated legendre polynomials.
      do 50 l=0,ltop
         do 40 m=0,l
            do 30 theta=1,ntheta
                  pl(theta,l,m)=plgndr(l,m,ut(theta))
                  if(debug) then
                     write(iout,1000) l,m,ut(theta),pl(theta,l,m)
1000                 format(' plm=',2i5,2d13.6)
                  endif
   30       continue
   40    continue
   50 continue
c     ----- calculate points and weight of chebyshev quadrature on phi.
c           note that we include the normalization factor in the weight.
      wphi=two*pi/float(nphi)
      wphi=wphi/(sqrt(two*pi))
      do 60 phi=1,nphi
         up(phi)=cos((float(phi)-half)*wphi)
   60 continue
      if(debug) then
         write(iout,*)' gauss-chebychev weight wphi=', wphi
         write(iout,*) ' points u=', (up(phi),phi=1,nphi)
      endif
c
c     ----- construct matrix of cos(m*phi) and sin(m*phi) -----
c     note that the real part of the integral is a polynomial
c     in cos(phi), while the imaginary part is polynomial in
c     sin(phi).  we do these as two separate quadratures, and
c     can use the same points and weights, but in one case the point
c     x is interpreted as cos(phi), while x=sin(phi) in the other.
c***** why is this wrong??
c     note the sign in csn which comes about from the dx=-sin(phi)d(phi),
c     whereas in x=sin(phi), dx=cos(phi)d(phi)
      do 80 phi=1,nphi
         do 70 m=0,ltop
            x=acos(up(phi))
            csn(phi,m)=cos(m*x)
            y=asin(up(phi))
            sn(phi,m)=sin(m*y)
            if(debug) then
               write(iout,1020) m,phi,csn(phi,m)
               write(iout,1030)m,phi,sn(phi,m)
 1020          format(' m,phi,cos(m*phi)=',2i5,d13.6)
 1030          format(' m,phi,sin(m*phi)=',2i5,d13.6)
            endif
   70    continue
   80 continue
c
c     ----- we're finally ready to roll. 
      int=0
      if(debug) then 
         write(iout,1040) 
 1040    format(1x,'  i  j  k  l  m')
      endif
c
      do 220 i=0,ltop
         do 210 j=0,ltop
            do 200 k=0,ltop
               if((i+j+k).le.ltop) then
c
c                 ----- evaluate the cartesian function x**i y**j z**k
c                    with radial dependence removed at each of the
c                    quadrature points.
c                    recall that ut is cos(theta) and up is interpreted
c                    as either cos(phi) or sin(phi) depending on whether
c                    we're doing the real or imaginary component.
                  do 100 theta=1,ntheta
                     do 90 phi=1,nphi
                        fre(theta,phi)=ut(theta)**k
     $                                *up(phi)**i
     $                                *(sqrt(one-up(phi)*up(phi)))**j
     $                       *(sqrt(one-ut(theta)*ut(theta)))**(i+j)
                        fim(theta,phi)=ut(theta)**k
     $                                *(sqrt(one-up(phi)*up(phi)))**i
     $                                *up(phi)**j
     $                       *(sqrt(one-ut(theta)*ut(theta)))**(i+j)
   90                continue
  100             continue
c
c                 ----- project out spherical harmonic components
c                       of function. the parity of l must be the same
c                       as i+j+k.
                  lbeg=mod(i+j+k,2)
                  ptcf1(i,j,k)=int
                  do 150 l=lbeg,ltop,2
                     do 140 m=0,l
c                        ----- the quadrature over cos(m*phi); 
c                              the real part
****must change dimensioning here to use this
                        ar=gauleg(ntheta,nphi,pl(1,l,m),wt,
     $                            fre,csn(1,m),wphi,pi)
                        ai=gauleg(ntheta,nphi,pl(1,l,m),wt,
     $                            fim,sn(1,m),wphi,pi)
c                       ----- handle some special cases in which
c                             the function cannot be written as a
c                             simple polynomial and so the quadrature
c                             will not be exact. fortunately, the result
c                             can be determined by symmetry -----
                        if (mod(j,2).ne.0) ar=zero
                        if(mod(i+m,2).eq.0) ai=zero
c                       ----- put this result in the output list -----
                        int=int+1
                        acoef1(int)=dcmplx(ar,ai)
                        if(debug) then
                           write(iout,*) i,j,k,l,m,ar,ai
                        endif
  140                continue
  150             continue
               endif
  200       continue
  210    continue
  220 continue
c
c
      return
      end
