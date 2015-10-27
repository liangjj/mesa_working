*deck @(#)recur1.f	2.1  10/10/91
      subroutine recur1 (nstart,nmax,lmax,q,alpha,rk,argab)
      implicit real*8(a-h,o-z)
c
c     controls computation of q(n,l) by recursion.
c
c     ----- arguments unchanged -----
      integer nstart,nmax,lmax
      real*8 alpha,rk
c     ----- arguments returned -----
      real*8 q(0:nmax,0:lmax),argab
c     integer and
c     ----- local variables -----
      real*8 zero,one,two,three,four
      real*8 half
      real*8 talph
c     ----- common -----
      real*8 pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00,three=3.0d+00)
      parameter (four=4.0d+00,half=0.5d+00)
c
c********** change this routine to recognize that lambda must have same
c           parity as ltot or else the angular integral will kill it???
c**** no, keep this as a general routine
c
c     ----- initialize the output array q -----
      do 10 n=0,nstart+lmax
         do 10 lambda=0,n
            q(n,lambda)=zero
   10    continue
   20 continue
c
c     ----- check for special case of zero argument -----
c           note the early return
      if (rk.eq.zero) then
         nhi=nstart+lmax
         argab=zero
         q(0,0)=half*sqrt(pi/alpha)
         q(1,0)=half/alpha
c        ----- recursion 38d in mcmurchie and davidson -----
         do 30 n=2,nhi
            q(n,0)=half*float(n-1)*q(n-2,0)/alpha
   30    continue
         return
      endif
c
c     ----- normal processing -----
      t=rk*rk/(four*alpha)
      argab=-t
      if (lmax.eq.0) then
         q(nstart,0)=qcomp(nstart,0,alpha,rk,t)
         return
      else if(lmax.eq.1) then
         q(nstart,0)=qcomp(nstart,0,alpha,rk,t)
         q(nstart+1,1)=qcomp(nstart+1,1,alpha,rk,t)
         return
      else
         talph=alpha+alpha
         if (nstart.eq.0) then
c           ----- n=0 terms in potential -----
c                 recursion 38a in mcmurchie and davidson 
            q(2,0)=qcomp(2,0,alpha,rk,t)
            do 50 l=0,lmax
               q(l+2,l)=rk*q(l+1,l-1)/talph
   50       continue
            nbeg=4
            if (t.le.three) then
c              ----- downwards recursion for (t.le.3.0) -----
c                    recursion 38c in mcmurchie and davidson
               q(lmax,lmax)=qcomp(lmax,lmax,alpha,rk,t)
               do 60 l=lmax-1,0,-1
                  q(l,l)=(talph*q(l+2,l)-rk*q(l+1,l+1))
     $                       /float(l+l+1)
   60          continue
               if (lmax.le.3)  return
            else
c              ----- upwards recursion for (t.gt.3.0) -----
c                    recursion 38b in mcmurchie and davidson
               q(0,0)=qcomp(0,0,alpha,rk,t)
               do 80 l=1,lmax
                  q(l,l)=(talph*q(l+1,l-1)-float(l+l-1)*q(l-1,l-1))/rk
   80          continue
               if (lmax.le.3) return
            endif
c
         else if(nstart.eq.1) then
c           ----- n=1 terms in potential -----
            nbeg=3
            if (t.le.three) then
c              ----- downwards recursion for (t.le.3.0) -----
c              equation 38f in mcmurchie and davdison
               q(lmax+1,lmax)=qcomp(lmax+1,lmax,alpha,rk,t)
               q(lmax,lmax-1)=qcomp(lmax,lmax-1,alpha,rk,t)
               do 100 l=lmax-2,0,-1
                  q(l+1,l)=(talph*q(l+3,l+2)
     $                        -(rk-float(l+l+3)*talph/rk)*q(l+2,l+1))
     $                        /float(l+l+2)
  100         continue
            else
c              ----- upwards recursion for (t.gt.3.0) -----
c                     recursion 38e in mcmurchie and davidson
               q(1,0)=qcomp(1,0,alpha,rk,t)
               q(2,1)=qcomp(2,1,alpha,rk,t)
               do 120 l=2,lmax
                  q(l+1,l)=(float(l+l-2)*q(l-1,l-2)
     $                       +(rk-float(l+l-1)*talph/rk)*q(l,l-1))
     $                       /talph
  120          continue
            endif
c
         else if(nstart.eq.2) then
c           ----- n=2 terms in potential -----
c                 recursion 38a in mcmurchie and davidson
            q(2,0)=qcomp(2,0,alpha,rk,t)
            do 140 l=1,lmax
               q(l+2,l)=rk*q(l+1,l-1)/talph
  140       continue
            nbeg=4
         endif
c
c
         nhi=lmax+nstart
         lbeg=0
         do 180 n=nbeg,nhi
            do 160 l=lbeg,n-nbeg,2
c******        if this is 38d it looks like there is a typo in the paper.
               q(n,l)=(float(n+l-1)*q(n-2,l)
     $                   +rk*q(n-1,l+1))/talph
  160       continue
c           lbeg=and(lbeg+1,1)
            lbeg=mod(lbeg+1,2)
  180    continue
      endif
c
c
      return
      end
