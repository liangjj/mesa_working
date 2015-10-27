*deck @(#)recur1.f	5.1 11/6/94
      subroutine recur1 (nstart,lmax)
      implicit real*8(a-h,o-z)
c
c     controls computation of q(n,l) by recurrence.
c
      integer and
      real*8 zero,one,two,three,four
      real*8 half
      common/qstore/q(13,11),alpha,rk,t
      common/argab/argab,expab
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      data nzero /0/
      save nzero
c
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00,three=3.0d+00)
      parameter (four=4.0d+00,half=0.5d+00)
c
c     ----- check for special case of zero argument -----
c           note the early return
      if (rk.eq.zero) then
         nmax=nstart+lmax
         argab=zero
         expab=one
         q(1,1)=half*sqrt(pi/alpha)
         q(2,1)=half/alpha
         do 10 i=2,nmax
            q(i+1,1)=half*float(i-1)*q(i-1,1)/alpha
   10    continue
         do 20 i=1,lmax
            do 20 j=i,nmax
               q(j+1,i+1)=zero
   20    continue
         return
      endif
c
c     ----- normal processing -----
      t=rk*rk/(four*alpha)
      if (lmax.eq.0) then
         q(nstart+1,1)=qcomp(nstart,0)
         return
      else if(lmax.eq.1) then
         q(nstart+1,1)=qcomp(nstart,0)
         q(nstart+2,2)=qcomp(nstart+1,1)
         return
      else
         talph=alpha+alpha
         if (nstart.eq.0) then
c           ----- n=0 terms in potential -----
            q(3,1)=qcomp(2,0)
            do 50 l=1,lmax-1
               q(l+3,l+1)=rk*q(l+2,l)/talph
   50       continue
            nbeg=4
            if (t.le.three) then
c              downwards recurrence for t .le. 3.
               q(lmax+1,lmax+1)=qcomp(lmax,lmax)
               lhi=lmax-1
               do 60 lr=nzero,lhi
                  l=lhi-lr
                  q(l+1,l+1)=(talph*q(l+3,l+1)-rk*q(l+2,l+2))
     $                       /float(l+l+1)
   60          continue
               if (lmax.le.3)  return
            else
c              upwards recurrence for t .gt. 3.
               q(1,1)=qcomp(0,0)
               do 80 l=1,lmax
                  q(l+1,l+1)=(talph*q(l+2,l)-float(l+l-1)*q(l,l))/rk
   80          continue
               if (lmax.le.3) return
            endif
c
         else if(nstart.eq.1) then
c           ----- n=1 terms in potential -----
            nbeg=3
            if (t.le.three) then
c              downwards recurrence for t .le. 3.
               q(lmax+2,lmax+1)=qcomp(lmax+1,lmax)
               q(lmax+1,lmax)=qcomp(lmax,lmax-1)
               lhi=lmax-2
               do 100 lr=nzero,lhi
                  l=lhi-lr
                  q(l+2,l+1)=(talph*q(l+4,l+3)
     $                        -(rk-float(l+l+3)*talph/rk)*q(l+3,l+2))
     $                        /float(l+l+2)
  100         continue
            else
c              upwards recurrence for t .gt. 3.
               q(2,1)=qcomp(1,0)
               q(3,2)=qcomp(2,1)
               do 120 l=2,lmax
                  q(l+2,l+1)=(float(l+l-2)*q(l,l-1)
     $                       +(rk-float(l+l-1)*talph/rk)*q(l+1,l))
     $                       /talph
  120          continue
            endif
c
         else if(nstart.eq.2) then
c           ----- n=2 terms in potential -----
            q(3,1)=qcomp(2,0)
            do 140 l=1,lmax
               q(l+3,l+1)=rk*q(l+2,l)/talph
  140       continue
            nbeg=4
         endif
c
c
         nhi=lmax+nstart
         lbeg=0
         do 180 n=nbeg,nhi
            lhi=n-nbeg
            do 160 l=lbeg,lhi,2
               q(n+1,l+1)=(float(n+l-1)*q(n-1,l+1)
     $                   +rk*q(n,l+2))/talph
  160       continue
            lbeg=and(lbeg+1,1)
c           lbeg=(lbeg+1).and.1
  180    continue
      endif
c
c
      return
      end
