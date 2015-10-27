*deck @(#)qpasy.f	5.1  11/6/94
      function qpasy (n,la1,lb1,alp,xka1,xkb1,iflag)
      implicit real*8(a-h,o-z)
c
c     partially asymptotic form for q(n,la1,lb1).......
c     result includes a factor exp(-xkb*xkb/(4.*alpha)) to prevent overflow
      common/qstore/q(13,11),alpha,rk,t
      common/dfac/dfac(30)
c
c first set up xkb as largest.
      if (iflag.eq.3) go to 10
      xka=xka1
      xkb=xkb1
      la=la1
      lb=lb1
      go to 20
   10 xka=xkb1
      xkb=xka1
      la=lb1
      lb=la1
   20 continue
c set up parameters for qcomp using xkb.
      alpha=1.
      rk=xkb/sqrt(alp)
      t=rk*rk/4.
c
c now run power series using xka,  obtaining initial q(n,l)'s
c    from qcomp then recurring upwards.
      tk=xka*xka/(2.*alp)
c  j=0 term in sum.
      qold1=qcomp(n+la,lb)
      sum=qold1/dfac(la+la+3)
      if (tk.eq.0) go to 50
c  j=1 term in sum.
      nprime=n+la+2
      qnew=qcomp(nprime,lb)
      term=qnew*tk/dfac(la+la+5)
      sum=sum+term
c  j=2  term
      j=2
      nprime=n+la+j+j
      coe=1./dfac(la+la+3)
      go to 40
c  increment j for next term
   30 j=j+1
      nprime=n+la+j+j
c    compute j-2 factor to be absorbed into qold1,qold2.
      coe=tk/float((j-2)*(la+la+j+j-3))
   40 qold2=qold1*coe
      qold1=qnew*coe
      qnew=(t+float(nprime+nprime-5)/2.)*qold1+float((lb-nprime+4)*(lb
     1 +nprime-3))*qold2/4.
c    now correct qnew for j-1 and j factors.
      term=qnew*tk*tk/(float((j-1)*(la+la+j+j-1))*float(j*(la+la+j+j+1))
     1 )
      sum=sum+term
      if (abs(term/sum).gt.1.d-13) go to 30
   50 if (la.eq.0) prefac=1./sqrt(alp**(n+la+1))
      if (la.ne.0) prefac=(xka**la)/sqrt(alp**(n+la+1))
      qpasy=prefac*sum
      return
      end
