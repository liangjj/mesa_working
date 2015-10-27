*deck @(#)fbgint.f	1.2  10/27/94
c***begin prologue     fbgint
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             @(#)fbgint.f	1.2 10/27/94 
c***purpose            greens function quadrature
c***
c***description        performs the integration of the radial
c***                   greens function and an inhomogeneity.
c               
c***references
c
c***routines called
c
c***end prologue       fbgint
c
      subroutine fbgint(flm,psilm,j,y,wt,pint,lmax,m,n,prnt)
      implicit integer (a-z)
      real*8 flm, psilm, wt, j, y, sumf, sumb, pint
      dimension flm(n,m:lmax), wt(n), psilm(n,m:lmax), j(n,0:lmax)
      dimension y(n,0:lmax), pint(n)
      logical prnt
      common /io/ inp, iout
      do 10 l=m,lmax
c----------------------------------------------------------------------c
c                in the following section of code we compute           c
c                the forward part of the integration                   c
c                                                                      c
c            f(i) = y(i) * sum ( k=1..i ) j(k)*flm(k)*wt(k)            c
c                                                                      c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c             initialize the integral at zero                          c
c----------------------------------------------------------------------c      
         sumf=0.d0
c----------------------------------------------------------------------c
c           step up on the integral and get it at all the remaining    c
c           points.                                                    c
c----------------------------------------------------------------------c
         do 20 i=1,n
            sumf=sumf+wt(i)*j(i,l)*flm(i,l)
            pint(i)=sumf
            psilm(i,l)=y(i,l)*sumf
   20    continue
         if (prnt) then
             write(iout,1) m, l, (pint(i),i=1,n)
         endif
c----------------------------------------------------------------------c
c                we are done with the forward part. begin the          c
c                backward sum                                          c
c                                                                      c
c            b(i) = j(i) * sum ( k=i+1...n ) y(k)*flm(k)*wt(k)         c
c                                                                      c
c                when we are done we have in psilm                     c
c                                                                      c
c            psilm(i) = y(i) * sum ( k=1..i ) j(k)*flm(k)*wt(k) +      c
c                       j(i) * sum ( k=i+1...n ) y(k)*flm(k)*wt(k)     c
c            which is the radial wavefunction at all of the quadrature c
c            points. this does not include the origin for a gauss      c
c                    type quadrature.                                  c
c            thus we have constructed a quadrature approximation to;   c
c                                                                      c
c            psilm(x) = y(x) * integral(0,x) [ j(x) * flm(x) ] +       c
c                       j(x) * integral(x,infinity) [ y(x)*flm(x0 ]    c
c                                                                      c
c----------------------------------------------------------------------c
         sumb=0.d0
         do 30 i=n-1,1,-1
            iplus=i+1
            sumb=sumb+wt(iplus)*y(iplus,l)*flm(iplus,l)
            pint(i)=sumb
            psilm(i,l)=psilm(i,l)+sumb*j(i,l)
   30    continue
         if (prnt) then
             write(iout,2) m, l, (pint(i),i=n-1,1,-1)
         endif
   10 continue
      
      return
    1 format(/,'forward partial integral for m = ',i3,1x,'l = ',i3,
     1          (/,5(e15.8,1x)))
    2 format(/,'backward partial integral for m = ',i3,1x,'l = ',i3,
     1          (/,5(e15.8,1x)))
      end















