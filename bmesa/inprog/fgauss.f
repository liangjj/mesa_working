*deck fgauss
c***begin prologue     fgauss
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            forward integration of integral equation
c***
c***description        performs the forward integration of the radial
c***                   integral equation involving a greens function and
c                      an inhomogeneity.
c               
c***references
c
c***routines called
c
c***end prologue       fgauss
c
      subroutine fgauss(flm,psilm,oldval,j,y,wt,lndex,nl,ltop,n)
      implicit integer (a-z)
      real*8 flm, psilm, wt, j, y, oldval
      dimension flm(n,nl), wt(n), psilm(n,nl), j(n,0:ltop), y(n,0:ltop)
      dimension lndex(nl), oldval(nl)
      common /io/ inp, iout
c----------------------------------------------------------------------c
c             initialize the integral at its first point using         c
c             the old value.                                           c
c----------------------------------------------------------------------c
      do 10 l=1,nl
         psilm(1,l)=oldval(l)+wt(1)*j(1,lndex(l))*flm(1,l)
   10 continue
c----------------------------------------------------------------------c
c           step up on the integral and get it at all the remaining    c
c           points.                                                    c
c----------------------------------------------------------------------c
      do 20 l=1,nl
         do 30 i=2,n
            psilm(i,l)=psilm(i-1,l)+wt(i)*j(i,lndex(l))*flm(i,l)
   30    continue
   20 continue
c----------------------------------------------------------------------c
c          store the last value in oldval so it can be used to         c
c          initialize the integral at the next call                    c
c----------------------------------------------------------------------c
      do 40 l=1,nl
         oldval(l)=psilm(n,l)
   40 continue
c----------------------------------------------------------------------c
c          now make the actual function by multiplying the intgral     c
c          by the prefactor in front of it.                            c
c----------------------------------------------------------------------c
      do 50 l=1,nl
         call vmul(psilm(1,l),psilm(1,l),y(1,lndex(l)),n)
   50 continue
      return
      end















