*deck @(#)addlag.f	1.1  11/30/90
      subroutine addlag(talag,ta,xlag,nob,noc,ndf)
      implicit integer(a-z)
c
      real*8 talag(*),ta(*),xlag(*)
c
c     talag(j,i) = talag(j,i) + ta(j,k)*xlag(k,i)
c
      ix=1
      jx=1
      nobnob=nob*nob
      nocnob=noc*nob
c
      do 10 i=1,ndf
c
         call apbtc(talag(jx),ta(ix),xlag,nob,nob,noc)
c
c1       write(6,7001)
7001     format(/,' ta_lagrangian ')
c1       call matout(talag(jx),nob,noc,nob,noc,6)
c
         ix=ix+nobnob
         jx=jx+nocnob
c
  10  continue
c
      return
      end
