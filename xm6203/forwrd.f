*deck @(#)forwrd.f	1.2  10/27/94
      subroutine forwrd(flm,j,y,wt,int0,psilm,temp,n)
c***begin prologue     forwrd
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             @(#)forwrd.f	1.2 10/27/94 
c***purpose            forward integration of indefinite integral
c***
c***description        performs the set of indefinite integrals which result
c***                   from integrating a function on [a,b] as [a,r(1)], [a,r(2)]
c                      [a,r(3)]....[a,b] where r(i) are the quadrature points in
c                      the interval. the integral is initialized as int0
c                      which is either its last value or zero depending on the
c                      interval. psilm(1)=[a,r1] and psilm(n-1)=[a,b]   
c               
c***references
c
c***routines called
c
c***end prologue       forwrd
c
      implicit integer (a-z)
      dimension flm(n), j(n), y(n), psilm(n-1), wt(n,n-1), temp(n-1)
      common /io/ inp, iout
      real*8 flm, j, y, psilm, wt, int0, sumf, temp
      sumf=int0
      do 10 pti=1,n-1
         do 20 ptj=1,n
            sumf=sumf+wt(ptj,pti)*flm(ptj)*j(ptj)
   20    continue
         psilm(pti)=sumf
   10 continue
      call copy(psilm,temp,n-1)
      int0=psilm(n-1)
      do 30 pti=1,n-1
         psilm(pti)=psilm(pti)*y(pti+1)
   30 continue
      return
      end















