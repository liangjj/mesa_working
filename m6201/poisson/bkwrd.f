*deck bkwrd
      subroutine bkwrd(flm,j,y,wt,int0,psilm,temp,n)
c***begin prologue     bkwrd
c***date written       930922   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            backward integration of indefinite integral
c***
c***description        performs the set of indefinite integrals which result
c***                   from integrating a function on [a,b] as [a,r(1)], [a,r(2)]
c                      [a,r(3)]....[a,b] where r(i) are the quadrature points in
c                      the interval. the integral is initialized as intgl(n)
c                      which is either its last value or zero depending on the
c                      interval. intgl(n-1)=[r(n-1),b]   
c                                             
c***references
c
c***routines called
c
c***end prologue       bkwrd
c
      implicit integer (a-z)
      dimension flm(n), j(n), y(n), wt(n,n-1), psilm(n-1), temp(n-1)
      common /io/ inp, iout
      real*8 flm, j, y, wt, psilm, int0, sumb, temp
c----------------------------------------------------------------------c
c      note that there is one less integral than grid point and that   c
c      the ith backward integral is needed to construct the            c
c      wavefunction at the ith grid point. the do 30 loop shows        c
c      that explicitly.                                                c
c----------------------------------------------------------------------c
      sumb=int0
      do 10 pti=n-1,1,-1
         do 20 ptj=1,n
            sumb=sumb+wt(ptj,pti)*flm(ptj)*y(ptj)
   20    continue
         psilm(pti)=sumb
   10 continue
      call copy(psilm,temp,n-1)
      int0=psilm(1)
      do 30 pti=1,n-1
         psilm(pti)=psilm(pti)*j(pti)
   30 continue 
      return
      end















