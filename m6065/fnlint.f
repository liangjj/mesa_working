*deck @(#)fnlint.f
c***begin prologue     fnlint
c***date written       xxxxxx   (yymmdd)
c***revision date      920409   (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (nsf)
c***source             m6060
c***purpose            optical potential integral
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       fnlint
      subroutine fnlint(sinnt,energy,eigval,nolam)
      implicit integer (a-z)
      real *8 energy 
      complex *16  sinnt, eigval, val
      dimension sinnt(nolam), eigval(nolam)
      common /io/ inp, iout
      val=dcmplx(0.d0,0.d0)
      do 10 lam=1,nolam
         val=val + sinnt(lam)*sinnt(lam)/(energy-eigval(lam))
   10 continue
      write(iout,20) val
   20 format(//,15x,'value of test integral',1x,e15.8,1x,e15.8)  
      return
      end





