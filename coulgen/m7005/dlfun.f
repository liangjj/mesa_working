*deck dlfun
      function dlfun(eta,l,fact)
c***begin prologue     dlfun
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***description        coulomb d  function
c                               l
c                      this function is trivially related to c
c                                                             l
c***references         NBS handbook
c
c***routines called    clfun(code)
c
c***end prologue       dlfun
c
      implicit integer (a-z)
      real*8 dlfun, eta, clfun, one, fact
      dimension fact(*)
      common/io/inp,iout
      data one/1.d0/
      dlfun=one/( (l+l+1)*clfun(eta,l,fact) )
      return
      end


