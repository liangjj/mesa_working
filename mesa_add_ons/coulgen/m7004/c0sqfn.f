*deck c0sqfn
      function c0sqfn(eta)
c***begin prologue     c0sqfn
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c                                         2
c***description        coulomb function c
c                                        0
c***references         NBS handbook
c
c***routines called
c
c***end prologue       c0sqfn
c
      implicit integer (a-z)
      real*8 c0sqfn, eta, pi, twpiet, one, two
      data pi /3.14159265358979323846d+00/
      data one, two/1.d0,2.d0/
      common/io/inp,iout
      twpiet=two*pi*eta
      c0sqfn=twpiet/(exp(twpiet)-one)
      return
      end
