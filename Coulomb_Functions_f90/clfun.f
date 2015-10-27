*deck  cl_fun
  REAL*8 Function cl_fun()
c***begin prologue     m2705
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c
c
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            compute c  function for positive energy coulomb
c***                            l                           functions. 
c***description
c***references         NBS handbook
c
c***routines called    cgamma(complex gamma function:clams)
c
c***end prologue       m2705
c
      implicit integer (a-z)
c
      real*8 eta, pi, fact, clfun
      real*8 zero, half, one, two
      complex *16 arg, cgamma, eye
      dimension fact(0:100)
      common/io/inp,iout
      data pi /3.14159265358979323846d+00/
      data zero, half, one, two/0.d0,.5d0,1.d0,2.d0/
      data eye/(0.d0,1.d0)/
c
c
      clfun=(two**l)*exp(-pi*eta*half)
      arg=l+1+eye*eta
      clfun=clfun*abs(cgamma(arg))/fact(2*l+1)
      return
      end
