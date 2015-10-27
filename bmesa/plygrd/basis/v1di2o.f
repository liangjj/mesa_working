*deck v1di2o.f
c***begin prologue     v1di2o
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***routines called    
c***end prologue       v1di2o
      subroutine v1di2o(vout,vin,p,n,nvc)
      implicit integer (a-z)
      real*8 vout, vin, p, fac
      dimension vout(n,nvc), vin(n,nvc), p(n,n)
      common/io/inp, iout
      do 10 i=1,n
         fac = 1.d0/p(i,i)
         do 20 j=1,nvc
            vout(i,j) = vin(i,j)*fac
 20      continue
 10   continue
      return
      end









