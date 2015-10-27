*deck v2di2o.f
c***begin prologue     v2di2o
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***routines called    
c***end prologue       v2di2o
      subroutine v2di2o(vout,vin,p1,p2,n1,n2,nvc)
      implicit integer (a-z)
      real*8 vout, vin, p1, p2, fac
      dimension vout(n2,n1,nvc), vin(n2,n1,nvc)
      dimension p1(n1,n1), p2(n2,n2)
      common/io/inp, iout
      do 10 i=1,n1
         do 20 j=1,n2
            fac = 1.d0/(p1(i,i)*p2(j,j))
            do 30 k=1,nvc
               vout(j,i,k) = vin(j,i,k)*fac
 30         continue
 20      continue
 10   continue
      return
      end









