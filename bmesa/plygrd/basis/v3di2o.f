*deck v3di2o.f
c***begin prologue     v3di2o
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***routines called    
c***end prologue       v3di2o
      subroutine v3di2o(vout,vin,p1,p2,p3,n1,n2,n3,nvc)
      implicit integer (a-z)
      real*8 vout, vin, p1, p2, p3, fac
      dimension vout(n3,n2,n1,nvc), vin(n3,n2,n1,nvc)
      dimension p1(n1,n1), p2(n2,n2), p3(n3,n3)
      common/io/inp, iout
      do 10 i=1,n1
         do 20 j=1,n2
            do 30 k=1,n3
               fac = 1.d0/(p1(i,i)*p2(j,j)*p3(k,k))
               do 40 l=1,nvc
                  vout(k,j,i,l) = vin(k,j,i,l)*fac
 40            continue   
 30         continue
 20      continue
 10   continue
      return
      end









