*deck v4di2o.f
c***begin prologue     v4di2o
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***routines called    
c***end prologue       v4di2o
      subroutine v4di2o(vout,vin,p1,p2,p3,p4,n1,n2,n3,n4,nvc)
      implicit integer (a-z)
      real*8 vout, vin, p1, p2, p3, p4, fac
      dimension vout(n4,n3,n2,n1,nvc), vin(n4,n3,n2,n1,nvc)
      dimension p1(n1,n1), p2(n2,n2), p3(n3,n3), p4(n4,n4)
      common/io/inp, iout
      do 10 i=1,n1
         do 20 j=1,n2
            do 30 k=1,n3
               do 40 l=1,n4
                  fac = 1.d0/(p1(i,i)*p2(j,j)*p3(k,k)*p4(l,l))
                  do 50 m=1,nvc
                     vout(l,k,j,i,m) = vin(l,k,j,i,m)*fac
 50               continue   
 40            continue   
 30         continue
 20      continue
 10   continue
      return
      end









