*deck vcprow.f
c***begin prologue     vcprow
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           copy
c***author             schneider, barry (nsf)
c***source             
c***purpose            vector copying
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vcprow
      subroutine vcprow(vout,vin,ni,nj,n,m,nc)
      implicit integer (a-z)
      real*8 vout, vin
      dimension vout(n,m,nc,nc), vin(m)
      common/io/inp, iout
      do 10 i=1,m
         do 20 j=1,n
            vout(j,i,ni,nj) = vout(j,i,ni,nj) + vin(i)
 20      continue
 10   continue
      return   
      end       


