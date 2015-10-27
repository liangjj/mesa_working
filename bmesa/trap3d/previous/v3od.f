*deck v3od.f
c***begin prologue     v3od
c***date written       970609   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           potential, three-dimension
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       v2od
      subroutine v3od(v,parts,n,nd)
      implicit integer (a-z)
      real*8 v, parts
      dimension v(*), parts(n,2), nd(2)
      common/io/inp, iout
      count=0
      do 10 i=1,nd(1)
         i1=nd(2)*(i-1)
         i2=nd(3)*(i-1)
         do 20 j=1,nd(2)
            i3=nd(3)*(j-1)
            ij=i1+j
            do 30 k=1,nd(3)
               ik=i2+k 
               jk=i3+k
               count=count+1
               v(count) = v(count) + parts(ij,1)
     1                             + parts(ik,2)
     2                             + parts(jk,3)
 30         continue
 20      continue
 10   continue
      return
      end       
