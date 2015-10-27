*deck v2od.f
c***begin prologue     v2od
c***date written       970609   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           potential, two-dimension
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       v2od
      subroutine v2od(v,parts,n,n1,n2)
      implicit integer (a-z)
      real*8 v, parts
      dimension v(*), parts(n,2)
      common/io/inp, iout
      count=0
      do 10 i=1,n1
         ii=n2*(i-1)
         do 20 j=1,n2
            ij=ii+j
            count=count+1
            v(count) = v(count) + parts(ij,1)
 20      continue
 10   continue
      return
      end       
