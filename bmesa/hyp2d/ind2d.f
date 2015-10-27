*deck ind2d.f
c***begin prologue     ind2d
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form two dimensional index array.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       ind2d
      subroutine ind2d(ind,n1,n2)
      implicit integer (a-z)
      dimension ind(n2,n1)
      common/io/inp, iout
      cnt=0 
      do 10 i=1,n1
         do 20 j=1,n2
            cnt=cnt+1
            ind(j,i)=cnt
 20      continue
 10   continue   
      return
      end       

