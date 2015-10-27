*deck setnd2.f
c***begin prologue     setnd2
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
c***end prologue       setnd2
      subroutine setnd2(ind,n1,n2,n)
      implicit integer (a-z)
      dimension ind(n,2)
      common/io/inp, iout
      cnt=0 
      do 10 i=1,n2
         do 20 j=1,n1
            cnt=cnt+1
            ind(cnt,1)=i
            ind(cnt,2)=j
 20      continue
 10   continue   
      return
      end       

