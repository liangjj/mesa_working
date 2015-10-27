*deck setnd3.f
c***begin prologue     setnd3
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form three dimensional index array.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       setnd3
      subroutine setnd3(ind,n1,n2,n3,n)
      implicit integer (a-z)
      dimension ind(n,3)
      common/io/inp, iout
      cnt=0
      do 10 i=1,n1
         do 20 j=1,n2
            do 30 k=1,n3
               cnt=cnt+1
               ind(cnt,1)=i
               ind(cnt,2)=j
               ind(cnt,3)=k
 30         continue
 20      continue
 10   continue 
      return
      end       

