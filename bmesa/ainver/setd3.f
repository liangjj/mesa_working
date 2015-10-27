*deck setd3.f
c***begin prologue     setd3
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
c***end prologue       setd3
      subroutine setd3(ind,n1,n2,n3,n)
      implicit integer (a-z)
      dimension ind(n1,n2,n3)
      common/io/inp, iout
      cnt=0
      do 10 i=1,n3
         do 20 j=1,n2
            do 30 k=1,n1
               ind(k,j,i)=cnt
 30         continue
 20      continue
 10   continue 
      return
      end       

