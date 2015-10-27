*deck setd4.f
c***begin prologue     setd4
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form four dimensional index array.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       setd4
      subroutine setd4(ind,n1,n2,n3,n4,n)
      implicit integer (a-z)
      dimension ind(n1,n2,n3,n4)
      common/io/inp, iout
      cnt=0
      do 10 i=1,n4
         do 20 j=1,n3
            do 30 k=1,n2
               do 40 l=1,n1
                  ind(l,k,j,i)=cnt
 40            continue   
 30         continue
 20      continue
 10   continue 
      return
      end       

