*deck setnd1.f
c***begin prologue     setnd1
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form one dimensional index array.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       setnd1
      subroutine setnd1(ind,n1,n)
      implicit integer (a-z)
      dimension ind(n,1)
      common/io/inp, iout 
      do 10 i=1,n1
         ind(i,1)=i
 10   continue   
      return
      end       

