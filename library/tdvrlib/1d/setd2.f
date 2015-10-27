*deck setd2.f
c***begin prologue     setd2
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form index array for a 2D hamiltonian with channels.
c***                   one of the dimensions may be time.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       setd2
      subroutine setd2(ind,n1,n2,nc)
      implicit integer (a-z)
      dimension ind(n1,n2,nc)
      common/io/inp, iout
      cnt=0 
      do 10 ic=1,nc
         do 20 i=1,n2
            do 30 j=1,n1
               cnt=cnt+1
               ind(j,i,ic)=cnt
 30         continue   
 20      continue
 10   continue   
      return
      end       





