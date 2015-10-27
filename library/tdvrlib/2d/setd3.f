*deck setd3.f
c***begin prologue     setd3
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form index array for 3D hamiltonian with channels.
c***                   one of the dimensions may be time
c***                   
c***references         
c
c***routines called    
c***end prologue       setd3
      subroutine setd3(ind,n1,n2,n3,nc)
      implicit integer (a-z)
      dimension ind(n1,n2,n3,nc)
      common/io/inp, iout
      cnt=0
      do 10 ic=1,nc
         do 20 i=1,n3
            do 30 j=1,n2
               do 40 k=1,n1
                  ind(k,j,i,ic)=cnt
 40            continue   
 30         continue
 20      continue
 10   continue 
      return
      end       

