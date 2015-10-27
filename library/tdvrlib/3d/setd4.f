*deck setd4.f
c***begin prologue     setd4
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form index array for 4D hamiltonian with channels.
c***                   one of the dimensions may be time.
c***                   
c***references         
c
c***routines called    
c***end prologue       setd4
      subroutine setd4(ind,n1,n2,n3,n4,nc)
      implicit integer (a-z)
      dimension ind(n1,n2,n3,n4,nc)
      common/io/inp, iout
      cnt=0
      do 10 ic=1,nc
         do 20 i=1,n4
            do 30 j=1,n3
               do 40 k=1,n2
                  do 50 l=1,n1
                     ind(l,k,j,i,ic)=cnt
 50               continue   
 40            continue   
 30         continue
 20      continue
 10   continue 
      return
      end       

