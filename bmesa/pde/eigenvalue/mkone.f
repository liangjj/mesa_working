*deck mkone.f
c***begin prologue     mkone
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           1-d hamiltonian
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate one body hamiltonian matrix.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       mkone
      subroutine mkone(h,v,n)
      implicit integer (a-z)
      real*8 h, v
      dimension h(n,n), v(n) 
      common/io/inp, iout
      do 10 i=1,n
         h(i,i) = h(i,i) + v(i)
 10   continue   
      return         
      end       
