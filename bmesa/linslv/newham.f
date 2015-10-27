*deck newham.f
c***begin prologue     newham
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            transformed hamiltonian
c***                   
c***description        
c***                     
c***references         
c
c***routines called    
c***end prologue       newham
      subroutine newham(ham,eig,v,n)
      implicit integer (a-z)
      real*8 ham, eig, v
      dimension ham(n,n), eig(n), v(n,n)
      common/io/inp, iout 
      call copy(v,ham,n*n)
      do 10 i=1,n
         ham(i,i) = ham(i,i) + eig(i)
 10   continue   
      return
      end       
