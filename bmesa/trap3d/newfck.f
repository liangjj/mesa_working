*deck newfck.f
c***begin prologue     newfck
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             tstdiis
c***purpose            form diis approximation to new fock matrix
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       newfck
      subroutine newfck(f,ham0,vnl,sol,n,iter)
      implicit integer (a-z)
      real*8 f, ham0, vnl, sol, csum
      dimension f(n,n), ham0(n,n), vnl(n,*), sol(iter)
      common/io/inp, iout
      csum=0.d0
      do 10 i=1,iter
         csum = csum + sol(i)
 10   continue
c
      do 20 i=1,n
         do 30 j=1,i
            f(i,j) = csum*ham0(i,j)
            f(j,i) = f(i,j)
 30      continue
 20   continue   
c
      do 40 i=1,n
         do 50 j=1,iter
            f(i,i) = f(i,i) + sol(j)*vnl(i,j)    
 50      continue   
 40   continue   
      return 
      end       


