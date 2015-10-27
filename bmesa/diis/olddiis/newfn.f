*deck newfn.f
c***begin prologue     newfn
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***description          
c***references         
c
c***routines called    
c***end prologue       newfn
      subroutine newfn(f,eig,psi,work,n)
      implicit integer (a-z)
      real*8 f, eig, psi, work
      dimension obj(n,n), eig(n), psi(n), work(n,n)
      common/io/ inp, iout
      call copy(f,work,n*n)
      call tred2(n,n,work,eig,psi,work)
      call tql2(n,n,eig,psi,work,ierr)
      call copy(work(1,1),psi,n)
      return
      end       

