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
      subroutine newfn(f,ham0,vnl,eig,psi,n)
      implicit integer (a-z)
      real*8 f, ham0, vnl, eig, psi
      dimension f(n,n), ham0(n,n), vnl(*), eig(n), psi(n)
      common/io/ inp, iout
      call copy(ham0,f,n*n)
      do 10 i=1,n
         f(i,i) = f(i,i) + vnl(i)
 10   continue   
      call tred2(n,n,f,eig,psi,f)
      call tql2(n,n,eig,psi,f,ierr)
      call copy(f,psi,n)
      return
      end       

