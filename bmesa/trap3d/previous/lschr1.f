*deck lschr1.f
c***begin prologue     lschr1
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           schroedinger equation
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for one dimensional schroedinger equation.
c***                   
c***references         
c
c***routines called    
c***end prologue       lschr1
      subroutine lschr1(ham,eig,v,n)
      implicit integer (a-z)
      real*8 ham, eig, v
      character*80 title
      dimension ham(n,n), eig(n), v(n)
      common/io/inp, iout
      write(iout,1)      
      do 10 i=1,n
         ham(i,i) = ham(i,i) + v(i) 
 10   continue   
      call tred2(n,n,ham,eig,v,ham)
      call tql2(n,n,eig,v,ham,ierr)
      return
 1    format(/,5x,'diagonalize time-independent one dimensional '
     1            'schroedinger equation')
      end       
