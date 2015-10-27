*deck fock.f
c***begin prologue     fock
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             tstdiis
c***purpose            form current approximation to fock matrix
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       fock
      subroutine fock(f,ham0,vnl,n,iter,prnt)
      implicit integer (a-z)
      real*8 f, ham0, vnl
      character*3 itoc
      character*80 title
      logical prnt
      dimension f(n,n), ham0(n,n), vnl(n)
      common/io/inp, iout
c    
c     copy trap hamiltonian to fock matrix
c
      call copy(ham0,f,n*n)
c
c     add in non-linear term
c
      do 10 i=1,n
         f(i,i) = f(i,i) + vnl(i)    
 10   continue   
      if(prnt) then
         write(iout,1) iter
         title='fock matrix iteration = '//itoc(iter)
c         call prntrm(title,f,n,n,n,n,iout)
         call prntfm(title,f,n,n,n,n,iout)
      endif
      return 
   1  format(/,5x,'fock matrix formation iteration =', i3)
      end       


