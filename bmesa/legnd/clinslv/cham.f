*deck cham.f
c***begin prologue     cham
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian, one-particle
c***author             schneider, barry (nsf)
c***source             
c***purpose            one particle hamiltonian.
c***                   
c***description        construct one particle hamiltonian in 
c***                   polynomial basis.  the hamiltonian is
c                             2    2
c***                  ( - hbar )  d        
c                     (   ___  )  -    + v
c                     (   2*m  )    2  
c                                 dr    
c
c***                   v can be one of serveral interactions.  see
c***                   the if statements in the code
c***references         
c
c***routines called    
c***end prologue       cham
      subroutine cham(p,dp,ddp,q,wt,hamil,eig,vecr,v,work,n,
     1                lftbc,rtbc,diag,prn)
      implicit integer (a-z)
      real*8 p, dp, ddp, q, wt
      complex*16 hamil, v
      logical prn, diag
      character*80 title
      complex*16 eig, vecr
      real*8 work
      dimension p(n,*), dp(n,*), ddp(n,*)
      dimension q(n), wt(n), hamil(n,n), v(n)
      dimension eig(n), vecr(n,n), work(3*n)
      common/io/inp, iout
      ilft=0
      irt=0
      if(lftbc.ne.0) then
         ilft=1
      endif
      if(rtbc.ne.0) then
         irt=1
      endif              
      call czero(hamil,n*n)
      call tcart(p,dp,ddp,wt,hamil,ilft,irt,n,nonsym)
      call cvscal(hamil,hamil,.5d0,n*n)
      do 10 i=1,n
         hamil(i,i) = hamil(i,i) + v(i)
 10   continue
c      title='hamiltonian'
c      call prntcm(title,hamil,n,n,n,n,iout)
      if(diag) then
         call cgeev(hamil,n,n,eig,vecr,n,work,1,ierr)
         title='eigenvalues'
         call prntcm(title,eig,n,1,n,1,iout)
      endif         
      return
      end       





