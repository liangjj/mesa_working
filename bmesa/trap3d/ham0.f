*deck ham0.f
c***begin prologue     ham0
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian, one-particle
c***author             schneider, barry (nsf)
c***source             
c***purpose            one particle hamiltonian, eigenvalues
c***                                and eigenvectors.
c***                   
c***description        diagonalize one particle hamiltonian in 
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
c***end prologue       ham0
      subroutine ham0(p,dp,ddp,wt,hamil,v,u,eig,tmp,hbar,mass,n,
     1                npts,prn)
c                             ,dir)
      implicit integer (a-z)
      real*8 p, dp, ddp, wt, eig, hamil, v, u, hbar, mass
      real*8 tmp, tfac
      logical prn
c     logical dir
      character*80 title
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1)
      dimension wt(npts), hamil(n,n), v(n,*)
      dimension u(n,n), eig(n), tmp(*) 
      common/io/inp, iout
      call rzero(hamil,n*n)
      tfac = hbar*hbar*.5d0/mass
      do 10 i=1,n
         do 20 j=1,i
            do 30 k=1,npts
               hamil(i,j) = hamil(i,j) - wt(k)*p(k,i-1)*ddp(k,j-1)
 30         continue
            hamil(i,j) = hamil(i,j) + p(npts,i-1)*dp(npts,j-1) 
     1                      - p(1,i-1)*dp(1,j-1)
            hamil(i,j) = tfac*hamil(i,j)
            hamil(j,i) = hamil(i,j)
 20      continue   
 10   continue
      call copy(hamil,u,n*n)
c      if(dir) then
c         do 40 i=1,n
c            do 50 j=1,n
c               u(i,j) = u(i,j) + v(i,j)
c 50         continue   
c 40      continue
c      else
       do 60 i=1,n
          u(i,i) = u(i,i) + v(i,1)
 60    continue
c      endif   
      if (prn) then
          title='unperturbed hamiltonian'
          call prntfm(title,u,n,n,n,n,iout)
      endif
      call rdiag(u,eig,tmp,n)
      return
      end       
