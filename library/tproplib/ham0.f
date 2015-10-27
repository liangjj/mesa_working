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
      subroutine ham0(p,dp,ddp,q,wt,hamil,v,u,eig,
     1                vec,work,hbar,mass,n,lftbc,rtbc,
     2                coord,comp,prn)
      implicit integer (a-z)
      real*8 p, dp, ddp, q, wt, eig, hamil, v, u 
      real*8 vec, hbar, mass, tfac
      logical prn
      character*80 title
      character*(*) coord, comp
      dimension p(n,n), dp(n,n), ddp(n,n)
      dimension q(n), wt(n), hamil(n,n), v(n)
      dimension u(n,n), eig(n), work(*) 
      dimension vec(n,n) 
      common/io/inp, iout
      ilft=0
      irt=0
      if(lftbc.ne.0) then
         ilft=1
      endif
      if(rtbc.ne.0) then
         irt=1
      endif              
      write(iout,1) coord, comp
      call rzero(hamil,n*n)
      tfac = hbar*hbar*.5d0/mass
      if(comp.eq.'x'.or.comp.eq.'y'.or.comp.eq.'z') then
         call tcart(p,dp,ddp,wt,hamil,chamil,ilft,irt,n)
      elseif(comp.eq.'rho') then
c
c        for cylindrical rho coordinate
c
         call tcylin(p,dp,ddp,q,wt,hamil,ilft,irt,n)
      elseif(comp.eq.'phi') then
c
c        for angular phi coordinate.  form is same as cartesian, just different
c        end points
c
         call tcart(p,dp,ddp,wt,hamil,ilft,irt,n)
      elseif(comp.eq.'r') then
c 
c        radial coordinate
c
         call tcart(p,dp,ddp,wt,hamil,ilft,irt,n)
      elseif(comp.eq.'theta') then
c        
c        spherical theta coordiante
c
         call tthet(p,dp,ddp,q,wt,hamil,work,ilft,irt,n)
      else
         call lnkerr('error in coordinate type')
      endif
      call vscale(hamil,hamil,tfac,n*n)
      do 20 i=1,n
         hamil(i,i) = hamil(i,i) + v(i)
 20   continue
      call copy(hamil,u,n*n)
      if (prn) then
          title='unperturbed hamiltonian'
          call prntrm(title,u,n,n,n,n,iout)
      endif
      call rdiag(u,eig,vec,work,n,n)
      return
 1    format(/,5x,'forming a component of the unperturbed hamiltonian',
     1       /,5x,'coordinate system = ',a16,/,5x,'coordinate = ',a8)
      end       





