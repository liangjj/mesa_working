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
      subroutine ham0(p,dp,ddp,q,wt,hamil,v,u,eig,tmp,t,hbar,mass,n,
     1                lftbc,rtbc,coord,comp,prn,parity,type,
     2                neven,nodd)
      implicit integer (a-z)
      real*8 p, dp, ddp, q, wt, eig, hamil, v, u, hbar, mass
      real*8 tmp, t, tfac
      logical prn, parity
      character*80 title
      character*(*) coord, comp, type
      dimension p(n,*), dp(n,*), ddp(n,*)
      dimension q(n), wt(n), hamil(n,n), v(*)
      dimension u(n,n), eig(n), tmp(*), t(*) 
      common/io/inp, iout
c
c     if the states are not good parity states the entire hamiltonian
c     is diagonalized.  this will be true even if parity is a good
c     quantum number but is not declared as such on the input data.
c     if parity is declared on the input, then either the even or the odd
c     parity states will be determined, not both.  this is done using the
c     variable type which is passed in to the routine.  the other parity
c     matrices and basis functions are not destroyed but sit behind the
c     matrices for the desired parity.  the output of this routine always
c     contains the matrices and vectors of the parity desired, stored 
c     first and as if the other information was not required.  In other words
c     all matrices are packed with no extraneous zeros and the vectors of the
c     desired parity are first in the list.
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
         call tcart(p,dp,ddp,wt,hamil,ilft,irt,n)
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
         call tthet(p,dp,ddp,q,wt,hamil,tmp,ilft,irt,n)
      else
         call lnkerr('error in coordinate type')
      endif
      call vscale(hamil,hamil,tfac,n*n)
      do 10 i=1,n
         hamil(i,i) = hamil(i,i) + v(i)
 10   continue
      call copy(hamil,u,n*n)
      if(parity) then
c
c        transform to the states of good parity
c      
         call trevod(t,neven,nodd,n)
         call ebc(tmp,hamil,t,n,n,n)
         call ebtc(hamil,t,tmp,n,n,n)
         call ebc(tmp,u,t,n,n,n)
         call ebtc(u,t,tmp,n,n,n)
         call ebtc(tmp,p,t,n,n,n)
         call copy(tmp,p,n*n)
         call ebtc(tmp,dp,t,n,n,n)
         call copy(tmp,dp,n*n)
         call ebtc(tmp,ddp,t,n,n,n)
         call copy(tmp,ddp,n*n)
         call pkmat(hamil,tmp,n,neven,nodd,type,'matrix')
         call pkmat(u,tmp,n,neven,nodd,type,'matrix')
         call pkmat(p,tmp,n,neven,nodd,type,'vector')
         call pkmat(dp,tmp,n,neven,nodd,type,'vector')
         call pkmat(ddp,tmp,n,neven,nodd,type,'vector')
      endif
      if (prn) then
          title='unperturbed hamiltonian'
          call prntfm(title,u,n,n,n,n,iout)
      endif
      if(parity) then
         if(type.eq.'even') then
            call rdiag(u,eig,tmp,neven,neven)
         elseif(type.eq.'odd') then
            call rdiag(u,eig,tmp,nodd,nodd)
         endif
      else
         call rdiag(u,eig,tmp,n,n)
      endif
      return
 1    format(/,5x,'forming a component of the unperturbed hamiltonian',
     1       /,5x,'coordinate system = ',a16,/,5x,'coordinate = ',a8)
      end       
