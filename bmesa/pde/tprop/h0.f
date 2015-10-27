*deck h0.f
c***begin prologue     h0
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
c***end prologue       h0
      subroutine h0(p,dp,ddp,q,wt,hamil,v,u,eig,vec,work,
     1              hbar,mass,n,lftbc,rtbc,coord,comp,prn)
      implicit integer (a-z)
      real*8 p, dp, ddp, q, wt, eig, hamil, v, u 
      real*8 vec, hbar, mass, work, tfac
      logical prn
      character*16 mattyp
      character*80 title
      character*(*) coord, comp
      dimension p(n,n), dp(n,n), ddp(n,n)
      dimension q(n), wt(n), hamil(n,n), v(n)
      dimension u(n,n), eig(n), vec(n,n), work(*) 
      data mattyp/'real-symmetric'/
      data ldum/0/
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
         call tcart(p,dp,ddp,wt,hamil,hamil,ilft,irt,n,mattyp)
      elseif(comp.eq.'rho') then
c
c        for cylindrical rho coordinate
c
         call tcylin(p,dp,ddp,q,wt,hamil,hamil,ilft,irt,n,mattyp)
      elseif(comp.eq.'phi') then
c
c        for angular phi coordinate.  form is same as cartesian, just different
c        end points
c
         call tcart(p,dp,ddp,wt,hamil,hamil,ilft,irt,n,mattyp)
      elseif(comp.eq.'r') then
c 
c        radial coordinate
c
         call tcart(p,dp,ddp,wt,hamil,hamil,ilft,irt,n,mattyp)
      elseif(comp.eq.'theta') then
c        
c        spherical theta coordiante
c
         call tthet(p,dp,ddp,q,wt,hamil,hamil,work,ilft,irt,
     1              n,mattyp)
      else
         call lnkerr('error in coordinate type')
      endif
      call vscale(hamil,hamil,tfac,n*n)
      do 10 i=1,n
         hamil(i,i) = hamil(i,i) + v(i)
 10   continue
      call copy(hamil,u,n*n)
      if (prn) then
          title='unperturbed hamiltonian'
          call prntrm(title,u,n,n,n,n,iout)
      endif
      call rdiag(u,u,eig,eig,vec,vec,vec,
     1           work,ldum,work,n,n,mattyp)
      return
 1    format(/,5x,'forming a component of the unperturbed hamiltonian',
     1       /,5x,'coordinate system = ',a16,/,5x,'coordinate = ',a8)
      end       





