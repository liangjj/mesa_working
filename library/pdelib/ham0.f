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
      subroutine ham0(p,dp,ddp,q,wt,hamil,chamil,v,u,cu,eig,
     1                ceig,vec,vecl,vecr,cvecl,cvecr,work,lwork,rwork,
     2                hbar,mass,n,lftbc,rtbc,coord,comp,prn,mattyp)
      implicit integer (a-z)
      real*8 p, dp, ddp, q, wt, eig, hamil, v, u 
      real*8 vec, hbar, mass, vecl, vecr, rwork, tfac
      complex*16 chamil, cu, ceig, cvecl, cvecr, work
      logical prn
      character*(*) mattyp
      character*80 title
      character*(*) coord, comp
      dimension p(n,n), dp(n,n), ddp(n,n)
      dimension q(n), wt(n), hamil(n,n), v(n)
      dimension u(n,n), eig(n), work(*), rwork(*) 
      dimension chamil(n,n), cu(n,n), ceig(n), vec(n,n) 
      dimension vecl(n,n), vecr(n,n), cvecl(n,n), cvecr(n,n)
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
      if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
         call czero(chamil,n*n)
      else
         call rzero(hamil,n*n)
      endif      
      tfac = hbar*hbar*.5d0/mass
      if(comp.eq.'x'.or.comp.eq.'y'.or.comp.eq.'z') then
         call tcart(p,dp,ddp,wt,hamil,chamil,ilft,irt,n,mattyp)
      elseif(comp.eq.'rho') then
c
c        for cylindrical rho coordinate
c
         call tcylin(p,dp,ddp,q,wt,hamil,chamil,ilft,irt,n,mattyp)
      elseif(comp.eq.'phi') then
c
c        for angular phi coordinate.  form is same as cartesian, just different
c        end points
c
         call tcart(p,dp,ddp,wt,hamil,chamil,ilft,irt,n,mattyp)
      elseif(comp.eq.'r') then
c 
c        radial coordinate
c
         call tcart(p,dp,ddp,wt,hamil,chamil,ilft,irt,n,mattyp)
      elseif(comp.eq.'theta') then
c        
c        spherical theta coordiante
c
         call tthet(p,dp,ddp,q,wt,hamil,chamil,work,ilft,irt,
     1              n,mattyp)
      else
         call lnkerr('error in coordinate type')
      endif
      if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
         call cvscal(chamil,chamil,tfac,n*n)
         do 10 i=1,n
            chamil(i,i) = chamil(i,i) + v(i)
 10      continue
         call cc2opy(chamil,cu,n*n)
         if (prn) then
             title='unperturbed hamiltonian'
             call prntcm(title,cu,n,n,n,n,iout)
         endif
         call rdiag(u,cu,eig,ceig,vec,cvecl,cvecr,
     1              work,lwork,rwork,n,n,mattyp)
      else
         call vscale(hamil,hamil,tfac,n*n)
         do 20 i=1,n
            hamil(i,i) = hamil(i,i) + v(i)
 20      continue
         call copy(hamil,u,n*n)
         if (prn) then
             title='unperturbed hamiltonian'
             call prntrm(title,u,n,n,n,n,iout)
         endif
         call rdiag(u,cu,eig,ceig,vec,cvecl,cvecr,
     1              work,lwork,rwork,n,n,mattyp)
      endif
      return
 1    format(/,5x,'forming a component of the unperturbed hamiltonian',
     1       /,5x,'coordinate system = ',a16,/,5x,'coordinate = ',a8)
      end       





