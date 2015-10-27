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
      subroutine h0(p,dp,ddp,q,wt,hamil,v,hbar,mass,n,fl,fr,
     1              qtyp,grd,prn)
      implicit integer (a-z)
      real*8 p, dp, ddp, q, wt, hamil, v
      real*8 hbar, mass, tfac
      real*8 y
      logical prn
      character*80 title
      character*(*) qtyp, grd
      dimension p(n,n), dp(n,n), ddp(n,n)
      dimension q(n), wt(n), hamil(n,n), v(n)
      data ldum/0/
      common/io/inp, iout
      pointer(py,y(1))
      u=1
      eig=u+n*n
      vec=eig+n
      sc=vec+n*n
      need=wpadti(sc+n*n)
      call memory(need,py,ngot,'diag',0)
      ilft=0
      irt=0
      if(fl.ne.0) then
         ilft=1
      endif
      if(fr.ne.0) then
         irt=1
      endif              
      write(iout,1) grd, qtyp
      call rzero(hamil,n*n)
      tfac = hbar*hbar*.5d0/mass
      if(qtyp.eq.'x'.or.qtyp.eq.'y'.or.qtyp.eq.'z') then
         call tcart(p,dp,ddp,wt,hamil,hamil,ilft,irt,n,'real-symmetric')
      elseif(qtyp.eq.'rho') then
c
c        for cylindrical rho coordinate
c
         call tcylin(p,dp,ddp,q,wt,hamil,hamil,ilft,irt,n,
     1               'real-symmetric')
      elseif(qtyp.eq.'phi') then
c
c        for angular phi coordinate.  form is same as cartesian, just different
c        end points
c
         call tcart(p,dp,ddp,wt,hamil,hamil,ilft,irt,n,'real-symmetric')
      elseif(qtyp.eq.'r') then
c 
c        radial coordinate
c
         call tcart(p,dp,ddp,wt,hamil,hamil,ilft,irt,n,'real-symmetric')
      elseif(qtyp.eq.'theta') then
c        
c        spherical theta coordinate
c
         call tthet(p,dp,ddp,q,wt,hamil,hamil,y(sc),ilft,irt,
     1              n,'real-symmetric')
      else
         call lnkerr('error in coordinate type')
      endif
      call vscale(hamil,hamil,tfac,n*n)
      call copy(hamil,y(u),n*n)
      call addd(y(u),y(u),v,n)
      if (prn) then
          title='unperturbed hamiltonian'
          call prntrm(title,y(u),n,n,n,n,iout)
      endif
      call rdiag(y(u),y(u),y(eig),y(eig),y(vec),y(vec),y(vec),
     1           y(sc),ldum,y(sc),n,n,'real-symmetric')
      title='eigenvalues'
      call prntrm(title,y(eig),n,1,n,1,iout)             
      call iosys('write real "ke eigenvalues-'//grd//' for '
     1           //qtyp//'" to bec',n,y(eig),0,' ')
      call iosys('write real "ke mtrx-'//grd//' for '//qtyp//
     1           '" to bec',n*n,hamil,0,' ')             
      call iosys('write real "trn mtrx-'//grd//' for '//qtyp//
     1           '" to bec',n*n,y(vec),0,' ')             
      call memory(-ngot,py,idum,'diag',idum)
      return
 1    format(/,5x,'forming a component of the one body hamiltonian',
     1       /,5x,'grid number = ',a3,/,5x,'coordinate = ',a8)
      end       





