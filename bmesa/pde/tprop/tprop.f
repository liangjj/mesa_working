*deck tprop.f
c***begin prologue     tprop
c***date written       990104   (yymmdd)
c***revision date               (yymmdd)
c***keywords           pde
c***                   
c***author             schneider, b. i.(nsf)
c***source             tprop
c***purpose            propagation of three dimensional, time-dependent
c***                   schroedinger equation.
c
c***description        the one, two or three dimensional, time-dependent
c***                   schroedinger equation is solved using a 
c***                   discrete variable representation in space and various
c***                   time-propagation schemes. 
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       tprop
      program tprop
c
      implicit integer (a-z)
      character*4096 ops
      character*8 prtflg
      character*2 itoc
      character*80 cpass, title, chrkey
      character*1600 card
      character*2 atom
      character*16 coord
      character*24 i0stat, tpert      
      character*8 qtyp, qdtyp
      logical dollar, logkey, prnply, prnwpt, prnh 
      logical prnht, prall, chko, check, nfix, prnh0
      logical toau, useau, phase, nospac
      real*8 z, amass, omegat, left, right, fpkey 
      real*8 pi, dleft, dright, t0, tf, delt
      real*8 hbar, eps
      real*8 massau, lenau, timau, dum, scale, sfac
      dimension nmax(3), npt(3)
      dimension lftbc(3), rtbc(3)
      dimension prnply(3), prnwpt(3), prnh(3), check(3)
      dimension q(3), wt(3), p(3), dp(3)
      dimension ddp(3), pn(3), dpn(3), ddpn(3)
      dimension eigc(3), wtc(3), work(6)
      dimension omegat(3), eig(3), left(3), right(3)
      dimension dleft(3), dright(3), eps(2)
      dimension ham(3), v(3), u(3), vec(3)
      dimension qtyp(3), qdtyp(3)
      dimension sfac(3)
      dimension nfix(2)
      pointer(pnt,z(1)),(pnt,ia(1))
      common/io/inp, iout      
      data pi/3.14159265358979323844d0/
c     hbar in joule-sec      
      data hbar/1.054592d-34/
c
c          mass in Kilograms length in meters
      data massau, lenau, timau / 9.109558d-31, 5.291771d-11, 
     1                            2.418884d-17 /
      call drum
      write(iout,*)
      write(iout,*) '    three dimensional time-propagation code       '
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
c
c                       set options
c
      numdim=intkey(ops,'number-of-dimensions',1,' ')   
      coord=chrkey(ops,'coordinate-system','cartesian',' ')
      prall=logkey(ops,'print=m7080=all',.false.,' ')
      nospac=logkey(ops,'no-space-hamiltonian',.false.,' ')
      write(iout,1) coord
      write(iout,2)
      if(prall) then
         prnply(1)=.true.
         prnply(2)=.true.
         prnply(3)=.true.
         prnwpt(1)=.true. 
         prnwpt(2)=.true.
         prnwpt(3)=.true.
         prnh0=.true. 
         prnht=.true. 
         prnh(1)=.true.
         prnh(2)=.true.
         prnh(3)=.true.
      else         
         prnply(1)=logkey(ops,'print=m7080=q1-polynomials',.false.,' ')
         prnply(2)=logkey(ops,'print=m7080=q2-polynomials',.false.,' ')
         prnply(3)=logkey(ops,'print=m7080=q3-polynomials',.false.,' ')
         prnwpt(1)=logkey(ops,'print=m7080=q1-points/weights',
     1                    .false.,' ')
         prnwpt(2)=logkey(ops,'print=m7080=q2-points/weights',
     1                    .false.,' ')
         prnwpt(3)=logkey(ops,'print=m7080=q3-points/weights',
     1                    .false.,' ')
         prnh0=logkey(ops,'print=m7080=h0',.false.,' ')
         prnht=logkey(ops,'print=m7080=h',.false.,' ')
         prnh(1)=logkey(ops,'print=m7080=q1-hamiltonian',.false.,' ')
         prnh(2)=logkey(ops,'print=m7080=q2-hamiltonian',.false.,' ')
         prnh(3)=logkey(ops,'print=m7080=q3-hamiltonian',.false.,' ')
      endif
      chko=logkey(ops,'m7080=check-orthogonality',.false.,' ')
      if(chko) then
         check(1)=.true.
         check(2)=.true.
         check(3)=.true.
      else                  
         check(1)=logkey(ops,'check-q1-orthogonality',.false.,' ')
         check(2)=logkey(ops,'check-q2-orthogonality',.false.,' ')
         check(3)=logkey(ops,'check-q3-orthogonality',.false.,' ')
      endif
      toau=logkey(ops,'to-atomic-units',.false.,' ')
      useau=logkey(ops,'use-atomic-units',.false.,' ')
      if(numdim.eq.1) then
         write(iout,3) prnply(1),
     1                 prnwpt(1), prnh0,
     2                 prnh(1), 
     3                 check(1)
      endif
      if(numdim.eq.2) then
         write(iout,4) prnply(1), prnply(2),
     1                 prnwpt(1), prnwpt(2), prnh0,
     2                 prnh(1), prnh(2), 
     3                 check(1), check(2)
      endif     
      if(numdim.eq.3) then
         write(iout,5) prnply(1), prnply(2), prnply(3),
     1                 prnwpt(1), prnwpt(2), prnwpt(3), prnh0,
     2                 prnh(1), prnh(2), prnh(3), 
     3                 check(1), check(2), check(3)
      endif                     
c
c              set various program options
c      
c
c         set trap configuration, numbers of atoms, various constants
c         and other physical parameters
c       
      if ( dollar('$trap',card,cpass,inp) ) then
           atom=chrkey(card,'atom','Cs',' ')
           call atmdat(atom,amass,omegat)
           if(numdim.eq.1) then
              write(iout,6) amass, omegat(1)
           elseif(numdim.eq.2) then
              write(iout,7) amass, (omegat(i),i=1,2)
           elseif(numdim.eq.3) then
              write(iout,8) amass, (omegat(i),i=1,3)
           endif
      endif
c
c               spatial basis set information
c       
      if ( dollar('$q1-polynomials',card,cpass,inp) ) then
           call bfdat(.false.,card,qtyp(1),qdtyp(1),left(1),right(1),
     1                dleft(1),dright(1),lftbc(1),rtbc(1),
     2                npt(1),nmax(1),nfix,1)
      endif     
      if(numdim.gt.1) then
         if ( dollar('$q2-polynomials',card,cpass,inp) ) then
              call bfdat(.false.,card,qtyp(2),qdtyp(2),left(2),right(2),
     1                   dleft(2),dright(2),lftbc(2),rtbc(2),
     2                   npt(2),nmax(2),nfix,1)
         endif
      endif
      if(numdim.gt.2) then
         if ( dollar('$q3-polynomials',card,cpass,inp) ) then
              call bfdat(.false.,card,qtyp(3),qdtyp(3),left(3),right(3),
     1                   dleft(3),dright(3),lftbc(3),rtbc(3),
     2                   npt(3),nmax(3),nfix,1)
         endif
      endif
      if( dollar('$time',card,cpass,inp) ) then
          i0stat=chrkey(card,'driver','state-vector',' ')
          tpert=chrkey(card,'time-perturbation','none',' ')
          state=intkey(card,'initial-state',0,' ')
          phase=logkey(card,'cosine-phase',.false.,' ')
          t0=fpkey(card,'initial-time',0.d0,' ')
          tf=fpkey(card,'final-time',5.d0,' ')
          delt=fpkey(card,'time-step',.01d0,' ')
          eps(1)=1.d-08
          eps(2)=1.e-08
          call fparr(card,'accuracy',eps,2,' ')      
          write(iout,9) t0, tf, delt, eps(1), eps(2), tpert,
     1                  i0stat, state
      endif      
      if(toau) then
c
c        converts everything to atomic units
c
         write(iout,*) 'converting to atomic units'
         hbar=1.d0
         scatl=scatl/lenau
         do 10 i=1,numdim
            left(i)=left(i)/lenau
            right(i)=right(i)/lenau
            omegat(i)=omegat(i)*timau
 10      continue   
         amass=amass/massau        
      endif
      if(useau) then
c
c        if this option is used then all is entered in atomic units
c
         write(iout,*) 'assuming atomic units'
         hbar=1.d0
         amass=1.d0
         do 20 i=1,numdim
            omegat(i)=1.d0
 20      continue            
      endif
      maxd=0
      do 30 i=1,numdim
         maxd=max(maxd,npt(i))
 30   continue
      n3d=1
      do 40 i=1,numdim
         n3d=n3d*nmax(i)
 40   continue
      neq=2*n3d 
      write(iout,11) n3d
      ioff=1
      call wrdadd(q(1),wt(1),eigc(1),wtc(1),
     1            p(1),dp(1),ddp(1),
     2            pn(1),dpn(1),ddpn(1),
     3            ham(1),eig(1),v(1),u(1),vec(1),
     4            npt(1),ioff)
      if(numdim.gt.1) then
         call wrdadd(q(2),wt(2),eigc(2),wtc(2),
     1               p(2),dp(2),ddp(2),
     2               pn(2),dpn(2),ddpn(2),
     3               ham(2),eig(2),v(2),u(2),vec(2),
     4               npt(2),ioff)
      endif
      if(numdim.gt.2) then
         call wrdadd(q(3),wt(3),eigc(3),wtc(3),
     1               p(3),dp(3),ddp(3),
     2               pn(3),dpn(3),ddpn(3),
     3               ham(3),eig(3),v(3),u(3),vec(3),
     4               npt(3),ioff)
      endif
      ind=wpadti(ioff)
      psi0=iadtwp(ind+2*(numdim+1)*n3d)
      vtot=psi0+neq
      work(1)=vtot+n3d
      work(2)=work(1)+max(2*maxd*maxd,3*n3d)
      work(3)=work(2)+max(2*maxd*maxd,n3d)
      work(4)=work(3)+max(2*maxd*maxd,n3d)
      need=wpadti(work(4)+2*n3d)
      call memory(need,pnt,ngot,'tprop',0)
c
c     in order to proceed it is necessary to be able to calculate integrals
c     between the orthogonal polynomials to be generated on the desired
c     integration interval.  the calculation of the points and weights is
c     done on the interval ( -1.,1. ) and then converted to ( left, right ).
c     the fixed points on the standard interval are the two end points.      
c
c     once that is done the time-independent parts of the hamiltonian
c     are generated.
c
      do 50 i=1,numdim
         call getqpt(z(q(i)),z(wt(i)),left(i),right(i),
     1                 qdtyp(i),'before',z(work(1)),
     2                 nfix,npt(i),npt(i),1,.false.)
         if(qdtyp(i).eq.'hermite') then
            call smul(z(q(i)),z(q(i)),sfac(i),npt(i))
            call smul(z(wt(i)),z(wt(i)),sfac(i),npt(i))
         endif   
         call lgrply(z(p(i)),z(dp(i)),z(ddp(i)),z(q(i)),
     1               z(q(i)),npt(i),npt(i),.false.)
         call prepfn(z(p(i)),z(dp(i)),z(ddp(i)),z(q(i)),
     1               z(wt(i)),npt(i),sfac(i),qdtyp(i))
         call copy(z(p(i)),z(pn(i)),npt(i)*npt(i))
         call copy(z(dp(i)),z(dpn(i)),npt(i)*npt(i))
         call copy(z(ddp(i)),z(ddpn(i)),npt(i)*npt(i))
         if(prnply(i)) then                  
            title='dvr polynomials '//itoc(i)//' dimension'
            call prntfm(title,z(pn(i)),npt(i),npt(i),
     1                  npt(i),npt(i),iout)
            title='first derivative of dvr polynomials '
     1             //itoc(i)//' dimension'
            call prntfm(title,z(dpn(i)),npt(i),npt(i),
     1                  npt(i),npt(i),iout)
            title='second derivative of dvr polynomials '
     1             //itoc(i)//' dimension'
            call prntfm(title,z(ddpn(i)),npt(i),npt(i),
     1                  npt(i),npt(i),iout)
         endif
         if(check(i)) then
            call chk(z(pn(i)),z(wt(i)),z(work(1)),z(work(2)),
     1               npt(i),npt(i))
            title='overlap matrix dvr polynomials '//itoc(i)
     1             //' dimension'
            call prntfm(title,z(work(2)),npt(i),npt(i),
     1                  npt(i),npt(i),iout)
         endif
         if(lftbc(i).eq.0.or.rtbc(i).eq.0) then
            call dropfn(z(p(i)),z(dp(i)),z(ddp(i)),
     1                  z(pn(i)),z(dpn(i)),
     2                  z(ddpn(i)),z(q(i)),
     3                  z(wt(i)),z(eigc(i)),z(wtc(i)),
     4                  lftbc(i),rtbc(i),npt(i),nmax(i))
         endif     
         if(prnply(i)) then                  
            title='dvr polynomials '//itoc(i)//' dimension'
            call prntfm(title,z(pn(i)),nmax(i),nmax(i),
     1                  nmax(i),nmax(i),iout)
            title='first derivative of dvr polynomials '
     1             //itoc(i)//' dimension'
            call prntfm(title,z(dpn(i)),nmax(i),nmax(i),
     1                  nmax(i),nmax(i),iout)
            title='second derivative of dvr polynomials '
     1             //itoc(i)//' dimension'
            call prntfm(title,z(ddpn(i)),nmax(i),nmax(i),
     1                  nmax(i),nmax(i),iout)
         endif
         if(check(i)) then
            call chk(z(pn(i)),z(wtc(i)),z(work(1)),z(work(2)),
     1               nmax(i),nmax(i),)
            title='overlap matrix dvr polynomials '//itoc(i)
     1             //' dimension'
            call prntfm(title,z(work(2)),nmax(i),nmax(i),
     1                  nmax(i),nmax(i),iout)
         endif
         call vmat1d(z(eigc(i)),z(v(i)),hbar,amass,scale,nmax(i),
     1               i,prnh0,ops)
c
c        the hamiltonian matrix contains the kinetic energy and any
c        one body potential.
c
         call h0(z(pn(i)),z(dpn(i)),z(ddpn(i)),z(eigc(i)),
     1           z(wtc(i)),z(ham(i)),z(v(i)),z(u(i)),z(eig(i)),
     2           z(vec(i)),z(work(1)),hbar,amass,nmax(i),
     3           lftbc(i),rtbc(i),coord,qtyp(i),prnh0)
         call vscale(z(eig(i)),z(eig(i)),scale,nmax(i))
         if(prnh0) then
            title='eigenvalues of unperturbed hamiltonian '//itoc(i)
     1             //' dimension'
            call prntfm(title,z(eig(i)),nmax(i),1,nmax(i),1,iout)
            call vscale(z(eig(i)),z(eig(i)),1.d0/scale,nmax(i))
         endif
         if(nospac) then
            call rzero(z(ham(i)),nmax(i)*nmax(i))
            call rzero(z(eig(i)),nmax(i))
         endif            
 50   continue
c
c     calculate the initial state wavefunction
c
      call setind(ia(ind),nmax,n3d,numdim)
      if(i0stat.eq.'state-vector') then
c
c        calculate the projection on to the dvr basis
c        for an initial state which is a pure stationary state 
c        of the unperturbed hamiltonian times an 
c        oscillatory factor.
c
         call srteig(z(vec(1)),z(vec(2)),z(vec(3)),
     1               z(eig(1)),z(eig(2)),z(eig(3)),z(work(1)),
     2               z(psi0),t0,ia(ind),nmax,numdim,n3d,
     3               state+1,.false.)
      elseif(i0stat.eq.'wavepacket') then         
c
c        calculate the projection on to the dvr basis
c        for an initial state which is a wavepacket.
c
c        first we calculate the projection onto the pure time part.  then
c        we do the spatial projection
         call gpaket(z(pn(1)),z(pn(2)),z(pn(3)),
     1               z(eigc(1)),z(eigc(2)),z(eigc(3)),
     2               z(wtc(1)),z(wtc(2)),z(wtc(3)),
     3               z(psi0),z(work(1)),z(work(2)),z(work(3)),
     4               z(work(4)),n3d,nmax(1),nmax(2),nmax(3),
     5               numdim,coord,phase)
      elseif(i0stat.eq.'one') then
             call vfill(z(psi0),1.d0,nmax(1))
      endif
c     
c             now that we have the initial state we can begin
c             the time propagation.
c
c             the wavefunction and hence the hamiltonian is written
c             in the dvr representation.  the time evolution of the
c             system is contained in the coefficients of the
c             expansion of the wavefunction in the basis.
c 
      call tdrver(z(pn(1)),z(pn(2)),z(pn(3)),
     1            z(ham(1)),z(ham(2)),z(ham(3)),
     2            z(eigc(1)),z(eigc(2)),z(eigc(3)),
     3            z(vtot),z(psi0),t0,tf,delt,eps,tpert,
     4            hbar,numdim,nmax,n3d)
      call memory(-ngot,pnt,idum,'tprop',idum)
      call chainx(0)               
      stop
 1    format(/,5x,'coordinate system    = ',a16)
 2    format(/,15x,'program options',/,/,5x,
     1             'diagonalize zeroth-order hamiltonian')  
 3    format(/,15x,'code options',/,
     1       /,5x,'print q1 polynomials                   = ',l1,
     2       /,5x,'print q1 points/weights                = ',l1,
     3       /,5x,'print spatial hamiltonian information  = ',l1, 
     4       /,5x,'print q1 eigenvectors                  = ',l1, 
     5       /,5x,'check q1 orthonormality                = ',l1)
 4    format(/,15x,'code options',/,
     1       /,5x,'print q1 polynomials                   = ',l1,
     2       /,5x,'print q2 polynomials                   = ',l1,
     3       /,5x,'print q1 points/weights                = ',l1,
     4       /,5x,'print q2 points/weights                = ',l1,
     5       /,5x,'print spatial hamiltonian information  = ',l1, 
     6       /,5x,'print q1 eigenvectors                  = ',l1, 
     7       /,5x,'print q2 eigenvectors                  = ',l1, 
     8       /,5x,'check q1 orthonormality                = ',l1, 
     9       /,5x,'check q2 orthonormality                = ',l1) 
 5    format(/,15x,'code options',/,
     1       /,5x,'print q1 polynomials                   = ',l1,
     2       /,5x,'print q2 polynomials                   = ',l1,
     3       /,5x,'print q3 polynomials                   = ',l1,
     4       /,5x,'print q1 points/weights                = ',l1,
     5       /,5x,'print q2 points/weights                = ',l1,
     6       /,5x,'print q3 points/weights                = ',l1,
     7       /,5x,'print spatial hamiltonian information  = ',l1, 
     8       /,5x,'print q1 eigenvectors                  = ',l1, 
     9       /,5x,'print q2 eigenvectors                  = ',l1, 
     x       /,5x,'print q3 eigenvectors                  = ',l1, 
     x       /,5x,'check q1 orthonormality                = ',l1, 
     x       /,5x,'check q2 orthonormality                = ',l1, 
     x       /,5x,'check q3 orthonormality                = ',l1) 
 6    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                      = ',e15.8,/,5x,
     2             'trap q1 frequency                = ',e15.8)
 7    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                       = ',e15.8,/,5x,
     2             'trap q1 frequency                 = ',e15.8,/,5x,
     3             'trap q2 frequency                 = ',e15.8)
 8    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                       = ',e15.8,/,5x,
     2             'trap q1 frequency                 = ',e15.8,/,5x,
     3             'trap q2 frequency                 = ',e15.8,/,5x,
     4             'trap q3 frequency                 = ',e15.8)
 9    format(/,5x,'time-dependent data',/,5x,
     1            'initial time             = ',e15.8,/,5x,
     2            'final time               = ',e15.8,/,5x,     
     3            'starting step size       = ',e15.8,/,5x,
     4            'relative error           = ',e15.8,/,5x,
     5            'absolute error           = ',e15.8,/,5x,
     6            'space-time perturbation  = ',a24,/,5x,
     7            'wavepacket type          = ',a24,/,5x,
     8            'initial-state            = ',i4)
 11   format(/,5x,'basis set size = ',i6)
      end
