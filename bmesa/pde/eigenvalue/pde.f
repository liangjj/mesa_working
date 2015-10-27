*deck pde.f 
c***begin prologue     pde
c***date written       970810   (yymmdd)
c***revision date               (yymmdd)
c***keywords           pde
c***                   
c***author             schneider, b. i.(nsf)
c***source             pde
c***purpose            three dimensional eigenvalue code
c
c***description        the one, two or three dimensional eigenvalue
c***                   problem is solved using a discrete variable
c***                   representation. 
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       pde
      program pde
c
      implicit integer (a-z)
      parameter ( ngmax=10 )
      character*4096 ops
      character*8 prtflg
      character*2 itoc
      character*80 cpass, title, chrkey
      character*1600 card
      character*128 filham
      character*2 atom
      character*24 mattyp
      character*16 coord
      character*8 qtyp, qdtyp
      logical dollar, logkey, prnply, prnwpt, prnh 
      logical prnht, prall, chko, check, nfix, prnh0
      logical toau, useau, hamd, incore, cgrid
      logical itdiag, prbufh
      logical prdvd, dvdall, frmdsk
      logical bypass, scatt, typot
      real*8 z, amass, omegat, left, right, energy 
      real*8 pi, dleft, dright
      real*8 thresh, cnverg, hbar
      real*8 massau, lenau, timau, dum, scale, sfac
      complex*16 cdum
      pointer(pnt,z(1)),(pnt,ia(1))
      dimension nmax(ngmax,3), npt(ngmax,3)
      dimension lftbc(3), rtbc(3), dim(3)
      dimension prnply(3), prnwpt(3), prnh(3), check(3)
      dimension q(ngmax,3), wt(ngmax,3), p(ngmax,3), dp(ngmax,3)
      dimension ddp(ngmax,3), pn(ngmax,3), dpn(ngmax,3), ddpn(ngmax,3)
      dimension eigc(ngmax,3), wtc(ngmax,3), work(4)
      dimension omegat(3), eig(3), left(3), right(3)
      dimension dleft(3), dright(3)
      dimension ham(3), v(3), u(3), vecl(3), vecr(3)
      dimension prdvd(11), qtyp(3), qdtyp(3)
      dimension sfac(3), energy(100)
      dimension nfix(2), tmp(5)
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
      write(iout,*) '    three dimensional pde eigenvalue code       '
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
c
c                       set options
c
      frmdsk=logkey(ops,'from-disk',.false.,' ')   
      ngrid=intkey(ops,'number-of-grids',1,' ')
      ug=intkey(ops,'grid',1,' ')
      numdim=intkey(ops,'number-of-dimensions',1,' ')   
      coord=chrkey(ops,'coordinate-system','cartesian',' ')
      scatt=logkey(ops,'scattering-calculation',.false.,' ')
      typot=logkey(ops,'coulomb-potential',.false.,' ')
      itdiag=logkey(ops,'iterative-diagonalization',.false.,' ')      
      mattyp=chrkey(ops,'matrix-type','real-symmetric',' ')
      prall=logkey(ops,'print=m7020=all',.false.,' ')
      bypass=logkey(ops,'bypass=on',.false.,' ')
c----------------------------------------------------------------------c
c   1. if task = diagonalize-h0, the code will diagonalize a hamiltonian
c      consisting of the kinetic energy plus one and two body local 
c      interactions.  the interaction terms are not assumed to depend on
c      the wavefunction.
c----------------------------------------------------------------------c 
      write(iout,1) coord
      write(iout,2)
      if(bypass) then
         call genmat(numdim,itdiag,mattyp)
         call chainx(0)               
         stop
      endif   
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
         prnply(1)=logkey(ops,'print=m7020=q1-polynomials',.false.,' ')
         prnply(2)=logkey(ops,'print=m7020=q2-polynomials',.false.,' ')
         prnply(3)=logkey(ops,'print=m7020=q3-polynomials',.false.,' ')
         prnwpt(1)=logkey(ops,'print=m7020=q1-points/weights',
     1                    .false.,' ')
         prnwpt(2)=logkey(ops,'print=m7020=q2-points/weights',
     1                    .false.,' ')
         prnwpt(3)=logkey(ops,'print=m7020=q3-points/weights',
     1                    .false.,' ')
         prnh0=logkey(ops,'print=m7020=h0',.false.,' ')
         prnht=logkey(ops,'print=m7020=h',.false.,' ')
         prnh(1)=logkey(ops,'print=m7020=q1-hamiltonian',.false.,' ')
         prnh(2)=logkey(ops,'print=m7020=q2-hamiltonian',.false.,' ')
         prnh(3)=logkey(ops,'print=m7020=q3-hamiltonian',.false.,' ')
      endif
      chko=logkey(ops,'m7020=check-orthogonality',.false.,' ')
      if(chko) then
         check(1)=.true.
         check(2)=.true.
         check(3)=.true.
      else                  
         check(1)=logkey(ops,'check-q1-orthogonality',.false.,' ')
         check(2)=logkey(ops,'check-q2-orthogonality',.false.,' ')
         check(3)=logkey(ops,'check-q3-orthogonality',.false.,' ')
      endif
      if(scatt) then
         nen=intkey(ops,'number-of-energies',1,' ')
         call fparr(ops,'energies',energy,nen,' ')
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
           call bfdat(frmdsk,card,qtyp(1),qdtyp(1),left(1),right(1),
     1                dleft(1),dright(1),lftbc(1),rtbc(1),
     2                npt(1,1),nmax(1,1),nfix,ngrid)
      endif     
      if(numdim.gt.1) then
         if ( dollar('$q2-polynomials',card,cpass,inp) ) then
              call bfdat(frmdsk,card,qtyp(2),qdtyp(2),left(2),right(2),
     1                   dleft(2),dright(2),lftbc(2),rtbc(2),
     2                   npt(1,2),nmax(1,2),nfix,ngrid)
         endif
      endif
      if(numdim.gt.2) then
         if ( dollar('$q3-polynomials',card,cpass,inp) ) then
              call bfdat(frmdsk,card,qtyp(3),qdtyp(3),left(3),right(3),
     1                   dleft(3),dright(3),lftbc(3),rtbc(3),
     2                   npt(1,3),nmax(1,3),nfix,ngrid)
         endif
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
      do 30 i=1,numdim
         dim(i)=0
         do 40 j=1,ngrid
            dim(i)=max(dim(i),npt(j,i))            
 40      continue
 30   continue
      maxd=max(dim(1),dim(2),dim(3))
      n3d=1
      do 50 i=1,numdim
         n3d=n3d*nmax(ug,i)
 50   continue 
      nroots=n3d         
c
c             set diagonalization procedures
c      
      if(itdiag) then
         call dvddat(card,cpass,n3d,nroots,ntrials,nattim,cnverg,
     1               thresh,niter,nvec,lenbuf,cgrid,prbufh,hamd,
     2               n0,prdvd,dvdall,filham)
         write(iout,9) nroots, nattim, thresh, cnverg, niter
      else
         nroots=intkey(ops,'number-of-roots',n3d,' ')
         nroots=min(nroots,n3d)
         write(iout,11) n3d, nroots
      endif
      ioff=1
c
c     note the alignment of variables.  it is important in passing these 
c     arrays later.
c 
      call wrdcnt(q(ug,1),wt(ug,1),eigc(ug,1),wtc(ug,1),
     1            p(ug,1),dp(ug,1),ddp(ug,1),
     2            pn(ug,1),dpn(ug,1),ddpn(ug,1),
     3            ham(1),eig(1),v(1),u(1),vecl(1),vecr(1),
     4            mattyp,dim(1),ngrid,ioff)
      if(numdim.gt.1) then
         call wrdcnt(q(ug,2),wt(ug,2),eigc(ug,2),wtc(ug,2),
     1               p(ug,2),dp(ug,2),ddp(ug,2),
     2               pn(ug,2),dpn(ug,2),ddpn(ug,2),
     3               ham(2),eig(2),v(2),u(2),vecl(2),vecr(2),
     4               mattyp,dim(2),ngrid,ioff)
      endif
      if(numdim.gt.2) then
         call wrdcnt(q(ug,3),wt(ug,3),eigc(ug,3),wtc(ug,3),
     1               p(ug,3),dp(ug,3),ddp(ug,3),
     2               pn(ug,3),dpn(ug,3),ddpn(ug,3),
     3               ham(3),eig(3),v(3),u(3),vecl(3),vecr(3),
     4               mattyp,dim(3),ngrid,ioff)
      endif
      ind=wpadti(ioff)
      idum=ind+(numdim+1)*n3d
      words=iadtwp(idum+numdim*n3d)            
      space=words    
      work(1)=words
      if(mattyp.eq.'complex'.or.
     1          mattyp.eq.'real-unsymmetric') then
         lwork=20*maxd
         work(2)=work(1)+max(4*maxd*maxd,6*n3d,lwork)
         words=work(2)+4*maxd*maxd
      elseif(mattyp.eq.'real-symmetric') then
         lwork=0
         work(2)=work(1)+max(2*maxd*maxd,3*n3d)
         words=work(2)+2*maxd*maxd
      endif                  
      if(.not.itdiag) then
         eigtot=words
         if(mattyp.eq.'complex'.or.
     1             mattyp.eq.'real-unsymmetric') then
            hamtot=eigtot+2*n3d
            vtot=hamtot+2*n3d*n3d
            bigvec=vtot+2*n3d
            bigvec2=bigvec+2*n3d*n3d
            words=bigvec2+2*n3d*n3d
         elseif(mattyp.eq.'real-symmetric') then
            hamtot=eigtot+n3d
            vtot=hamtot+n3d*n3d
            bigvec=vtot+n3d
            bigvec2=bigvec+n3d*n3d
            words=bigvec2+n3d*n3d
         endif     
      else
         if(mattyp.eq.'complex'.or.
     1             mattyp.eq.'real-unsymmetric') then
             hbuf=words
             ibuf=wpadti(hbuf+2*lenbuf)
             diag=iadtwp(ibuf+2*lenbuf)
             vtot=diag+n3d*2
             etrial=vtot+n3d*2
             trials=etrial+n3d*2
             diagv=trials+n3d*ntrials*2
             psi=diagv+n3d*2
             pvec=psi+n3d*nroots*2
             hpvec=pvec+n3d*nvec*2
             vec=hpvec+nvec*n3d*2
             bmat=vec+n3d*nvec*2
             bmatm=bmat+nvec*nvec*2
             svec=bmatm+nvec*nvec*2
             sveca=svec+nvec*nvec*2
             eigtot=sveca+nvec*nvec*2
             etmp=eigtot+n3d*2
             work(3)=etmp+n3d*2
             mwork=10*nvec*2
             work(4)=work(3)+max(n3d*2*nvec,mwork)
             words=work(4)+max(n3d*2*nvec,4*nvec)
         elseif(mattyp.eq.'real-symmetric') then
             hbuf=words
             ibuf=wpadti(hbuf+lenbuf)
             diag=iadtwp(ibuf+2*lenbuf)
             vtot=diag+n3d
             etrial=vtot+n3d
             trials=etrial+n3d
             diagv=trials+n3d*ntrials
             psi=diagv+n3d
             pvec=psi+n3d*nroots
             hpvec=pvec+n3d*nvec
             vec=hpvec+nvec*n3d
             bmat=vec+n3d*nvec
             bmatm=bmat+nvec*nvec
             svec=bmatm+nvec*nvec
             sveca=svec+nvec*nvec
             eigtot=sveca+nvec*nvec
             etmp=eigtot+n3d
             work(3)=etmp+n3d
             mwork=0
             work(4)=work(3)+n3d*nvec
             words=work(4)+max(n3d*nvec,2*nvec)
         endif
      endif      
      words=wpadti(words)
      call manmem(0,idum,idum,'pde',idum)
      call manmem(words,pnt,ngot,'pde',0)
      lwork=max(lwork,mwork)
c
c     in order to proceed it is necessary to be able to calculate integrals
c     between the orthogonal polynomials to be generated on the desired
c     integration interval.  the calculation of the points and weights is
c     done on the interval ( -1.,1. ) and then converted to ( left, right ).
c     the fixed points on the standard interval are the two end points.      
c
      do 80 i=1,numdim
c      
         write(iout,12) i, ug
         if(.not.frmdsk) then
c             call getqpt(z(q(ug,i)),z(wt(ug,i)),left(i),right(i),
c     1                   qdtyp(i),'before',z(work(1)),
c     2                   nfix,npt(ug,i),npt(ug,i),1,.false.)
             call getqpt(z(q(ug,i)),z(wt(ug,i)),left(i),right(i),
     1                   qdtyp(i),'before',z(work(1)),
     2                   nfix,npt(ug,i),npt(ug,i),1,.false.)
             if(qdtyp(i).eq.'hermite') then
                call smul(z(q(ug,i)),z(q(ug,i)),sfac(i),npt(ug,i))
                call smul(z(wt(ug,i)),z(wt(ug,i)),sfac(i),npt(ug,i))
             endif   
             call lgrply(z(p(ug,i)),z(dp(ug,i)),z(ddp(ug,i)),z(q(ug,i)),
     1                   z(q(ug,i)),npt(ug,i),npt(ug,i),.false.)
             call prepfn(z(p(ug,i)),z(dp(ug,i)),z(ddp(ug,i)),z(q(ug,i)),
     1                   z(wt(ug,i)),npt(ug,i),sfac(i),qdtyp(i))
             call copy(z(p(ug,i)),z(pn(ug,i)),npt(ug,i)*npt(ug,i))
             call copy(z(dp(ug,i)),z(dpn(ug,i)),npt(ug,i)*npt(ug,i))
             call copy(z(ddp(ug,i)),z(ddpn(ug,i)),npt(ug,i)*npt(ug,i))
             if(prnply(i)) then                  
                title='dvr polynomials '//itoc(i)//' dimension'
                call prntfm(title,z(pn(ug,i)),npt(ug,i),npt(ug,i),
     1                      npt(ug,i),npt(ug,i),iout)
                title='first derivative of dvr polynomials '
     1                 //itoc(i)//' dimension'
                call prntfm(title,z(dpn(ug,i)),npt(ug,i),npt(ug,i),
     1                      npt(ug,i),npt(ug,i),iout)
                title='second derivative of dvr polynomials '
     1                 //itoc(i)//' dimension'
                call prntfm(title,z(ddpn(ug,i)),npt(ug,i),npt(ug,i),
     1                      npt(ug,i),npt(ug,i),iout)
             endif
             if(check(i)) then
                call chk(z(pn(ug,i)),z(wt(ug,i)),z(work(1)),z(work(2)),
     1                   npt(ug,i),npt(ug,i))
                title='overlap matrix dvr polynomials '//itoc(i)
     1                 //' dimension'
                call prntfm(title,z(work(2)),npt(ug,i),npt(ug,i),
     1                      npt(ug,i),npt(ug,i),iout)
             endif
             if(lftbc(i).eq.0.or.rtbc(i).eq.0) then
                call dropfn(z(p(ug,i)),z(dp(ug,i)),z(ddp(ug,i)),
     1                      z(pn(ug,i)),z(dpn(ug,i)),
     2                      z(ddpn(ug,i)),z(q(ug,i)),
     3                      z(wt(ug,i)),z(eigc(ug,i)),z(wtc(ug,i)),
     4                      lftbc(i),rtbc(i),npt(ug,i),nmax(ug,i))
             endif     
         else
             call iosys('read real "points grid = '//itoc(ug)
     1                  //'" from rwf',nmax(ug,i),z(eigc(ug,i)),0,' ')
             call iosys('read real "weights grid = '//itoc(ug)
     1                 //'" from rwf',nmax(ug,i),z(wtc(ug,i)),0,' ')
             call iosys('read real "lgrbf = '//itoc(ug)//
     1                  ' grid = '//itoc(ug)//'" from rwf',
     2                   nmax(ug,i)*nmax(ug,i),z(pn(ug,i)),0,' ')
             call iosys('read real "dlgrbf = '//itoc(ug)//
     1                 ' grid = '//itoc(ug)//'" from rwf',
     2                   nmax(ug,i)*nmax(ug,i),z(dpn(ug,i)),0,' ')
             call iosys('read real "ddlgrbf = '//itoc(ug)//
     3                  ' grid = '//itoc(ug)//'" from rwf',
     4                   nmax(ug,i)*nmax(ug,i),z(ddpn(ug,i)),0,' ')
         endif              
         if(prnply(i)) then                  
            title='dvr polynomials '//itoc(i)//' dimension'
            call prntfm(title,z(pn(ug,i)),nmax(ug,i),nmax(ug,i),
     1                  nmax(ug,i),nmax(ug,i),iout)
            title='first derivative of dvr polynomials '
     1             //itoc(i)//' dimension'
            call prntfm(title,z(dpn(ug,i)),nmax(ug,i),nmax(ug,i),
     1                  nmax(ug,i),nmax(ug,i),iout)
            title='second derivative of dvr polynomials '
     1             //itoc(i)//' dimension'
            call prntfm(title,z(ddpn(ug,i)),nmax(ug,i),nmax(ug,i),
     1                  nmax(ug,i),nmax(ug,i),iout)
         endif
         if(check(i)) then
            call chk(z(pn(ug,i)),z(wtc(ug,i)),z(work(1)),z(work(2)),
     1               nmax(ug,i),nmax(ug,i),)
            title='overlap matrix dvr polynomials '//itoc(i)
     1             //' dimension'
            call prntfm(title,z(work(2)),nmax(ug,i),nmax(ug,i),
     1                  nmax(ug,i),nmax(ug,i),iout)
         endif
         call vmat1d(z(eigc(ug,i)),z(v(i)),hbar,amass,scale,nmax(ug,i),
     1               i,prnh0,ops)
c
c        the hamiltonian matrix contains the kinetic energy and any
c        one body potential.
c
         call ham0(z(pn(ug,i)),z(dpn(ug,i)),z(ddpn(ug,i)),z(eigc(ug,i)),
     1             z(wtc(ug,i)),z(ham(i)),z(ham(i)),z(v(i)),
     2             z(u(i)),z(u(i)),z(eig(i)),z(eig(i)),
     3             z(vecl(i)),z(vecl(i)),z(vecr(i)),
     4             z(vecl(i)),z(vecr(i)),z(work(1)),lwork,z(work(2)),
     5             hbar,amass,nmax(ug,i),lftbc(i),rtbc(i),coord,
     6             qtyp(i),prnh0,mattyp)
         if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
            call cvscal(z(eig(i)),z(eig(i)),scale,nmax(ug,i))
c            if(prnh0) then
            title='eigenvalues of unperturbed hamiltonian '//itoc(i)
     1             //' dimension'
            call prntcm(title,z(eig(i)),nmax(ug,i),1,nmax(ug,i),
     1                  1,iout)
c            endif
            call cvscal(z(eig(i)),z(eig(i)),1.d0/scale,nmax(ug,i))
         elseif(mattyp.eq.'real-symmetric') then
            call vscale(z(eig(i)),z(eig(i)),scale,nmax(ug,i))
c           if(prnh0) then
            title='eigenvalues of unperturbed hamiltonian '//itoc(i)
     1             //' dimension'
            call prntfm(title,z(eig(i)),nmax(ug,i),1,
     1                  nmax(ug,i),1,iout)
            call vscale(z(eig(i)),z(eig(i)),1.d0/scale,nmax(ug,i))
c            endif            
         endif
 80   continue
      if(.not.itdiag) then 
         nroots=n3d
      endif          
      call setind(ia(ind),nmax,n3d,ug,numdim,ngmax)
      if(numdim.eq.1) then
         call vmat(z(eigc(ug,1)),dum,dum,z(vtot),z(work(1)),
     1             hbar,amass,scale,n3d,nmax(ug,1),nmax(ug,2),
     2             nmax(ug,3),numdim,prnh0,ops)   
      elseif(numdim.eq.2) then
         call vmat(z(eigc(ug,1)),z(eigc(ug,2)),dum,z(vtot),z(work(1)),
     1             hbar,amass,scale,n3d,nmax(ug,1),nmax(ug,2),
     2             nmax(ug,3),numdim,prnh0,ops)   
      elseif(numdim.eq.3) then
         call vmat(z(eigc(ug,1)),z(eigc(ug,2)),z(eigc(ug,3)),z(vtot),
     1             z(work(1)),hbar,amass,scale,n3d,
     2             nmax(ug,1),nmax(ug,2),nmax(ug,3),
     3             numdim,prnh0,ops)     
      endif
c
c     subtract the one body interactions included in the zeroth
c     order hamiltonian.
c
      call vpert(z(vtot),z(v(1)),z(v(2)),z(v(3)),
     1           n3d,nmax(ug,1),nmax(ug,2),nmax(ug,3),numdim,mattyp)
      if(itdiag) then
         call iosys('write integer "length of davidson vector" '//
     1              'to ham',1,n3d,0,' ')
         if(numdim.eq.1) then
            call vtrial(z(ham(1)),z(ham(1)),dum,cdum,dum,cdum,
     1                  z(vecl(1)),z(vecl(1)),dum,cdum,dum,cdum,
     2                  z(eig(1)),z(eig(1)),dum,cdum,dum,cdum,
     3                  z(trials),z(trials),z(etrial),z(etrial),
     4                  ia(idum),nmax(ug,1),nmax(ug,2),nmax(ug,3),
     5                  numdim,n3d,ntrials,.false.,mattyp)
            title='buffered hamiltonian for '//itoc(numdim)
     1             //'d dimensional hamiltonian'
            call hamil(z(ham(1)),z(ham(1)),
     1                 dum,cdum,
     2                 dum,cdum,
     3                 z(vtot),z(hbuf),z(hbuf),ia(ibuf),
     4                 ia(ind),z(diag),z(diag), lenbuf,
     5                 nmax(ug,1),nmax(ug,2),nmax(ug,3),
     6                 numdim,n3d,prbufh,nel,incore,title,mattyp)
            if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
               call cdvd(z(hbuf),ia(ibuf),z(diag),z(trials),
     1                   z(eig(1)),dum,dum,
     2                   z(vecl(1)),dum,dum,
     3                   z(vecr(1)),dum,dum,     
     4                   z(pvec),z(hpvec),z(vec),
     5                   z(bmat),z(bmatm),z(etmp),z(svec),z(sveca),
     6                   z(vec),z(eigtot),z(work(3)),z(work(4)),
     7                   lwork,scale,cnverg,thresh,
     8                   nmax(ug,1),nmax(ug,2),nmax(ug,3),
     9                   numdim,n3d,nroots,ntrials,niter,nvec,
     x                   lenbuf,hamd,incore,nel,prdvd,mattyp)
            elseif(mattyp.eq.'real-symmetric') then
               call rsdvd(z(hbuf),ia(ibuf),z(diag),z(trials),
     1                    z(eig(1)),dum,dum,
     2                    z(vecl(1)),dum,dum,z(psi),
     3                    z(pvec),z(hpvec),z(vec),
     4                    z(bmat),z(bmatm),z(etmp),z(svec),z(vec),
     5                    z(eigtot),z(work(3)),z(work(4)),
     6                    scale,cnverg,thresh,
     7                    nmax(ug,1),nmax(ug,2),nmax(ug,3),
     8                    numdim,n3d,nroots,ntrials,nattim,niter,
     9                    nvec,lenbuf,hamd,incore,nel,prdvd,mattyp)
            endif          
         elseif(numdim.eq.2) then
            call vtrial(z(ham(1)),z(ham(1)),z(ham(2)),z(ham(2)),
     1                  dum,cdum,
     2                  z(vecl(1)),z(vecl(1)),z(vecl(2)),z(vecl(2)),
     3                  dum,cdum,
     4                  z(eig(1)),z(eig(1)),z(eig(2)),z(eig(2)),
     5                  dum,cdum,
     6                  z(trials),z(trials),z(etrial),z(etrial),
     7                  ia(idum),nmax(ug,1),nmax(ug,2),nmax(ug,3),
     8                  numdim,n3d,ntrials,.false.,mattyp)                 
            title='buffered hamiltonian for '//itoc(numdim)
     1             //'d dimensional hamiltonian'
            call hamil(z(ham(1)),z(ham(1)),
     1                 z(ham(2)),z(ham(2)),
     2                 dum,cdum,
     3                 z(vtot),z(hbuf),z(hbuf),ia(ibuf),
     4                 ia(ind),z(diag),z(diag), lenbuf,
     5                 nmax(ug,1),nmax(ug,2),nmax(ug,3),numdim,n3d,
     6                 prbufh,nel,incore,title,mattyp)
            if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
               call cdvd(z(hbuf),ia(ibuf),z(diag),z(trials),
     1                   z(eig(1)),z(eig(2)),dum,
     2                   z(vecl(1)),z(vecl(2)),dum,
     3                   z(vecr(1)),z(vecr(2)),dum,
     4                   z(pvec),z(hpvec),z(vec),
     5                   z(bmat),z(bmatm),z(etmp),z(svec),z(sveca),
     6                   z(vec),z(eigtot),z(work(3)),z(work(4)),
     7                   lwork,scale,cnverg,thresh,
     8                   nmax(ug,1),nmax(ug,2),nmax(ug,3),
     9                   numdim,n3d,nroots,ntrials,niter,nvec,lenbuf,
     x                   hamd,incore,nel,prdvd,mattyp)
            elseif(mattyp.eq.'real-symmetric') then     
               call rsdvd(z(hbuf),ia(ibuf),z(diag),z(trials),
     1                    z(eig(1)),z(eig(2)),dum,
     2                    z(vecl(1)),z(vecl(2)),dum,z(psi),
     3                    z(pvec),z(hpvec),z(vec),
     4                    z(bmat),z(bmatm),z(etmp),z(svec),z(vec),
     5                    z(eigtot),z(work(3)),z(work(4)),scale,
     6                    cnverg,thresh,
     7                    nmax(ug,1),nmax(ug,2),nmax(ug,3),
     8                    numdim,n3d,nroots,ntrials,nattim,niter,
     9                    nvec,lenbuf,hamd,incore,nel,prdvd,mattyp)
            endif      
         elseif(numdim.eq.3) then
            call vtrial(z(ham(1)),z(ham(1)),z(ham(2)),z(ham(2)),
     1                  z(ham(3)),z(ham(3)),
     2                  z(vecl(1)),z(vecl(1)),z(vecl(2)),z(vecl(2)),
     3                  z(vecl(3)),z(vecl(3)),
     4                  z(eig(1)),z(eig(1)),z(eig(2)),z(eig(2)),
     5                  z(eig(3)),z(eig(3)),
     6                  z(trials),z(trials),z(etrial),z(etrial),
     7                  ia(idum),nmax(ug,1),nmax(ug,2),nmax(ug,3),
     8                  numdim,n3d,ntrials,.false.,mattyp)         
            title='buffered hamiltonian for '//itoc(numdim)
     1             //'d dimensional hamiltonian'
            call hamil(z(ham(1)),z(ham(1)),
     1                 z(ham(2)),z(ham(2)),
     2                 z(ham(3)),z(ham(3)),
     3                 z(vtot),z(hbuf),z(hbuf),ia(ibuf),
     4                 ia(ind),z(diag),z(diag), lenbuf,
     5                 nmax(ug,1),nmax(ug,2),nmax(ug,3),numdim,n3d,
     6                 prbufh,nel,incore,title,mattyp)               
            if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
               call cdvd(z(hbuf),ia(ibuf),z(diag),z(trials),
     1                   z(eig(1)),z(eig(2)),z(eig(2)),
     2                   z(vecl(1)),z(vecl(2)),z(vecl(3)),
     3                   z(vecr(1)),z(vecr(2)),z(vecr(3)),
     4                   z(pvec),z(hpvec),z(vec),
     5                   z(bmat),z(bmatm),z(etmp),z(svec),z(sveca),
     6                   z(vec),z(eigtot),z(work(3)),z(work(4)),
     7                   lwork,scale,cnverg,thresh,
     8                   nmax(ug,1),nmax(ug,2),nmax(ug,3),
     9                   numdim,n3d,nroots,ntrials,niter,nvec,lenbuf,
     x                   hamd,incore,nel,prdvd,mattyp)
            elseif(mattyp.eq.'real-symmetric') then
               call rsdvd(z(hbuf),ia(ibuf),z(diag),z(trials),
     1                    z(eig(1)),z(eig(2)),z(eig(2)),
     2                    z(vecl(1)),z(vecl(2)),z(vecl(3)),z(psi),
     3                    z(pvec),z(hpvec),z(vec),
     4                    z(bmat),z(bmatm),z(etmp),z(svec),z(vec),
     5                    z(eigtot),z(work(3)),z(work(4)),
     6                    scale,cnverg,thresh,
     7                    nmax(ug,1),nmax(ug,2),nmax(ug,3),
     7                    numdim,n3d,nroots,ntrials,nattim,niter,
     8                    nvec,lenbuf,hamd,incore,nel,prdvd,mattyp)
            endif     
         endif
      else    
        if(numdim.eq.1) then
            call lschr(z(ham(1)),z(ham(1)),
     1                 dum,cdum,
     2                 dum,cdum,
     3                 z(vtot),z(hamtot),z(hamtot),
     4                 z(eigtot),z(eigtot),
     5                 z(bigvec),z(bigvec),z(bigvec2),
     6                 z(bigvec),z(bigvec2),
     7                 z(work(1)),lwork,z(work(2)),
     8                 ndum,n3d,nmax(ug,1),nmax(ug,2),nmax(ug,3),
     9                 numdim,.true.,.true.,prnht,mattyp)
         elseif(numdim.eq.2) then
            call lschr(z(ham(1)),z(ham(1)),
     1                 z(ham(2)),z(ham(2)),
     2                 dum,cdum,
     3                 z(vtot),z(hamtot),z(hamtot),
     4                 z(eigtot),z(eigtot),
     5                 z(bigvec),z(bigvec),z(bigvec2),
     6                 z(bigvec),z(bigvec2),
     7                 z(work(1)),lwork,z(work(2)),
     8                 ia(ind),n3d,nmax(ug,1),nmax(ug,2),nmax(ug,3),
     9                 numdim,.true.,.true.,prnht,mattyp)                   
         elseif(numdim.eq.3) then
            call lschr(z(ham(1)),z(ham(1)),
     1                 z(ham(2)),z(ham(2)),
     2                 z(ham(3)),z(ham(3)),
     3                 z(vtot),z(hamtot),z(hamtot),
     4                 z(eigtot),z(eigtot),
     5                 z(bigvec),z(bigvec),z(bigvec2),
     6                 z(bigvec),z(bigvec2),
     7                 z(work(1)),lwork,z(work(2)),
     8                 ia(ind),n3d,nmax(ug,1),nmax(ug,2),nmax(ug,3),
     9                 numdim,.true.,.true.,prnht,mattyp)                      
         endif
         if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
            call cvscal(z(eigtot),z(eigtot),scale,nroots)
            title='eigenvalues of hamiltonian'
            call prntcm(title,z(eigtot),nroots,1,nroots,1,iout)
         elseif(mattyp.eq.'real-symmetric') then
            call vscale(z(eigtot),z(eigtot),scale,nroots)
            title='eigenvalues of hamiltonian'
            call prntfm(title,z(eigtot),nroots,1,nroots,1,iout)
         endif
         if(scatt) then
            do 200 i=1,nen
               call rmtrx(z(pn(ug,1)),z(eigc(ug,1)),z(eigtot),
     1                    z(bigvec),energy(i),typot,n3d)
  200       continue         
         endif
      endif
      call manmem(-ngot,pnt,idum,'pde',idum)
      call chainx(0)               
      stop
 1    format(/,15x,'coordinate system = ',a16)
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
 9    format(/,15x,'iterative diagonalization information',/,/,5x,
     1             'number of roots                    = ',i3,/,5x,
     2             'number of roots at a time          = ',i3,/,5x
     3             'overlap tolerance                  = ',e15.8,/,5x,
     4             'convergence criterion              = ',e15.8,/,5x,
     5             'maximum number of iterations       = ',i6)
 11   format(/,15x,'direct diagonalization',/,/,5x,
     1             'size of matrix  = ',i3,/,5x,
     2             'number of roots = ',i3)
 12   format(/,5x,'dimension = ',i2,2x,'grid = ',i3)
      end
