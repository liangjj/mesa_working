*deck timdvr.f 
c***begin prologue     timdvr
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, orthogonal polynomial
c***                   bose-einstein
c***author             schneider, b. i.(nsf)
c***source             tprop
c***purpose            driver for solution of time dependent
c***                   schroedinger equation using time based dvr.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       timdvr
      program timdvr
c
      implicit integer (a-z)
      parameter ( maxint=200, mgr=20 )
      character*4096 ops
      character*2 itoc, ic
      character*8 prtflg, qtyp
      character*8 type
      character*16 key
      character*80 cpass, titphr, title, chrkey, precon
      character*24 coord
      character*800 card
      character*24 i0stat
      character*128 filbec, filham
      character*24 timh0, sptimv
      character*2 atom
      character*3 timkey, trapon, hydro
      logical dollar, logkey
      logical prnt, mprnt
      logical toau, test, useau, totrap
      logical itsolv, itdiag, mltgrd
      logical spac, h0test, hdiag
      real*8 basis, array, soln, driver, scr
      real*8 fpkey, omegat, amass 
      real*8 scatl, pi, u0
      real*8 hbar, omega, thresh, cnverg, eps
      real*8 scale
      real*8 massau, lenau, timau
      real*8 energy, bigomg, lenscl, enscl
      real*8 sctlen, eta
      dimension nmax(mgr,4), npt(mgr,4), nsubg(4), ug(4)
      dimension nptr(4,mgr), nmaxr(4,mgr), n3d(mgr), ntot(mgr), m(mgr)
      dimension maxvec(mgr), maxit(mgr), energy(mgr)
      dimension fl(4), fr(4), ngot(20)
      dimension q(4), qm(4), wtm(4)
      dimension pm(4), dpm(4), ddpm(4)
      dimension ham(4), v(4), trn(4), eig(4), fwadd(2,4)
      dimension prnt(30), mprnt(30), key(4)
      dimension omegat(4), qtyp(4)
      dimension type(2)
      common/io/inp, iout      
      common/punch/pun
      pointer(pdsk,dskp(1))
      pointer(pfun,funp(1))
      pointer(pone,onep(1))
      pointer(pbasis,basis(1))
      pointer(parray,array(1))
      pointer(pscr,scr(1))
      pointer(pt,pnt(1))
      pointer(psoln,soln(1))
      pointer(pdrver,driver(1))
      pointer(phamil,hamil(1)), (phamil,ihamil(1))
      pointer(plufac,lufac(1)), (plufac,indlu(1))
      data pi/3.1415926535897932384d0/
      data massau, lenau, timau / 9.109558d-31, 5.291771d-11, 
     1                            2.41888d-17 /
c     hbar in joule-sec      
      data hbar/1.054592d-34/
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
      pnch=logkey(ops,'punch=on',.false.,' ')
c
c         set spatial dimensionality of problem and program options 
c
      spdim=intkey(ops,'number-of-space-dimensions',1,' ')
      ntreg=intkey(ops,'number-of-time-regions',1,' ') 
      dim=spdim+1
      coord=chrkey(ops,'coordinate-system','cartesian',' ')
      toau=logkey(ops,'to-atomic-units',.false.,' ')
      trapon=chrkey(ops,'trap','off',' ')
      timkey=chrkey(ops,'m6295=time','on',' ')
      hydro=chrkey(ops,'trap=hydro','off',' ')
      mltgrd=logkey(ops,'multigrid-procedure=on',.false.,' ')
      topg=intkey(ops,'top-grid',1,' ')
      key(1)='$trap=off'
      if(trapon.eq.'on') then
         key(1)='$trap=on'
      endif         
      key(2)='$time=off'
      if(timkey.eq.'on') then
         key(2)='$time=on'
         ug(1)=intkey(ops,'time-grid-number',1,' ')
      endif     
      call intarr(ops,'spatial-grid-number',ug(2),spdim,' ')
      totrap=logkey(ops,'use-trap-units',.false.,' ')
      useau=logkey(ops,'use-atomic-units',.false.,' ')
      prnt(1)=logkey(ops,'print=m6295=all',.false.,' ')
      prnt(2)=logkey(ops,'print=m6295=grids',.false.,' ')
      prnt(3)=logkey(ops,'print=m6295=one-body-hamiltonian',.false.,' ')
      prnt(4)=logkey(ops,'print=m6295=initial-psi',.false.,' ')
      prnt(5)=logkey(ops,'print=m6295=driving-terms',.false.,' ')
      prnt(6)=logkey(ops,'print=m6295=full-hamiltonian',.false.,' ')
      prnt(7)=logkey(ops,'print=m6295=solution',.false.,' ')
      prnt(8)=logkey(ops,'print=m6295=block-factors',.false.,' ')
      do 10 i=2,8
         if(prnt(1)) then
            prnt(i)=.true.
         endif
 10   continue   
      itsolv=logkey(ops,'iterative-linear-system-solve',.false.,' ')
      h0test=logkey(ops,'test-h0',.false.,' ')
      hdiag=logkey(ops,'diagonalize-hamiltonian',.false.,' ')
      type(1)='unknown'
      type(2)=chrkey(ops,'open-ham','unknown',' ')
      write(iout,1) spdim
c
      write(iout,2)
      call iosys ('read character "bec filename" from rwf',
     1             -1,0,0,filbec)
      call iosys ('open bec as '//type(1),0,0,0,filbec)
      call iosys('rewind all on bec read-and-write',0,0,0,' ')
c
      qtyp(2)=chrkey(ops,'space-dimension-1','x',' ')
      qtyp(3)=chrkey(ops,'space-dimension-2','y',' ')
      qtyp(4)=chrkey(ops,'space-dimension-3','z',' ')
c
c               [ get the pointers and the memory needed for]
c               [           the arrays                      ]
c
      call pntmem(mgr,dim,pdsk,pfun,pone,ngot)
      call setpnt(q,qm,wtm,pm,dpm,ddpm,ham,v,trn,eig,mgr,dim)
c
c  
      call spamem(qtyp,npt,nmax,mgr,nsubg,dim,ntreg,pbasis,
     1            parray,ngot(4),fwadd)
c
c     calculate the memory pointers for all space dimensions
c
c
c             read in the disk pointers to all the arrays.
c
      ic=itoc(2)
      titphr=qtyp(2)//' pointer'
      call rdpnt(qtyp(2),fl(2),fr(2),npt(1,2),nmax(1,2),
     1           nsubg(2),dskp(q(2)), dskp(qm(2)),dskp(wtm(2)),
     2           dskp(pm(2)),dskp(dpm(2)),dskp(ddpm(2)),
     3           mgr,titphr,'space')
c
c             read in spatial variables
c
      call fillz(funp(q(2)),funp(qm(2)),funp(wtm(2)),
     1           funp(pm(2)),funp(dpm(2)),funp(ddpm(2)),
     2           onep(ham(2)),onep(v(2)),onep(trn(2)),onep(eig(2)),
     3           nsubg(2),npt(1,2),nmax(1,2),mgr,fwadd(1,2))
      titphr=qtyp(2)//' title'
      call rdfun(qtyp(2),npt(1,2),nmax(1,2),nsubg(2),
     1           mgr,dskp(q(2)),dskp(qm(2)),dskp(wtm(2)),
     2           dskp(pm(2)),dskp(dpm(2)),dskp(ddpm(2)),
     3           funp(q(2)),funp(qm(2)),funp(wtm(2)),
     4           funp(pm(2)),funp(dpm(2)),funp(ddpm(2)),
     5           pbasis,titphr,prnt(2))
c
c             read in the one body operators
c   
      call rdke(qtyp(2),nmax(1,2),nsubg(2),mgr,
     1          onep(ham(2)),onep(v(2)),onep(trn(2)),onep(eig(2)),
     2          parray,'space',prnt(3))
c
c               [    second coordinate    ]
c
      if(spdim.gt.1) then
c
c             get the disk pointers to all the arrays.
c
         ic=itoc(3)
         titphr=qtyp(3)//' pointer'
         call rdpnt(qtyp(3),fl(3),fr(3),npt(1,3),nmax(1,3),
     1              nsubg(3),dskp(q(3)),dskp(qm(3)),dskp(wtm(3)),
     2              dskp(pm(3)),dskp(dpm(3)),dskp(ddpm(3)),
     3              mgr,titphr,'space')
c
c             read in spatial variables
c
         call fillz(funp(q(3)),funp(qm(3)),funp(wtm(3)),
     1              funp(pm(3)),funp(dpm(3)),funp(ddpm(3)),
     2              onep(ham(3)),onep(v(3)),onep(trn(3)),onep(eig(3)),
     3              nsubg(3),npt(1,3),nmax(1,3),mgr,fwadd(1,3))
         titphr=qtyp(3)//' title'
         call rdfun(qtyp(3),npt(1,3),nmax(1,3),nsubg(3),
     1              mgr,dskp(q(3)),dskp(qm(3)),dskp(wtm(3)),
     2              dskp(pm(3)),dskp(dpm(3)),dskp(ddpm(3)),
     3              funp(q(3)),funp(qm(3)),funp(wtm(3)),
     4              funp(pm(3)),funp(dpm(3)),funp(ddpm(3)),
     5              pbasis,titphr,prnt(2))
c
c             read in the one body operators
c   
         call rdke(qtyp(3),nmax(1,3),nsubg(3),mgr,
     1             onep(ham(3)),onep(v(3)),onep(trn(3)),
     2             onep(eig(3)),parray,'space',prnt(3))
      endif
c
c
c               [    third coordinate    ]
c
      if(spdim.gt.2) then
c
c             get the disk pointers to all the arrays.
c
         ic=itoc(4)
         titphr=qtyp(4)//' pointer'
         call rdpnt(qtyp(4),fl(4),fr(4),npt(1,4),nmax(1,4),
     1              nsubg(4),dskp(q(4)),dskp(qm(4)),dskp(wtm(4)),
     2              dskp(pm(4)),dskp(dpm(4)),dskp(ddpm(4)),
     3              mgr,titphr,'space')
c
c             read in spatial variables
c       
         call fillz(funp(q(4)),funp(qm(4)),funp(wtm(4)),
     1              funp(pm(4)),funp(dpm(4)),funp(ddpm(4)),
     2              onep(ham(4)),onep(v(4)),onep(trn(4)),onep(eig(4)),
     3              nsubg(4),npt(1,4),nmax(1,4),mgr,fwadd(1,4))
         titphr=qtyp(4)//' title'
         call rdfun(qtyp(4),npt(1,4),nmax(1,4),nsubg(4),
     1              mgr,dskp(q(4)),dskp(qm(4)),dskp(wtm(4)),
     2              dskp(pm(4)),dskp(dpm(4)),dskp(ddpm(4)),
     3              funp(q(4)),funp(qm(4)),funp(wtm(4)),
     4              funp(pm(4)),funp(dpm(4)),funp(ddpm(4)),
     5              pbasis,titphr,prnt(2))
c
c             read in the one body operators
c   
         call rdke(qtyp(4),nmax(1,4),nsubg(4),mgr,
     1             onep(ham(4)),onep(v(4)),onep(trn(4)),
     2             onep(eig(4)),parray,'space',prnt(3))
      endif
c
c         set trap configuration, numbers of atoms, various constants
c         and other physical parameters
c       
      if(toau) then
         hbar=1.d0
      endif
      if(useau) then
         hbar=1.d0
      endif
      amass=1.d0
      if ( dollar(key(1),card,cpass,inp) ) then
           atom=chrkey(card,'atom','cs',' ')
           sctlen=fpkey(card,'scattering-length',52.d0,' ')
           natom=intkey(card,'number-of-atoms',0,' ')
           call atmdat(atom,amass,omegat(2))
           if(toau) then
              write(iout,*) '          converting to atomic units'
              omegat(2)=omegat(2)*timau
              omegat(3)=omegat(3)*timau
              omegat(4)=omegat(4)*timau
              amass=amass/massau
           endif
           if(useau) then
              write(iout,*) '          assuming atomic units'
              amass=1.d0
              omegat(2)=1.d0
              omegat(3)=1.d0
              omegat(4)=1.d0
           endif
           if(totrap) then
              bigomg=max(omegat(2),omegat(3),omegat(4))
              omegat(2)=omegat(2)/bigomg
              omegat(3)=omegat(3)/bigomg
              omegat(4)=omegat(4)/bigomg
              if(omegat(2).eq.1.d0) then
                 ibig=2
              endif
              if(omegat(3).eq.1.d0) then
                 ibig=3
              endif
              if(omegat(4).eq.1.d0) then
                 ibig=4
              endif
              lenscl=sqrt(hbar/(amass*bigomg))
              escal=hbar*bigomg
              eta=sctlen*natom/lenscl
           endif
           if(spdim.eq.1) then
              write(iout,3) amass, omegat(2), natom, sctlen
           elseif(spdim.eq.2) then
              write(iout,4) amass, ( omegat(i), i=2,3 ), natom, sctlen
           elseif(spdim.eq.3) then
              write(iout,5) amass, ( omegat(i), i=2,4 ), natom, sctlen
           endif                           
      endif
      if( dollar(key(2),card,cpass,inp) ) then
          timh0=chrkey(card,'pure-time-perturbation','none',' ')
          spac=logkey(card,'no-spatial-hamiltonian',.false.,' ')
          i0stat=chrkey(card,'driver','state-vector',' ')
          state=intkey(card,'initial-state',0,' ')
          sptimv=chrkey(card,'space-time-perturbation','none',' ')
          write(iout,6) spac, timh0, sptimv, i0stat, state, ntreg
          state=state+1  
      endif
c
c             Begin Spatial Code
c     
      maxd=0
      do 20 i=2,dim
         maxd=max(nmax(ug(i),i),maxd)
 20   continue
      if(h0test) then
         scrwd=wpadti(1+2*maxd+maxd*maxd)   
         call memory(scrwd,pscr,got,'scratch',0)
      endif
      do 30 i=2,dim
         call hamxyz(onep(ham(i)),onep(v(i)),scr(1),
     1               scr(1+maxd*maxd),nmax(1,i),ug(i),
     2               mgr,maxd,h0test,parray,prnt(3),spac)
 30   continue   
      if(h0test) then
         call memory(-got,pscr,idum,'scratch',idum)
      endif
c
c         End Spatial Code
c
c         Begin Time Calculation
c
c
c             set linear system procedures
c      
      if(itsolv) then
         if( dollar('$gmres',card,cpass,inp) ) then
             call lindat(card,cpass,cnverg,thresh,eps,precon,nblck,
     1                   prnt(9),filham,type(2))
         endif
      endif
      if(mltgrd) then
         call memory(4*mgr,pt,nwpt,'tpnt',0)
         pvt=1
         pdrv=pvt+mgr
         psln=pdrv+mgr
         prhs=psln+mgr
      endif
      do 40 i=2,dim
         do 50 j=1,nsubg(i)
            nptr(i,j)=npt(j,i)
            nmaxr(i,j)=nmax(j,i)
 50      continue
 40   continue   
      do 60 tim=1,ntreg
c
c             get the disk pointers to all the arrays.
c
         ic=itoc(tim)
         qtyp(1)='t-'//ic
         call iosys('read integer "no. subgrids for '
     1              //qtyp(1)//'" from bec',1,nsubg(1),0,' ')
         call iosys('read integer "no. points for '
     1              //qtyp(1)//'" from bec',nsubg(1),
     2                npt(1,1),0,' ')
         call iosys('read integer "mod. no. points '//
     1              'for '//qtyp(1)//'" from bec',nsubg(1),
     2               nmax(1,1),0,' ')
         do 70 i=1,nsubg(1)
            nptr(1,i)=npt(i,1)
            nmaxr(1,i)=nmax(i,1)
 70      continue   
         titphr='t pointer region '//ic
         call rdpnt(qtyp(1),fl(1),fr(1),npt(1,1),nmax(1,1),
     1              nsubg(1),dskp(q(1)),dskp(qm(1)),dskp(wtm(1)),
     2              dskp(pm(1)),dskp(dpm(1)),dskp(ddpm(1)),
     3              mgr,titphr,'time')
c
c             read in spatial variables
c
         titphr='t title region '//ic
         call fillz(funp(q(1)),funp(qm(1)),funp(wtm(1)),
     1              funp(pm(1)),funp(dpm(1)),funp(ddpm(1)),
     2              onep(ham(1)),onep(v(1)),onep(trn(1)),onep(eig(1)),
     3              nsubg(1),npt(1,1),nmax(1,1),mgr,fwadd(1,1))   
         call rdfun(qtyp(1),npt(1,1),nmax(1,1),nsubg(1),
     1              mgr,dskp(q(1)),dskp(qm(1)),dskp(wtm(1)),
     2              dskp(pm(1)),dskp(dpm(1)),dskp(ddpm(1)),
     3              funp(q(1)),funp(qm(1)),funp(wtm(1)),
     4              funp(pm(1)),funp(dpm(1)),funp(ddpm(1)),pbasis,
     5              titphr,prnt(2))
c
c             read in the one body operators
c   
         call rdke(qtyp(1),nmax(1,1),nsubg(1),mgr,
     1             onep(ham(1)),onep(v(1)),onep(trn(1)),onep(eig(1)),
     2             parray,'time',prnt(3))
         call hamt(onep(ham(1)),onep(v(1)),nmax(1,1),ug(1),mgr,
     1             h0test,parray,prnt(3))
c
c            either use standard or multigrid approach
c
         if(.not.mltgrd) then
            call psit(funp(q(1)),funp(qm(1)),funp(wtm(1)),
     1                funp(pm(1)),onep(ham(1)),onep(v(1)),
     2                onep(trn(1)),onep(eig(1)),pbasis,parray,spdim,
     3                dim,npt,nmax,ug,mgr,spac,tim,timh0,i0stat,
     4                sptimv,state,coord,itsolv,cnverg,thresh,eps,
     5                precon,nblck,card,prnt,hdiag,pnch)
         else
c            call mgset(funp(q(1)),funp(qm(1)),funp(wtm(1)),
c     1                 funp(pm(1)),onep(ham(1)),onep(v(1)),onep(trn(1)),
c     2                 onep(eig(1)),pnt(psln),pnt(pvt),pnt(pdrv),
c     3                 pnt(prhs),pbasis,parray,pvoft,pdrver,
c     4                 psoln,pprhs,phamil,plufac,energy,spdim,dim,
c     5                 nptr,nmaxr,n3d,ntot,m,topg,mgr,spac,
c     6                 cnverg,thresh,eps,maxit,maxvec,nblck,
c     7                 mcycle,precon,tim,timh0,i0stat,sptimv,
c     8                 state,coord,prnt,mprnt,ngot(6))
c
c           the one body potentials are no longer needed since they
c           have been added into the one body hamiltonian and voft.
c
c            call fullmg(funp(pm(1)),onep(ham(1)),onep(trn(1)),
c     1                  onep(eig(1)),pnt(pvt),pnt(pdrv),pnt(prhs),
c     2                  pnt(psln),pbasis,parray,pvoft,pdrver,
c     3                  psoln,pprhs,phamil,plufac,energy,
c     4                  spdim,dim,nmaxr,n3d,ntot,m,topg,mgr,
c     5                  cnverg,thresh,eps,maxit,maxvec,nblck,mcycle,
c     6                  precon,timh0,i0stat,prnt(9),mprnt,ngot(6))
         endif
 60   continue   
      if(mltgrd) then
         call memory(-nwpt,pt,idum,'tpnt',idum)
      endif
      call chainx(0)               
      stop
 1    format(/,20x,'time-dependent basis function code',//,20x,
     1             'number of spatial dimensions = ',i1)      
 2    format(/,15x,'calculation = solve time-dependent schrodinger'
     1             ' equation')
 3    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                      = ',e15.8,/,5x,
     2             'trap frequency                   = ',e15.8,/,5x,
     3             'number of atoms in trap          = ',i8,/,5x,
     4             'scattering length                = ',e15.8)
 4    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                      = ',e15.8,/,5x,
     2             'q1-trap frequency                = ',e15.8,/,5x,
     3             'number of atoms in trap          = ',i8,/,5x,
     4             'scattering length                = ',e15.8)
 5    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                      = ',e15.8,/,5x,
     2             'q1-trap frequency                = ',e15.8,/,5x,
     3             'q2-trap frequency                = ',e15.8,/,5x,
     4             'q3-trap frequency                = ',e15.8,/,5x,
     5             'number of atoms in trap          = ',i8,/,5x,
     6             'scattering length                = ',e15.8)
 6    format(/,5x,'time-dependent data',/,5x,
     1            'no spatial hamiltonian   = ',l1,/,5x,
     2            'pure-time perturbation   = ',a24,/,5x,
     3            'space-time-perturbation  = ',a24,/,5x,
     4            'driver                   = ',a24,/,5x,
     5            'state                    = ',i2,/,5x,
     6            'number of time intervals = ',i3)
      end
