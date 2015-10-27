*deck mgdrv.f 
c***begin prologue     mgdrv
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, orthogonal polynomial
c***                   bose-einstein
c***author             schneider, b. i.(nsf)
c***source             tprop
c***purpose            solve time dependent schroedinger equation 
c***                   using a multigrid method based on dvr basis sets.
c***description        
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       mgdrv
      subroutine mgdrv(q,qm,wtm,pm,dpm,ddpm,ham,v,trn,eig,energy,
     1                 pbasis,parray,pvt,pdrv,psln,ppsi0,spdim,
     2                 dim,npt,nmax,n3d,ntot,m,topg,mgr,spac,
     3                 maxit,maxvec,tim,timh0,i0stat,sptimv,state,
     4                 coord,prnt)
c
      implicit integer (a-z)
      real*8 basis, array, voft, driver, soln, psi0, resid, lufac
      real*8 hamil, fac, t0, gett0, scr
      real*8 hbar, omega, thresh, cnverg, eps
      real*8 energy, bigomg, lenscl, enscl
      real*8 sctlen, eta, scale, fpkey, dum
      character*(*) timh0, sptimv, i0stat, coord
      character*80 title, cpass, chrkey, precon
      character*320 card
      character*128 filham
      character*3 ok, itoc
      logical spac, prnt, dollar, logkey, mprnt
      dimension q(mgr,dim), qm(mgr,dim), wtm(mgr,dim)
      dimension pm(mgr,mgr,dim), dpm(mgr,mgr,dim), ddpm(mgr,mgr,dim)    
      dimension npt(4,mgr), nmax(4,mgr), maxit(mgr), maxvec(mgr)
      dimension ham(mgr,dim), v(mgr,dim), trn(mgr,dim), eig(mgr,dim)
      dimension pvt(mgr), pdrv(mgr), psln(mgr), ppsi0(mgr)
      dimension n3d(mgr), ntot(mgr), m(mgr), ngot(10)
      dimension energy(mgr)
      dimension prnt(*), mprnt(10)
      common/io/inp, iout      
      pointer(pbasis,basis(1))
      pointer(parray,array(1))
      pointer(pvoft,voft(1))
      pointer(pdrver,driver(1))
      pointer(psoln,soln(1))
      pointer(pres,resid(1))  
      pointer(pppsi0,psi0(1))
      pointer(phamil,hamil(1)), (phamil,ihamil(1))
      pointer(pind,ind(1))
      pointer(plufac,lufac(1)), (plufac,indlu(1))
      pointer(pscr,scr(1)), (pscr,iscr(1))
c
c
c
      if( dollar('$multigrid',card,cpass,inp) ) then
          call lindat(card,cpass,cnverg,thresh,eps,cpass,nblck,
     1                prnt(9),filham,'scratch')
          mprnt(1)=logkey(card,'print=mgrid=coarse-grid-solution',
     1                    .false.,' ')
          mprnt(2)=logkey(card,'print=mgrid=interpolated-initial-'
     1                         //'guess',.false.,' ')
          mprnt(3)=logkey(card,'print=mgrid=interpolated-right-'
     1                         //'hand-side',.false.,' ')
          mprnt(4)=logkey(card,'print=mgrid=interpolated-solution',
     1                    .false.,' ')
          mprnt(5)=logkey(card,'print=mgrid=all',.false.,' ')
          if(mprnt(5)) then
             mprnt(1)=.true.
             mprnt(2)=.true.
             mprnt(3)=.true.
             mprnt(4)=.true.
          endif
      endif
      call intarr(card,'maximum-number-of-iterations',maxit,mgr,' ')
      call intarr(card,'maximum-number-of-vectors',maxvec,mgr,' ')
      nrhs=intkey(card,'number-of-right-hand-sides',1,' ')
      ntrial=intkey(card,'number-of-trial-vectors',1,' ')
      mcycle=intkey(card,'number-of-multigrid-cycles',1,' ')
      precon=chrkey(card,'preconditioner','none',' ')
      if(precon.ne.'none') then
         nblck=intkey(card,'maximum-size-of-preconditioning-block',
     1                      200,' ')                      
      endif
      write(iout,1)
      bigm=0
      bign=0
      do 10 grid=1,topg
c
c        find out the size of the problem for this grid
c
         n3d(grid)=1
         do 20 ndim=2,dim
            n3d(grid)=n3d(grid)*nmax(ndim,grid)
 20      continue   
         ntot(grid)=n3d(grid)*nmax(1,grid)   
         m(grid)=2*ntot(grid)
         bigm=max(bigm,m(grid))
         bign=max(bign,ntot(grid))
 10   continue
c
c
c     check to see if the problem can be solved exactly on grid one.
c
      call chkmem(m(1),ntot(1),ok)
      if(ok.ne.'ok') then
         write(iout,2)                  
         call lnkerr('quit multigrid')         
      endif
      nwvoft=1
      nwdrver=1
      nwsol=1
      nwpsi0=1
      do 30 grid=1,topg
         pvt(grid)=nwvoft      
         pdrv(grid)=nwdrver
         psln(grid)=nwsol
         ppsi0(grid)=nwpsi0
         nwvoft=nwvoft+ntot(grid)
         nwdrver=nwdrver+m(grid)
         nwsol=nwsol+m(grid)
         nwpsi0=nwpsi0+m(grid)
 30   continue   
c
c     get memory for all potentials, solutions, driving terms and residuals
c
      nwvoft=wpadti(nwvoft)
      nwdrver=wpadti(nwdrver)
      nwsol=wpadti(nwsol)
      call memory(nwvoft,pvoft,ngot(1),'voft',0)   
      call memory(nwdrver,pdrver,ngot(2),'driver',0)   
      call memory(nwsol,psoln,ngot(3),'soln',0)
      call memory(nwpsi0,pppsi0,ngot(4),'psi0',0)
      call memory(nwsol,pres,ngot(5),'resid',0)
c 
c     calculate potentials and driving terms on all grids.
c
      do 40 grdi=1,topg
c
c        make the one body hamiltonian
c
         do 50 i=2,dim
            hami=ham(grdi,i)
            poti=v(grdi,i)
            ni=nmax(i,grdi)
            call addd(array(hami),array(hami),array(poti),ni)
 50      continue   
c
         nxyz=nmax(2,grdi)
         nt=nmax(1,grdi)
         msze=m(grdi) 
         nsze=ntot(grdi)
         n=n3d(grdi)
         hx=ham(grdi,2)
         hy=ham(grdi,3)
         hz=ham(grdi,4)
         ht=ham(grdi,1)
         eigx=eig(grdi,2)
         eigy=eig(grdi,3)
         eigz=eig(grdi,4)
         wtx=wtm(grdi,2)
         wty=wtm(grdi,3)
         wtz=wtm(grdi,4)
         wtt=wtm(grdi,1)
         trnx=trn(grdi,2)
         trny=trn(grdi,3)
         trnz=trn(grdi,4)
         x=qm(grdi,2)
         y=qm(grdi,3)
         z=qm(grdi,4)
         t=qm(grdi,1)
         px=pm(grdi,grdi,2)
         py=pm(grdi,grdi,3)
         pz=pm(grdi,grdi,4)
         pt=pm(grdi,grdi,1)
         vtot=pvt(grdi)
         vt=v(grdi,1)
         drvi=pdrv(grdi)
         p0=ppsi0(grdi)
c
c----------------------------------------------------------------------c
c        Calculate the time dependent perturbation.  It consists of a  c
c        space and a time part for this grid. store it                 c
c----------------------------------------------------------------------c
         call vpert(voft(vtot),array(vt),basis(x),basis(y),basis(z),
     1              basis(t),n,nxyz,nt,spdim,sptimv,timh0)
c----------------------------------------------------------------------c
c
c----------------------------------------------------------------------c
c        Make the proper driving term                                  c
c----------------------------------------------------------------------c
         call mkpsi0(array(trnx),array(trny),array(trnz),
     1               array(eigx),array(eigy),array(eigz),
     2               basis(px),basis(py),basis(pz),
     3               basis(x),basis(y),basis(z),
     4               basis(wtx),basis(wty),basis(wtz),psi0(p0),
     5               energy(grdi),n,nxyz,state,spdim,
     6               coord,tim,i0stat,prnt(4))
         call mkhpsi(array(hx),array(hy),array(hz)),
     1               driver(drvi),psi0(p0),voft(vtot),
     2               basis(pt),basis(wtt),n,nx,ny,nz,nt,spac,
     3               spdim,prnt(5))
 40   continue   
c
c     memory for preconditioning
c
      if(precon.ne.'none') then
         if(nblck.ge.bign) then
            nblck=bign
         endif
         luh=1
         lujnd=wpadti(hsub+2*nblck*nblck)
         nwlu=lujnd+bign
         lutrial=1
         lutmp=lutrial+2*bign
         luind=wpadti(lutmp+bign)
         scrwds=luind+(spdim+1)*bign
         call memory(nwlu,plufac,ngot(6),'lufac',0)
         call memory(scrwds,pscr,nscr,'scratch',0)
         do 60 grd=2,topg
            nx=nmax(2,grd)
            ny=nmax(3,grd)
            nz=nmax(4,grd)
            nt=nmax(1,grd)
            nsze=ntot(grd)
            n=n3d(grd)
            hx=ham(grd,2)
            hy=ham(grd,3)
            hz=ham(grd,4)
            ht=ham(grd,1)
            vtot=pvt(grd)
            drvcur=pdrv(grd)
            call factors(lufac(luh),driver(drvcur),scr(lutrial),
     1                   scr(lutmp),array(hx),array(hy),array(hz),
     2                   array(ht),voft(vtot),iscr(luind),indlu(lujnd),
     3                   nx,ny,nz,nt,nsze,nblck,spdim,
     4                   grd,prnt(8).true.)
 60      continue   
         call memory(-nscr,pscr,idum,'scratch',idum)
      endif
c
c     calculate solution on the first grid exactly.  save factorization.
c
      nx=nmax(2,1)
      ny=nmax(3,1)
      nz=nmax(4,1)
      nt=nmax(1,1)
      msze=m(1) 
      nsze=ntot(1)
      n=n3d(1)
      hx=ham(1,2)
      hy=ham(1,3)
      hz=ham(1,4)
      ht=ham(1,1)
      x=qm(1,2)
      y=qm(1,3)
      z=qm(1,4)
      t=qm(1,1)
      px=pm(1,1,2)
      py=pm(1,1,3)
      pz=pm(1,1,4)
      pt=pm(1,1,1)
      vtot=pvt(1)
      drvcur=pdrv(1)
      solcur=psln(1)
      p0=ppsi0(1)
      t0=basis(q(1,1))
      h0=1
      h0pvt=wpadti( h0 + msze*msze ) 
      need=h0pvt+msze       
      call memory(need,phamil,ngot(7),'direct-solve',0)
      call memory(nsze,pind,scrwds,'scratch',0)
      call hamtot(hamil(h0),array(hx),array(hy),array(hz),array(ht),
     1            voft(vtot),ind,nx,ny,nz,nt,nsze,msze,spdim,prnt(6))
      call memory(-scrwds,pind,idum,'scratch',idum)
      call copy(driver(drvcur),soln(solcur),msze)
      call hamfac(hamil(h0),ihamil(h0pvt),msze)
      call hamslv(hamil(h0),soln(solcur),ihamil(h0pvt),msze,.false.)
      if(mprnt(1)) then
         title='coarse grid solution on grid = '//itoc(1)
         call prntrm(title,soln(solcur),msze,1,msze,1,iout)
      endif
      call pltwfn(soln(solcur),basis(x),basis(y),basis(z),basis(t),
     1            basis(px),basis(py),basis(pz),basis(pt),
     2            psi0(p0),t0,energy(1),spdim,nsze,
     3            nx,ny,nz,nt,timh0,i0stat)
c
c
      scale=1.d0    
      do 60 grdi=2,topg
         write(iout,3) grdi
c
c               get starting solution by interpolating from previously 
c               calculated solution.
c
c        call grdj the previous grid
c
         grdj=grdi-1
         nxi=nmax(2,grdi)
         nyi=nmax(3,grdi)
         nzi=nmax(4,grdi)
         nti=nmax(1,grdi)
         nszei=ntot(grdi)
         mi=m(grdi)
         nxj=nmax(2,grdj)
         nyj=nmax(3,grdj)
         nzj=nmax(4,grdj)
         ntj=nmax(1,grdj)
         nszej=ntot(grdj)
         pxij=pm(grdi,grdj,2)
         pyij=pm(grdi,grdj,3)
         pzij=pm(grdi,grdj,4)
         ptij=pm(grdi,grdj,1)
         pxii=pm(grdi,grdi,2)
         pyii=pm(grdi,grdi,3)
         pzii=pm(grdi,grdi,4)
         ptii=pm(grdi,grdi,1)
         soli=psln(grdj)
         solj=psln(grdi)
c
         call gi2gj(soln(solj),soln(soli),            
     1              basis(ptij),basis(pxij),basis(pyij),basis(pzij),
     2              nti,nxi,nyi,nzi,ntj,nxj,nyj,nzj,nszei,nszej,2,dim)
         call vivo(soln(soli),soln(soli),
     1             basis(ptii),basis(pxii),basis(pyi),basis(pzii),
     2             nti,nxi,nyi,nzi,nszei,2,dim)
         if(mprnt(2)) then
            title='interpolated initial guess on grid = '//itoc(grdi)
            call prntrm(title,soln(soli),mi,1,mi,1,iout)
         endif
c        
         do 70 grdj=grdi,2,-1
c
c           initiate the davidson procedure to reduce the error
c
            maxit(grdj)=min(maxit(grdj),m(grdj))
            maxvec(grdj)=min(maxvec(grdj),m(grdj))
c
c
            nxj=nmax(2,grdj)
            nyj=nmax(3,grdj)
            nzj=nmax(4,grdj)
            ntj=nmax(1,grdj)
            mj=m(grdj)
            nj=n3d(grdj)
            nszej=ntot(grdj)
            hxj=ham(grdj,2)
            hyj=ham(grdj,3)
            hzj=ham(grdj,4)
            htj=ham(grdj,1)
            vtotj=pvt(grdj)
            eigxj=eig(2,grdj)
            eigyj=eig(3,grdj)
            eigzj=eig(4,grdj)
            trnxj=trn(2,grdj)
            trnyj=trn(3,grdj)
            trnzj=trn(4,grdj)
            solj=psln(grdj)
            drvj=pdrv(grdj)            
            title='entering main davidson routine on grid = '
     1             //itoc(grdj)
            call qcklin(soln(solj),driver(drvj),resid(solj),
     1                  array(hxj),array(hyj),array(hzj),array(htj),
     2                  eig(eigxj),eig(eigyj),eig(eigzj),
     3                  array(trnxj),array(trnyj),array(trnzj),
     4                  voft(vtotj),lufac(luh),indlu(lujnd),
     5                  scale,cnverg,thresh,eps,nxj,nyj,nzj,ntj,
     6                  spdim,mj,nj,nszej,maxit,maxvec,precon,title,
     7                  prnt(9),mprnt)
         call pltwfn(soln(psln(grdi)),
     1               basis(qm(grdi,2)),basis(qm(grdi,3)),
     2               basis(qm(grdi,4)),basis(qm(grdi,1)),
     3               basis(pm(grdi,grdi,2)),basis(pm(grdi,grdi,3)),
     4               basis(pm(grdi,grdi,4)),basis(pm(grdi,grdi,1)),
     5               psi0(ppsi0(grdi)),t0,energy(grdi),spdim,
     6               ntot(grdi),n3d(grdi),nmax(2,grdi),nmax(3,grdi),
     7               nmax(4,grdi),nmax(1,grdi),timh0,i0stat)
c
 60   continue   
      call memory(-ngot(1),pvoft,idum,'voft',idum)   
      call memory(-ngot(2),pdrver,idum,'driver',idum)   
      call memory(-ngot(3),psoln,idum,'soln',idum)
      call memory(-ngot(4),pppsi0,idum,'psi0',idum)
      call memory(-ngot(5),pres,idum,'resid',idum)      
      call memory(-ngot(6),plufac,idum,'lufac',idum)
      call memory(-ngot(7),phamil,idum,'direct-solve',idum)
c      call initpsi(y(driver),y(psi0),basis(zpm(1)),n3d,nmax(1))
c
      return
 1    format(/,30x,'multigrid dvr code')
 2    format(/,5x,'problem may not be solved exactly on first grid')
 3    format(/,15x,'starting solution on grid = ',i3)
 4    format(/,5x,'interpolating initial wavefunction')
      end







