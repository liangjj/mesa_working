*deck fullmg.f 
c***begin prologue     fullmg
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, orthogonal polynomial
c***                   bose-einstein
c***author             schneider, b. i.(nsf)
c***source             tprop
c***purpose            multigrid solution of the time dependent 
c***                   schroedinger equation. 
c***                   
c***description        
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       fullmg
      subroutine fullmg(pm,ham,trn,eig,pvt,pdrv,prhs,psln,
     1                  pbasis,parray,pvoft,pdrver,psoln,pprhs,
     2                  phamil,plufac,energy,spdim,dim,nmax,n3d,
     3                  ntot,m,topg,mgr,cnverg,thresh,eps,
     4                  maxit,maxvec,nblck,mcycle,precon,
     5                  timh0,i0stat,prnt,mprnt,ngot)
c
      implicit integer (a-z)
      real*8 basis, array, driver, voft, rhs, soln, hamil, lufac
      real*8 t0
      real*8 thresh, cnverg, eps
      real*8 energy, scale
      character*(*) timh0, i0stat, precon
      character*80 title
      character*3 itoc
      logical prnt, mprnt
      dimension pm(mgr,mgr,dim)    
      dimension ham(mgr,dim), trn(mgr,dim), eig(mgr,dim)
      dimension pvt(mgr), pdrv(mgr), psln(mgr), prhs(mgr)
      dimension nmax(4,mgr), maxit(mgr), maxvec(mgr)
      dimension n3d(mgr), ntot(mgr), m(mgr), ngot(*)
      dimension energy(mgr)
      dimension prnt(*), mprnt(10)
      common/io/inp, iout      
      pointer(pbasis,basis(1))
      pointer(parray,array(1))
      pointer(pvoft,voft(1))     
      pointer(pdrver,driver(1))
      pointer(pprhs,rhs(1))
      pointer(psoln,soln(1))
      pointer(phamil,hamil(1)), (phamil,ihamil(1))
      pointer(plufac,lufac(1)), (plufac,indlu(1))
c
c     begin calculation
c
      scale=1.d0    
      do 10 grdi=2,topg
         write(iout,3) grdi
c
c        get starting solution by interpolating from previously 
c        calculated solution.
c
         call interp(pm,psln,pbasis,psoln,nmax,ntot,m,grdi,grdi-1,
     1               mgr,2,dim)
         if(mprnt(2)) then
            title='interpolated initial guess on grid = '//itoc(grdi)
            call prntrm(title,soln(soli),m(grdi),1,m(grdi),1,iout)
         endif
         title='entering main davidson routine on grid = '//itoc(grdi)
c
c        the output of the call to qcklin is a vector which is an
c        approximate solution to the equation;
c
c                     A|soln> = |Driver>
c 
         call qcklin(psoln,pdrver,pbasis,parray,pvoft,
     1               phamil,plufac,pm,psln,ham,trn,eig,pvt,pdrv,
     2               scale,cnverg,thresh,eps,maxit,maxvec,dim,
     3               spdim,nmax,m,ntot,n3d,mgr,grdi,precon,nblck,
     4               title,prnt(9),mprnt)
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
      call memory(-ngot(1),psoln,idum,'soln',idum)
      call memory(-ngot(2),pvoft,idum,'voft',idum)   
      call memory(-ngot(3),pdrver,idum,'driver',idum)   
      call memory(-ngot(4),pprhs,idum,'prhs',idum)
      call memory(-ngot(5),plufac,idum,'lufac',idum)
      call memory(-ngot(6),phamil,idum,'direct-solve',idum)
      return
 1    format(/,30x,'multigrid dvr code')
 2    format(/,5x,'problem may not be solved exactly on first grid')
 3    format(/,15x,'starting solution on grid = ',i3)
 4    format(/,5x,'interpolating initial wavefunction')
      end







