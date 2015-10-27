*deck @(#)mcaugh.f	1.2  7/30/91
      subroutine mcaugh(nmix,nf41,lenb,neig,thre,nviter,innorm,
     $     del, sqcdf, smlld, shftd,
     $     cr,icr,ncore,
     $     osqcdf,icutah,nmaxah,thrcah)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      880130   (yymmdd)
c   30 january 1988    bhl at brl
c      definition of iblk altered and routines davson and bigdav
c      corrected by changing nvec to nroot in most subroutine calls
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcaugh.f	1.2   7/30/91
c
c***purpose
c
c***description
c
c             this routine forms the augmented hessian matrix
c                     and diagonalizes it.
c
c --- input
c
c     nmix            no of orbital rotations
c     nf41            fortran no of the dataset containing the hessian
c                     and the gradient elements.
c     lenb            logical record length of nf41.
c     neig            desired number of eigenstates of the augmented
c                     hessian.  when neig .ne. 1, the eigenstate with
c                     the largest first element is selected.
c     thre            convergence threshold for eigenvalues
c     nviter          no of inverse iterations to be performed in
c                     determining the eigen vectors
c
c --- output
c
c     del(nmix)       the solution vector
c
c --- working storage
c
c     cr(ncor)
c
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8 (a-h,o-z)
      integer wpadti
      logical debug
cc
cmp   extended dummy cr,del,icr
cc
      dimension cr(2),icr(2)
      dimension del(2),thrv(20)
c
      data debug/.false./
      common /io/ inp,iout
c
      if (debug) then
         write (iout,9001) nf41
 9001    format(/'0************ augmented hessian section',
     $           ' ************',i8)
      end if
c
cc      call jtime(it1)
c-3     write(iout,*)'  icutah ',icutah
ctest
c       lenb=4*nmx
ctest
c
      nmx = nmix + 1
      if(nmix.gt.icutah)go to 301
c
      nah = nmx * (nmx + 1) / 2
c
c     allocate core.
      lah=1
      levl = lah + nah
      levc = levl + neig
      iscr = levc + neig * nmx
      last1 = iscr + 5*nmx
      lbuf = levl
      last2 = lbuf + lenb
c
      last=wpadti(max(last1,last2))
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(last,cr,ncore,'augmented hessian',0)
c
c-----------------------------
c     build augmented hessian
c-----------------------------
 199  continue
cc
c     if (debug) then
c        write(iout,9002) (cr(i+lah-1),i=0,nah)
c9002    format(/,'  mcaugh: lah'/4(2x,f14.8))
c     end if
c
      call rdtrah(cr(lah),cr(lbuf),lenb,nmix,nmx)
      call rdtgrd(cr(lah),cr(lbuf),lenb,nmix,nmx)
      call dampt(cr(lah),smlld,shftd,nmx)
cc
c---------------------------------------------------
c     diagonalize the augmented hessian matrix
c---------------------------------------------------
cc
cps      write(iout,*)'  givens-incore diagonalization '
c---------------------------
c     incore diagonalization
c---------------------------
      call givens(nmx,neig,nmx,cr(lah),cr(iscr),cr(levl),cr(levc))
c
      go to 3302
 301  continue
c
      if(icutah.lt.0)go to 302
c
c     allocate core.
      nah = nmx * nmx
      lah=1
      levl = lah + nah
      levc = levl + neig
      iscr = levc + neig * nmx
      last1 = iscr + 5*nmx
      lbuf = levl
      last2 = lbuf + lenb
c
      last=wpadti(max(last1,last2))
c
      nsml=neig+10
      msml=neig+3
      nroot=neig+1
      nlook=nmx
      thrn=1.d-6
      thrs=1.d-14
      nmax=nmaxah
c
      mdim=max(nsml,nmax*msml+msml)
      itr=(mdim*(mdim+1))/2
      isqr=mdim*mdim
      iblk=nmax*nmx*nroot+msml*nmx
c
      ib=iadtwp(last)
      ic=ib+iblk
      ish=ic+iblk+msml*nmx
      idiag=ish+itr
      itemp=idiag+nmx
      ieval=itemp+itr
      ievec=ieval+mdim
      ipt=ievec+isqr
      iscr=ipt+nsml
      ineed=wpadti(iscr+5*nmx)
c
      call getscm(0,cr,maxr,' big augmented hessian',0)
c
      if(ineed.gt.maxr) go to 302
c
c
      call getscm(last,cr,ncore,'augmented hessian',0)
cc
      call rdsqah(cr(lah),cr(lbuf),lenb,nmix,nmx)
      call rdsgrd(cr(lah),cr(lbuf),lenb,nmix,nmx)
      call damps(cr(lah),smlld,shftd,nmx)
cc
cc
c-3      write(iout,*)'  davidson-incore diagonalization  '
c-----
c     davidson diagonalization
c-----
c
      ncorr=ncore-last
c
      if(neig+2.gt.20) call lnkerr('mcaugh')
c
      if(osqcdf.gt.thrcah)go to 1303
      write(iout,1302) osqcdf,thrcah
 1302 format(/,20('*'),/,'  warning  osqcdf lt thrcah  ',2(2x,d15.7),/,
     $     '  decreasing thrcah  in davidson diagonalization of ',
     $     'aug. hes. ')
      thrcah=thrcah*1.d-3
 1303 continue
c
      do 303 ii=1,neig
         thrv(ii)=thrcah
 303  continue
      thrv(neig+1)=1.d-1
      thrv(neig+2)=1.d-1
c
      if(ncorr.lt.0) then
         call lnkerr(' ncorr lt 0  stop in mcaugh')
      endif
c
      call davsdr(cr(levc),cr(levc),ncorr,cr(lah),
     $     thrv,thrn,thrs,nmx,nmax,nroot,nsml,nlook,msml)
c
      go to 3302
c
 302  continue
c
      ld=1
      lah=ld+nmx
      levl=lah+lenb
      levc=levl+neig
      iscr=levc+neig*nmx
      last1=iscr+5*nmx
      lbuf=levl
      lg=lbuf+lenb
      last2=lg+nmx
c
      last=wpadti(max(last1,last2))
c
      nsml=neig+10
      msml=neig+3
      nroot=neig+1
      nlook=nmx
      thrn=1.d-6
      thrs=1.d-14
      nmax=nmaxah
c
      mdim=max(nsml,nmax*msml+msml)
      itr=(mdim*(mdim+1))/2
      isqr=mdim*mdim
      iblk=nmax*nmx*msml+msml*nmx
c
      ib=iadtwp(last)
      ic=ib+iblk
      ish=ic+iblk+msml*nmx
      idiag=ish+itr
      itemp=idiag+nmx
      ieval=itemp+itr
      ievec=ieval+mdim
      ipt=ievec+isqr
      iscr=ipt+nsml
      ineed=iscr+5*nmx
c
      call getscm(0,cr,maxr,' big augmented hessian',0)
c
      if(ineed.gt.maxr) then
         write(iout,*)'  ineed  maxr  ',ineed,maxr
         call lnkerr(' increase core for out-of-core augh ')
      endif
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(last,cr,ncore,'augmented hessian',0)
cc
cc
c-3      write(iout,*)'  davidson out-of-core diagonalization  '
c-----
c     davidson diagonalization
c-----
c
      ncorr=ncore-last
c
      if(neig+2.gt.20) call lnkerr(' mcaugh ')
c
      if(osqcdf.gt.thrcah)go to 2303
      write(iout,1302) osqcdf,thrcah
 2303 continue
c
      do 2304 ii=1,neig
         thrv(ii)=thrcah
 2304 continue
      thrv(neig+1)=1.d-1
      thrv(neig+2)=1.d-1
c
      if(ncorr.lt.0) then
         call lnkerr(' ncorr lt 0  stop in mcaugh')
      endif
c
      call formah(cr(lah),cr(lbuf),lenb,cr(lg),cr(ld),
     $            nmix,nmx,smlld,shftd)
      call bigdrv(cr(levc),cr(levc),ncorr,cr(lah),
     $            thrv,thrn,thrs,nmx,nmax,nroot,nsml,nlook,
     $            msml,cr(ld),lenb)
c
c
 3302 continue
c
      lv = levc - 1
c
      if (debug) then
         write (iout,9124) (cr(lv+i),i=1,nmx)
 9124    format('0*mcaugh eigenvectors'/4(1x,f16.8))
      end if
c-------------------------------------------
c     search for the desired eigenstate
c-------------------------------------------
      mvc = levc
      lvc = mvc
      t = 0.d+00
      do 600 i = 1, neig
         if (abs(cr(mvc)) .lt. t) go to 600
         lvc = mvc
         t = abs(cr(mvc))
 600  mvc = mvc + nmx
c--------------------------------------------------------------------
c     intermediate normalize the eigenvector and put correction
c        vector in del, with the correct phase
c--------------------------------------------------------------------
      t = cr(lvc)
      if (innorm .eq. 0) t = t / abs(t)
      sqcdf = 0.d0
      txt = 1.d0/t
      do 800 i = 1, nmix
         del(i) = cr(lvc+i) * txt
         sqcdf = sqcdf + del(i) * del(i)
 800  continue
c
      if (debug) then
         write (iout,9125) (del(i),i=1,nmix)
 9125    format(/,' ***mcaugh delta'/4(1x,f16.8))
      end if
c
c      call jtime(it2)
c      tim  = (it2 - it1) / 100.d0
c      write (iout,9099) tim
 9099 format('0   ****** tim  in this section ',10x,f10.2,' secs')
c
      return
      end
