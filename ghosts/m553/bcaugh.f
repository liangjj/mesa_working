      subroutine bcaugh(nmix,nf41,lenb,neig,thre,nviter,innorm,
     $     del, sqcdf, smlld, shftd,
     $     cr,icr,ncore,
     $     osqcdf,icutah,nmaxah,thrcah)
C
C***Begin prologue     bcaugh
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
c
c-----------------------------------------------------------------------
c
c --- description     this routine forms the augmented hessian matrix
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
c-----------------------------------------------------------------------
c
C
C***References
C
C***Routines called    (none)
C
C***End prologue       bcaugh
C
      implicit real*8 (a-h,o-z)
c
      dimension cr(2),icr(2)
      dimension del(2),iwrk1(16),thrv(20)
c
      common /io/ inp,iout
c
cps      write (iout,9001) nf41
cps 9001 format(/'0************ augmented hessian section',
cps     1      ' ************',i8)
c
cc      call jtime(it1)
c
c-3      write(iout,*)'  ICUTAH ',icutah
c
      nmx = nmix + 1
c
      if(nmix.gt.icutah)go to 301
c
      nah = nmx * (nmx + 1) / 2
      lah=1
      levl = lah + nah
      levc = levl + neig
      iscr = levc + neig * nmx
      last1 = iscr + 5*nmx
      lbuf = levl
      last2 = lbuf + lenb
c
      last=max(last1,last2)
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
      call rdtrah(cr(lah),cr(lbuf),lenb,nmix,nmx)
c.2      write(iout,*)' tri. hmo '
c.2      call printm(cr(lah),5,1)
      call rdtgrd(cr(lah),cr(lbuf),lenb,nmix,nmx)
c.2      write(iout,*)' tri. hmo + grad '
c.2      call printm(cr(lah),5,1)
      call dampt(cr(lah),smalld,shftd,nmx)
cc
c---------------------------------------------------
c     diagonalize the augmented hessian matrix
c---------------------------------------------------
cc
      write(iout,*)'  incore diagonalization '
c---------------------------
c     incore diagonalization
c---------------------------
      call givens(nmx,neig,nmx,cr(lah),cr(iscr),cr(levl),cr(levc))
c
      go to 302
 301  continue
c
      nah = nmx * nmx
      lah=1
      levl = lah + nah
      levc = levl + neig
      iscr = levc + neig * nmx
      last1 = iscr + 5*nmx
      lbuf = levl
      last2 = lbuf + lenb
c
      last=max(last1,last2)
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(last,cr,ncore,'augmented hessian',0)
cc
      call rdsqah(cr(lah),cr(lbuf),lenb,nmix,nmx)
      call rdsgrd(cr(lah),cr(lbuf),lenb,nmix,nmx)
c.2      write(iout,*)' square augmented hessian '
c.2      call printm(cr(lah),nmx,2)
      call damps(cr(lah),smalld,shftd,nmx)
cc
      write(iout,*)'  davidson diagonalization  '
c-----
c     davidson diagonalization
c-----
c
      ncorr=ncore-last
c
      nsml=neig+10
      msml=neig+3
      nroot=neig+1
      nlook=nmx
      thrn=1.d-6
      thrs=1.d-14
      nmax=nmaxah
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
      call davsdr(cr(levc),icr(levc),ncorr,cr(lah),
     $     thrv,thrn,thrs,nmx,nmax,nroot,nsml,nlook,msml)
c
 302  continue
c
      lv = levc - 1
c     write (iout,9124) (cr(lv+i),i=1,nmx)
c9124 format('0*mcaugh eigenvectors'/4(1x,f16.8))
c-------------------------------------------
c     search for the desired eigenstate
c-------------------------------------------
      mvc = levc
      lvc = mvc
      t = 0.d0
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
c      write (iout,9125) (del(i),i=1,nmix)
c 9125 format(/,' ***mcaugh delta'/4(1x,f16.8))
c
c      call jtime(it2)
c      tim  = (it2 - it1) / 100.d0
c      write (iout,9099) tim
 9099 format('0   ****** tim  in this section ',10x,f10.2,' secs')
c
      return
      end
