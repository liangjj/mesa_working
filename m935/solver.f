*deck @(#)solver.f	5.1  11/6/94
      subroutine solver(hss,b,c,t,tm,scrtch,c0,
     $                   hess,thess,ogl,gl,tgl,grad,td,
     $                   diag,rmsd,r,mdim,nrhs,nrow,incore,npvec,
     $                   mxiter,stopj,lnbuf,energy,
     $                   nbins,binsiz,nblks,blksiz,tsize,mxblk,lsize,
     $                   bptr,xbin,ibin,tvec,iu,freeze,npdim,filtyp)
c
c***begin prologue
c***date written       890124   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (llnl)
c***source             @(#)solver.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
c
c
      real*8 c0(mdim,npvec),hss(mdim,*),b(mdim,*),grad(mdim,nrhs)
      real*8 diag(mdim),ogl(*),gl(*),scrtch(*),rmsd(nrhs)
      real*8 t(mdim,*),c(mdim,*),thess(*),hess(*),tvec(*)
      real*8 td(mdim),tm(*),tgl(*),r(*),xbin(*)
      integer ibin(*),bptr(*)
      real*8 det(2) 
c
      integer binsiz,binsz2,blksiz,tsize,tblks
      integer freeze,twalks
      character*(*) filtyp
      logical debug
c
      parameter (debug=.false.)
      data eps/1.d-07/, mode/0/
      save eps,mode
c
      common /io/ inp,iout
c
c
      binsz2=binsiz+2
      maxpt=mxiter+nrhs+1
      twalks=mdim-npdim
c
      do 5 k=1,mdim
         td(k)=1.d0/diag(k)
  5   continue
c
c
      do 7 j=1,nrhs
         do 6 k=1,mdim
            b(k,j)=grad(k,j)*td(k)
 6       continue
 7    continue
c
c     ----- project the starting guess -----
c
      if(npvec.ne.0) then
         if(freeze.eq.0) then
            call  project(b,mdim,nrhs,c0,npvec,scrtch,t)
         else
            call zap(b,mdim,twalks,nrhs)
         end if
      end if
c
c     ----- normalize the vectors -----
c
      call schmit(b,mdim,0,nrhs,mmvec,scrtch,eps,ierr)
c
      if(debug) then
         write(iout,*)' solver:starting vectors '
         call matout(b,mdim,nrhs,mdim,nrhs,iout)
      endif
c
      mvec=nrhs
c
      itot=0
      iter=0
      jiter=0
c
c
 10   continue
c
      jiter=jiter+1
      if(debug) then
         write(iout,*)' iteration: ',jiter
      endif
c
      liter=iter
      is=iter+1
      iter=iter+mvec
      lastv=mvec
      nx=iter+1
c
c     iii-- project the gradient into the krylov space ---
c
      call ebtc(tgl,b(1,is),grad,mvec,mdim,nrhs)
      call movegrd(gl,ogl,tgl,liter,mvec,iter,nrhs)
c
      if(debug) then
         write(iout,*)'  current b-vector '
         call matout(b(1,is),mdim,nrhs,mdim,nrhs,iout)
      endif
c
c     ----- matrix multiplication -----
c
      call tmult(bptr,nbins,xbin,ibin,hss,binsiz,
     $           binsz2,tblks,nblks,lsize,blksiz,mxblk,
     $           iu,iprt,b(1,is),t,tvec,mdim,mvec,
     $           incore,energy,filtyp)
c
      if(debug) then
         call hmult(b(1,is),mdim,mvec,t,hss,nrow,incore)
         write(iout,*)' unprojected t-vector '
         call matout(t,mdim,nrhs,mdim,nrhs,iout)
      endif
c
c     ----- project the resultant vector -----
c
      if(npvec.ne.0) then
         if(freeze.eq.0) then
            call project(t,mdim,mvec,c0,npvec,scrtch,b(1,maxpt))
         else
            call zap(t,mdim,twalks,mvec)
         end if
      end if
c
      if (debug) then
         write(iout,*)' projected t-vector '
         call matout(t,mdim,nrhs,mdim,nrhs,iout)
      endif
c
c     -- get new correction vector --
c
c     -- old code        b(k,iter+j)=(t(k,j)-grad(k,j))*td(k)
c
c
      do 31 j=1,mvec
         do 30 k=1,mdim
            b(k,iter+j)=t(k,j)*td(k)
 30      continue
 31   continue
c
c     -- project the correction vector
c
      if(npvec.ne.0) then
         if(freeze.eq.0) then
            call project(b(1,nx),mdim,mvec,c0,npvec,
     $                   scrtch,b(1,maxpt))
         else
            call zap(b(1,nx),mdim,twalks,mvec)
         end if
      end if
c
      if(jiter.gt.1) then
         itot=liter*(liter+1)/2
         call scopy(itot,hess,1,thess,1)
      end if
c
c     form the small matrix
c
      call ebtc(r,b,t,iter,mdim,mvec)
      call makem(thess(itot+1),r,liter,mvec,iter)
c
      if(iter+nrhs.gt.mdim)go to 200
c
      call schmit(b,mdim,iter,mvec,mmvec,scrtch,eps,ierr)
      if(mmvec.ne.0) then
         lvec=mmvec
         call schmit(b,mdim,iter,lvec,mmvec,scrtch,eps,ierr)
         mvec=mmvec
      else
c        call lnkerr(' no new vectors ')
         ierr=1
      end if
c
c
      itot=iter*(iter+1)/2
c
c     -- save a triangular form of the small matrix
c
      call scopy(itot,thess,1,hess,1)
c
c     -- save the current projection of the rhs
c
      call scopy(iter*nrhs,gl,1,ogl,1)
c
      if(jiter.eq.1)go to 10
c
      icnvrg=0
c
c     -- square-up the small matrix
c
      call trtosq(tm,thess,iter,itot)
c
c     -- augment the small matrix with the rhs
c
      it=iter*iter+1
      call scopy(iter*nrhs,gl,1,tm(it),1)
c
c     -- solve the set of linear equations
c
      if(debug) then
         write(iout,*)'  tm '
         call matout(tm,iter,iter+nrhs,iter,iter+nrhs,iout)
      endif
c
      call rminv(iter,iter,tm,nrhs,tm(it),det,'all')
c
c     -- move the solutions to c(iter,nrhs)
c
      call scopy(iter*nrhs,tm(it),1,c,1)
c
c     -- test convergence
c
      call tcon(c,iter,nrhs,icnvrg,stopj)
c
      if(debug) then
         write(iout,*)' jiter ',jiter
         write(iout,*)' ierr icnvrg iter ',ierr,icnvrg,iter
      endif
c
      if(ierr.ne.0 .and. icnvrg.eq.0) go to 195
      if(icnvrg.eq.0.and.iter.lt.mxiter)go to 10
 195  continue
c
c     form final vectors
c
      call ebc(t,b,c,mdim,iter,nrhs)
c
      if(debug) then
         write(iout,*)' solution vectors '
         call matout(t,mdim,nrhs,mdim,nrhs,iout)
      endif
c
      call tmult(bptr,nbins,xbin,ibin,hss,binsiz,
     $             binsz2,tblks,nblks,lsize,blksiz,
     $             mxblk,iu,iprt,t,b(1,maxpt),tvec,
     $             mdim,mvec,incore,energy,filtyp)
c
c.      call hmult(t,mdim,nrhs,b(1,maxpt),hss,nrow,incore)
c
      if(npvec.ne.0) then
         if(freeze.eq.0) then
            call project(b(1,maxpt),mdim,nrhs,c0,npvec,scrtch,b)
         else
            call zap(b(1,maxpt),mdim,twalks,nrhs)
         end if
      end if
c
      if(debug) then
         write(iout,*)' approximate rhs '
         call matout(b(1,maxpt),mdim,nrhs,mdim,nrhs,iout)
      endif
c
      call tstcon(b(1,maxpt),grad,rmsd,mdim,nrhs,stopj,icnvrg)
c
      if(icnvrg.ne.1) then
         call lnkerr(' problem in linear equations ')
      endif
c
c
      return
c
c
 200  continue
      write(iout,210)
 210  format(/,'  stop  no solution to linear eqn.s')
      call lnkerr('no solution to linear equations ')
c
c
      end
