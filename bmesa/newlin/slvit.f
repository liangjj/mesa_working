      subroutine slvit (gr1,gr2,fopto,flago,foptn,flagn,varr,rhs,exvc
     1 ,exiter,vopt,wts,aold,anew,bold,bnew,ipvt,t1,t2,nbk,nbkpls,brk,c
     2 ,sc,spln,nopt,nlag,jind,nsts,ncst,ncmax,msize,ntchn,nptmx,iter
     3 ,maxvec,nfopt,nflag,nsol,numdo,ldim,nds1,nds2,nds3,nds4,
     4 convg,ovtol,ien,ops,soldir,nowgt)
      implicit integer(a-z)
      logical logkey, nowgt
      real *8gr1, gr2, fopto, flago, foptn, flagn, varr, rhs, exvc
      real *8vopt, wts, aold, anew, bold, bnew, t1, t2, convg, ovtol
      real *8exiter, brk, c, sc, spln
      dimension gr1(nptmx,ntchn), gr2(nptmx,ntchn)
      dimension fopto(nptmx,ldim,nfopt), flago(nptmx,ldim,nflag)
      dimension foptn(nptmx,ldim,nfopt), flagn(nptmx,ldim,nflag)
      dimension varr(ntchn,ntchn,nptmx), rhs(msize,nsol)
      dimension exvc(msize,iter), exiter(msize,iter), vopt(nfopt,nfopt)
      dimension wts(msize), aold(iter,iter), anew(iter,iter)
      dimension bold(iter,nsol), bnew(iter,nsol), ipvt(iter)
      dimension t1(msize,nsol), t2(nfopt,nsol)
      dimension nopt(nsts), nlag(nsts)
      dimension jind(ncmax,nsts), ncst(nsts)
      dimension brk(nds1), c(nds2), sc(nds3)
      dimension spln(nds4)
      common /io/ inp, iout
      character *30 status
      character *(*) ops, soldir
      character *80 wrtv, wrtspn
      character *3 itoc
      wrtv(n1)='write real "solutions'//itoc(n1)//'" to wavefn'
      wrtspn(n1)='write real "spline vecs'//itoc(n1)//'" to wavefn witho
     1ut rewinding'
*
*     ----- note that nbk,nbkpls,brk,c and sc are equivalenced to spln -
*     -----      in the calling sequence to this routine -----
*     ----- this is necessary in order to make write statement -----
*     -----        contiguous in memory -----
*
      ilen=nsol*ntchn*1203
*
*
*
      iwrit=0
      if (logkey(ops,'print=lam=convergence',.false.,' ')) iwrit=-1
      if (logkey(ops,'print=lam=iterations',.false.,' ')) iwrit=1
*
*
*
*
*     ----- the free wave right hand sides are already in rhs array ----
*     ----- if there are orbitals needed for orthogonalization -----
*     ----- put them in rhs from ntchn+1 to nsol -----
      if (nsol.gt.ntchn) then
      do 20 isol=ntchn+1,nsol
      do 10 i=1,msize
   10 rhs(i,isol)=0.d+00
   20 continue
      slncnt=ntchn
      bndcnt=0
      point=0
      do 60 is=1,nsts
      nl=nlag(is)
      if (nl.eq.0) go to 60
      nsch=ncst(is)
      do 50 l=1,nl
      slncnt=slncnt+1
      bndcnt=bndcnt+1
      pcnt=point
      do 40 nc=1,nsch
      lp=jind(nc,is)
      do 30 k=1,nptmx
      pcnt=pcnt+1
      rhs(pcnt,slncnt)=flagn(k,lp,bndcnt)
   30 continue
   40 continue
   50 continue
   60 point=point+nptmx*ncst(is)
      if (slncnt.ne.nsol) then
      write (iout,140) slncnt
      stop 'slncnt'
      endif
      endif
*
*
      if (soldir.ne.'all') then
      iloop=numdo
      nsp=1
      else
      iloop=1
      nsp=nsol
      endif
      do 90 loop=1,iloop
      vecpnt=1
      if (soldir.ne.'all') vecpnt=loop
      call linslv ('prepare for solution',status,exvc,rhs(1,vecpnt)
     1 ,ovtol,convg,0,0,0,wts,iter,msize,iout,nsp,nvec,iwrit,nowgt)
   70 continue
      nold=nvec
      testno=nold+nsp
      if (testno.ge.iter) then
      call linslv ('final best solution',status,exvc,0,0,0,bnew,0,rhs(1
     1 ,vecpnt),wts,0,msize,iter,iter,nvec,iwrit,nowgt)
      write (iout,130) status,nvec
      go to 90
      endif
      call linslv ('generate vectors',status,exvc,0,0,0,bnew,0,0,wts
     1 ,ntrial,msize,iter,iter,nvec,iwrit,nowgt)
      if (status.eq.'hit iteration limit') then
      write (iout,130) status,nvec
      go to 90
      endif
      if (ntrial.eq.0) go to 80
*
*
*
      iold=nold+1
      call genvec (gr1,gr2,fopto,foptn,varr,exvc(1,iold),exiter(1,iold)
     1 ,vopt,t1,t2,wts,ncst,jind,nopt,nsts,ncmax,ldim,msize,ntchn,nptmx
     2 ,nfopt,ntrial)
*
*
      call linslv ('solve equations',status,exvc,exiter,anew,aold,bnew
     1 ,bold,rhs(1,vecpnt),wts,ipvt,msize,iter,iter,nvec,iwrit,
     1     nowgt)
      if (status.eq.'continue iterations') go to 70
      if (status.eq.'starting to diverge') then
      write (iout,120) nold
      endif
      go to 90
   80 write (iout,150) nvec
   90 continue
*
*
*          write solution vectors to channel
*
      call iosys (wrtv(ien),nsol*msize,rhs,0,0)
      if (nflag.ne.0) return
      if (logkey(ops,'lam=output=wavefunctions',.false.,' ')) then
      call iosys ('create real "spline vecs'//itoc(ien)//'" on wavefn'
     1 ,ilen,0,0,0)
      do 110 isol=1,nsol
      icnt=0
      do 100 ii=1,nsts
      nsch=ncst(ii)
      do 100 jj=1,nsch
      call bsplin (nptmx,rgs,rhs(icnt+1,isol),3,nbk,brk,c,iflag,sc)
      call iosys (wrtspn(ien),1203,spln,0,0)
  100 icnt=icnt+nptmx
  110 continue
      endif
      return
c
  120 format (/,5x,'solution starting to diverge. will quit at',1x,i4,1x
     1 ,'vectors')
  130 format (/,5x,'warning msg from solver',2x,a40,//,5x,'no. vectors',
     1 ' used',1x,i5)
  140 format (/,5x,'error in solution count.',1x,i4)
  150 format (//,5x,'***** no more vectors. finished *****',//,10x,'no.
     1vectors used',1x,i5)
      end
