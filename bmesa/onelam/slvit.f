      subroutine slvit (gr1,gr2,varr,rhs,vguess,exvc,exiter,wts,aold,
     1                  anew,bold,bnew,ipvt,t1,msize,iter,
     2                  convg,ovtol,ops,nowgt,guess)
      implicit integer(a-z)
      logical logkey, nowgt
      real *8 gr1, gr2, varr, rhs, exvc, exiter, vguess
      real *8 wts, aold, anew, bold, bnew, t1, convg, ovtol
      dimension gr1(msize), gr2(msize)
      dimension varr(msize), rhs(msize), vguess(msize)
      dimension exvc(msize,iter), exiter(msize,iter)
      dimension wts(msize), aold(iter,iter), anew(iter,iter)
      dimension bold(iter), bnew(iter), ipvt(iter)
      dimension t1(msize)
      common /io/ inp, iout
      character *30 status
      character *(*) ops, guess
      character *3 itoc
      character *30 passin
*
      iwrit=0
      if (logkey(ops,'print=lam=convergence',.false.,' ')) iwrit=-1
      if (logkey(ops,'print=lam=iterations',.false.,' ')) iwrit=1
      if (logkey(ops,'non-orthogonal',.false.,' ')) then
          passin='generate vectors'
      else
           passin='generate orthogonal vectors'
      endif
      if (guess.eq.'rhs') then
          call copy(rhs,vguess,msize)
      endif
*
*
*
*
*     ----- the free wave right hand sides are already in rhs array ----

      slncnt=1
      bndcnt=0
      point=0
      iloop=1
      nsp=1
      call linslv ('prepare for solution',status,exvc,vguess,ovtol,
     1             convg,0,0,0,wts,iter,msize,iout,nsp,nvec,iwrit,
     2             nowgt)
   70 continue
      nold=nvec
      testno=nold+nsp
      if (testno.ge.iter) then
      call linslv ('final best solution',status,exvc,0,0,0,bnew,0,rhs,
     1             wts,0,msize,iter,iter,nvec,iwrit,nowgt)
      write (iout,130) status,nvec
      go to 90
      endif
      call linslv (passin,status,exvc,0,0,0,bnew,0,0,wts,ntrial,msize,
     1             iter,iter,nvec,iwrit,nowgt)
      if (status.eq.'hit iteration limit') then
      write (iout,130) status,nvec
      go to 90
      endif
      if (ntrial.eq.0) go to 80
*
*
*
      iold=nold+1
      call genvec (gr1,gr2,varr,exvc(1,iold),exiter(1,iold),t1,wts,
     1             msize)
*
*
      call linslv ('solve equations',status,exvc,exiter,anew,aold,
     1             bnew,bold,rhs,wts,ipvt,msize,iter,iter,nvec,iwrit,
     1             nowgt)
      if (status.eq.'continue iterations') go to 70
      if (status.eq.'starting to diverge') then
      write (iout,120) nold
      endif
      go to 90
   80 write (iout,150) nvec
   90 continue
*
*
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
