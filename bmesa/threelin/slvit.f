c $Header: slvit.f,v 1.2 92/12/12 09:35:09 bis Exp $
*deck slvit
      subroutine slvit (diag,sudiag,spdiag,v,f,rhs,x,s4,rdel,energy,
     1                  refe,temp,guess,exvc,exiter,aold,anew,bold,
     2                  bnew,ipvt,points,nfd,iter,convg,ovtol,gtype,
     3                  slndir,bcond,urefe,ops)
      implicit integer(a-z)
      logical logkey, nowgt, urefe
      real*8  diag, sudiag, spdiag, v, f, rhs, x, s4, exvc, exiter
      real*8 rdel, energy, temp, guess, aold, anew, bold, bnew, convg
      real*8 value, ovtol, wts, refe
      dimension diag(0:points), sudiag(0:points), spdiag(0:points)
      dimension v(0:points), f(0:points), rhs(0:points), x(0:points)
      dimension s4(4), exvc(points,iter), exiter(points,iter)
      dimension guess(0:points), aold(iter,iter), anew(iter,iter)
      dimension bold(iter), bnew(iter), ipvt(iter), temp(0:points,6)
      common /io/ inp, iout
      character *30 status
      character *(*) ops, gtype, slndir, bcond
*
      nowgt=.true.
      iwrit=0
      if (logkey(ops,'print=convergence',.false.,' ')) iwrit=-1
      if (logkey(ops,'print=iterations',.false.,' ')) iwrit=1
*
*
*
*     ----- initialize routine and pass in the guess vectors ----
      status=slndir
      nsp=1
      call linslv ('prepare for solution',status,exvc,guess(1),
     1              ovtol,convg,0,0,0,wts,iter,nfd,iout,nsp,nvec,
     2              iwrit,nowgt,points)
      maxtrp=iter/nsp+1
      do 20 trips=1,maxtrp
         nold=nvec
         testno=nold+nsp
         if (testno.lt.iter) then
             call linslv ('generate vectors',status,exvc,
     1                     0,0,0,bnew,0,0,wts,ntrial,nfd,iter,
     2                     iter,nvec,iwrit,nowgt,points)
             if (status.eq.'hit iteration limit') then
                 write (iout,130) status,nvec
                 return 
             endif
             if (ntrial.ne.0) then
                 iold=nold+1
                 if (slndir.eq.'one-plus-matrix'.or.
     1               slndir.eq.'one-minus-matrix') then
                     call genvpt(diag,sudiag,spdiag,v,f,x,s4,rdel,
     1                           energy,refe,temp,exvc(1,iold),
     2                           exiter(1,iold),value,
     3                           bcond,points,nfd,urefe)
                 else
                     call genvec (diag,sudiag,spdiag,
     1                            exvc(1,iold),exiter(1,iold),
     2                            points,nfd,gtype)
                 endif
                 call linslv ('solve equations',status,exvc,
     1                         exiter,anew,aold,bnew,bold,rhs(1),
     2                         wts,ipvt,nfd,iter,iter,nvec,iwrit,
     3                         nowgt,points)
                 if (status.eq.'starting to diverge') then
                     write (iout,120) nold
                     return
                 endif
             else
                 write (iout,150) nvec
                 return 
             endif
         else
             call linslv ('final best solution',status,exvc,0,0,0,
     1                     bnew,0,rhs(1),wts,0,nfd,iter,iter,nvec,
     2                     iwrit,nowgt,points)
             write (iout,130) status,nvec
             return
         endif
   20 continue
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

