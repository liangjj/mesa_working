      subroutine genvec (gr1,gr2,fopto,foptn,varr,veco,vecn,vopt,t1,t2
     1 ,wts,ncst,jind,nopt,nsts,ncmax,ldim,msize,ntchn,nptmx,nfopt
     2 ,ntrial)
      implicit integer(a-z)
      real *8gr1, gr2, fopto, foptn, varr, veco, vecn, vopt, t1, t2, wts
      real *8sumf, sumb
      dimension varr(ntchn,ntchn,nptmx)
      dimension veco(msize,ntrial), vecn(msize,ntrial), wts(msize)
      dimension gr1(nptmx,ntchn), gr2(nptmx,ntchn)
      dimension fopto(nptmx,ldim,nfopt), foptn(nptmx,ldim,nfopt)
      dimension vopt(nfopt,nfopt)
      dimension t1(msize,ntrial), t2(nfopt,ntrial)
      dimension ncst(nsts), nopt(nsts)
      dimension jind(ncmax,nsts)
*
*          subroutine to generate a new iteration in the vector
*          sequence to solve the integral equations.
*          the method used is a direct technique in which the
*          matrix is not actually formed. its effect on the old vector
*          is computed directly from the regular and
*          irregular solutions , the potential matrix and the
*          previous vector in the chain.
*
*
*
      do 10 itrl=1,ntrial
      do 10 i=1,msize
      vecn(i,itrl)=0.d+00
      t1(i,itrl)=0.d+00
   10 continue
*
*
*
*         calculate sum over channels at point i of the potential
*         multiplied by the old vector.
*
      do 20 itrl=1,ntrial
      do 20 i=1,nptmx
      call sgemv (ntchn,ntchn,varr(1,1,i),ntchn,veco(i,itrl),nptmx,t1(i
     1 ,itrl),nptmx,1)
   20 continue
*
*          now begin main loops.
*
*         do forward recursion
*
      do 50 itrl=1,ntrial
      vecnt=0
      do 40 i=1,ntchn
      sumf=0.d+00
      do 30 j=1,nptmx
      vecnt=vecnt+1
      sumf=sumf+gr1(j,i)*wts(vecnt)*t1(vecnt,itrl)
      vecn(vecnt,itrl)=vecn(vecnt,itrl)+gr2(j,i)*sumf
   30 continue
   40 continue
   50 continue
*
*         do backward recursion
*
      ptmnus=nptmx-1
      do 80 itrl=1,ntrial
      vecnt=0
      do 70 i=1,ntchn
      sumb=0.d+00
      do 60 j=ptmnus,1,-1
      jpt=j+1
      jvcnt=vecnt+jpt
      jvcntm=jvcnt-1
      sumb=sumb+gr2(jpt,i)*wts(jvcnt)*t1(jvcnt,itrl)
      vecn(jvcntm,itrl)=vecn(jvcntm,itrl)+gr1(j,i)*sumb
   60 continue
   70 vecnt=vecnt+nptmx
   80 continue
*
*         add in optical potential contribution if needed
      if (nfopt.eq.0) return
*
*          form overlap matrix of current iterate and old
*                         orbitals
*
      bndcnt=0
      point=0
      do 130 is=1,nsts
      nmoi=nopt(is)
      nsch=ncst(is)
      if (nmoi.eq.0) go to 130
      do 120 j=1,nmoi
      bndcnt=bndcnt+1
      do 110 itrl=1,ntrial
      t1(bndcnt,itrl)=0.d+00
      pcnt=point
      do 100 nc=1,nsch
      lp=jind(nc,is)
      do 90 k=1,nptmx
      pcnt=pcnt+1
      t1(bndcnt,itrl)=t1(bndcnt,itrl)+fopto(k,lp,bndcnt)*wts(pcnt)*veco
     1 (pcnt,itrl)
   90 continue
  100 continue
  110 continue
  120 continue
  130 point=point+nsch*nptmx
*           now muliply the optical potential matrix
*           with the overlap vector
*
      call sgemm (nfopt,nfopt,ntrial,vopt,nfopt,t1,msize,t2,nfopt,0,1)
      kcnt=0
      bndcnt=0
      do 180 is=1,nsts
      nmoi=nopt(is)
      nsch=ncst(is)
      if (nmoi.eq.0) go to 180
      do 170 nc=1,nsch
      lp=jind(nc,is)
      do 160 lpt=1,nptmx
      kcnt=kcnt+1
      orbcnt=bndcnt
      do 150 nl=1,nmoi
      orbcnt=orbcnt+1
      do 140 itrl=1,ntrial
      vecn(kcnt,itrl)=vecn(kcnt,itrl)-foptn(lpt,lp,orbcnt)*t2(orbcnt
     1 ,itrl)
  140 continue
  150 continue
  160 continue
  170 continue
  180 bndcnt=bndcnt+nmoi
      return
      end
