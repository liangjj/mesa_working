      program solver
c note if npvec .ne. nbtot, q must be zero
      implicit real*8(a-h,o-z)
      parameter (nbig=800, nsmall=100,nchnl=20)
      complex*16  hdenom(nbig,nbig),smat(nsmall,nsmall),work(nbig)
      complex*16 tmat(nsmall,nsmall),zdotu,cov(nbig,nsmall)
      complex*16 htop(nbig,nsmall),htopp(nbig,nsmall),det
      real*8 kchan(nchnl),hpp(nbig,nbig),hqq(nbig,nbig),
     1    hpq(nbig,nsmall),hfix(nbig),xsecmat(nchnl,nchnl)
      integer nlm(nchnl),kpvt(nbig)
      equivalence (hqq,hpq,hpp),(hdenom,htop)
      character*8 istat,istatic
      data istatic/6hstatic/
      open(5,file='insolve')
      open(6,file='outsolve')
      open(7,file='bmesa',form='unformatted')
      open(9,file='bcoef',form='unformatted')
      open(77,file='btmat',form='unformatted')
      open(8,file='bsolve',form='unformatted')
c      call second(time)
      read(5,*)iprint,iopen,iestart
      if(iopen.eq.1)then
      read(5,*)dele,coef
      write(6,*)' energy defect and ci coefficient',dele,coef
      endif
      read(5,55)istat
 55   format(a8)
      if(istat.eq.istatic)then
      open(99,file='bstatic',form='unformatted')
      endif
      read(7)nchan,nbtot,nfree,(nlm(i),i=1,nchan)
      read(8)npvec,npdim
      write(6,*)' npdim ',npdim
      if(npvec.ne.nbtot)then
      write(6,*)' npvec and nbtot are unequal!', nbtot, npvec
      endif
      write(6,112)nchan,nbtot,nfree,(nlm(i),i=1,nchan)
 112  format(' number of channels=',i3/
     ^ ' dimension of P-space is ',i5/
     ^ ' number of free (lm) channels is =',i4/
     ^ ' number of lm terms per channel:'/(20i3))
      ntot=nbtot+nfree
      if(nchan.gt.nchnl.or.nfree.gt.nsmall.or.
     ^ ntot.gt.nbig)then
      write(6,*)' check dimensions'
      stop
      endif
      nall=ntot+npdim
      write(9)nchan,nall,nfree
      iw=nbtot+1
      read(7)nener
      do 1 ie=1,nener
      do 100 i=1,nbig
      do 101 j=1,nbig
 101     hdenom(i,j)=0.
      do 102 j=1,nsmall
 102  htopp(i,j)=0.
 100     continue
      read(7)(kchan(ic),ic=1,nchan)
      write(9)(kchan(ic),ic=1,nchan)
      write(6,111) (kchan(i),i=1,nchan)
111   format(///,' ********************************',
     # ' new energy ********************************',
     # //,' channel momenta:',6e12.5)
      if(nbtot.eq.npvec)then
      read(7)((hdenom(i,j),i=1,ntot),j=iw,ntot)
      do 114 i=1,ntot
      do 114 j=iw,ntot
 114  hdenom(j,i)=hdenom(i,j)   
      else
      read(7)((hdenom(i,j),i=1,ntot),j=1,ntot)
      endif
      call rdbinsqr(hpp,nbig,npvec,8)
      if(istat.eq.istatic)then
      call rdbinsqr(hpp,nbig,nbtot,99)
      endif
      if(nbtot.eq.npvec.or.istat.eq.istatic)then
      do 2 i=1,nbtot
      do 2 j=1,nbtot
 2       hdenom(i,j)=hpp(i,j)
      endif
      if(npdim.gt.0)then
      call rdbinsqr(hqq,nbig,npdim,8)
      do 3 i=1,npdim
      do 3 j=1,npdim
      ii=i+ntot
      jj=j+ntot
 3    hdenom(ii,jj)=hqq(i,j)
      if(iopen.eq.1)then
c      hfix(1)=hqq(1,1)+etarget+kchan(1)**2/2.
      hfix(1)=hqq(1,1)+dele
      if(npdim.gt.1)then
      do 7 i=2,npdim
 7       hfix(i)=hqq(1,i)
      endif
      endif
      do 4 i=1,npvec
      read(8)(hpq(j,i),j=1,npdim)
      do 5 j=1,npdim
      hdenom(i,j+ntot)=hpq(j,i)
 5    hdenom(j+ntot,i)=hpq(j,i)
 4    continue
      endif
      if(iprint.ne.0) then
      write(6,107)
107   format(//' denominator matrix of (h-e)')
      do 60 i=1,nall
60    write(6,108) i,(hdenom(j,i),j=1,nall)
108   format(1x,i3,6e12.5,/,(4x,6e12.5))
      endif
      read(7)((htopp(i,j),i=1,ntot),j=1,nfree)
      read(7)nread
      write(6,*)' no. occupied plus scattering terms', nread
      read(7)((cov(i,j),i=1,nread),j=1,nfree)
      do 6 k=1,nbtot
      do 6 i=1,npdim
      do 6 j=1,nfree
 6       htopp(i+ntot,j)=htopp(i+ntot,j)+hpq(i,k)*cov(k,j)
c
c add extra terms to RHS for correlated target case
c
      if(iopen.eq.1)then
      do 8 i=1,npdim
      do 8 j=1,nfree
 8       htopp(i+ntot,j)=htopp(i+ntot,j)+coef*hfix(i)
     # *cov(nread,j)
c**********new**********
      do 88 i=1,nbtot
      do 88 j=1,nfree
 88      htopp(i,j)=htopp(i,j)+coef*hpq(1,i)*cov(nread,j)
c**********************
      endif
      if(iprint.ne.0) then
      write(6,109)
109   format(//' numerator matrix of (h-e)')
      do 61 i=1,nfree
61    write(6,108) i,(htopp(j,i),j=1,ntot)
      endif
c      call second(tnow)
      elapse=tnow-time
      write(6,777)elapse
777   format(" elapsed time for setup is ",e12.4)
      time=tnow      
      if(ie.ge.iestart)then
c      call matinv(hdenom,nall,htopp,nfree,det,nbig)
c      call csifa(hdenom,nbig,nall,kpvt,info)
      call cgefs(hdenom,nbig,nall,htopp(1,1),1,ind,work,kpvt)
c      do 666 i=1,nfree
      do 666 i=2,nfree
c 666     call csisl(hdenom,nbig,nall,kpvt,htopp(1,i))
 666  call cgefs(hdenom,nbig,nall,htopp(1,i),2,ind,work,kpvt)
      endif
      if(iprint.ne.0) then
      write(6,187)
187   format(//' solution matrix ')
      do 79 i=1,nall
79    write(6,108) i,(htopp(i,j),j=1,nfree)
      endif
      write(9)((htopp(i,j),i=1,nall),j=1,nfree)
c      call second(tnow)
      elapse=tnow-time
      write(6,778)elapse
      time=tnow
 778  format(" time for matinv is ",e12.4)
      do 103 i=1,nbig
      do 103 j=1,nsmall
 103  htop(i,j)=0.
      read(7)((htop(i,j),i=1,ntot),j=1,nfree)
      if(iprint.ne.0) then
      write(6,110)
110   format(//' conjugate numerator matrix of (h-e)')
      do 62 i=1,nfree
62    write(6,108) i,(htop(j,i),j=1,ntot)
      endif
      do 500 ilm=1,nfree
      do 500 jlm=1,nfree
      tmat(ilm,jlm)=zdotu(nall,htopp(1,jlm),1,htop(1,ilm),1)
500   continue
      read(7)((smat(i,j),i=1,nfree),j=1,nfree)
      if(ie.lt.iestart)go to 1
      if(iprint.ne.0) then
      write(6,113)
 113  format(//' first Born matrix')
      do 63 i=1,nfree
63    write(6,108) i,(smat(j,i),j=1,nfree)
      endif
      do 501 ilm=1,nfree
      do 501 jlm=1,nfree
      tmat(ilm,jlm)=-2.*(-tmat(ilm,jlm)+smat(ilm,jlm))
 501  continue
      write(6,115)
 115  format(//' T-Matrix')
      do 64 i=1,nfree
64    write(6,108) i,(tmat(j,i),j=1,nfree)
      istart = 0
      do 879 ic=1,nchan
      jstart = 0
      ni = nlm(ic)
      do 878 jc=1,nchan
      nj=nlm(jc)
      ist=istart+1
      ifin=istart+ni
      jst=jstart+1
      jfin=jstart+nj
      write(77)ic,jc,ni,nj,kchan(ic),kchan(jc)
c      write(88,177)ic,jc,ni,nj,kchan(ic),kchan(jc)
 177  format(4i5,2f20.10)
      write(77)((tmat(ii,jj),ii=ist,ifin),jj=jst,jfin)
c      write(88,277)((tmat(ii,jj),ii=ist,ifin),jj=jst,jfin)
 277  format(4e20.10)
      summod = 0.
      do 877 ilm=1,ni
      do 877 jlm=1,nj
      isub = istart + ilm
      jsub = jstart + jlm
877   summod = summod + abs(tmat(isub,jsub))**2
      xsecmat(ic,jc) = 4.0*3.141592654*summod/kchan(ic)**2
      jstart = jstart + nlm(jc)
878   continue
      istart = istart + nlm(ic)
879   continue
       write(6,874)
874   format(//,' total cross sections: row index = initial chnl,',
     # ' column index = final chnl')
      write(6,871) (i,i=1,nchan)
871   format(4x,6(6x,i2,4x))
      do 873 i=1,nchan
      write(6,872) i, (xsecmat(i,j),j=1,nchan)
873   continue
872   format(1x,i2,1x,6e12.5/(4x,6e12.5))
c      call second(tnow)
      elapse=tnow-time
      write(6,779)elapse
      time=tnow
 779  format("  time at end is ",e12.4)
 1    continue
      call exit
      end
