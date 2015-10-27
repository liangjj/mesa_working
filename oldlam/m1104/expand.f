      subroutine expand (ops,nsts,ntchn,nptmx,neng,nmotot,nfopt,nflag
     1 ,iter,maxvec,maxlex,totlex,lammax,ncmax,lplsmx,uperl,orbctl
     2 ,optctl,charge,rmtrad,convg,ovtol,nsol,numdo,msize,memmax,vloc,
     3 orbs,cntpsi,optfle,soldir,nowgt)
*
      implicit integer(a-z)
      common z(1)
      common /io/ inp, iout
      dimension a(1)
      dimension tim(3)
      real *8z, charge, rmtrad, convg, ovtol, tim
      equivalence (a,z)
      character *3 orbctl, optctl
      character *8 vloc, orbs, cntpsi, optfle
      character *(*) ops, soldir
      logical nowgt
      data maxcor /1/
*
*
*     ----- outlay memory for channel and target arrays, energies etc --
*
      if (orbctl.ne.'orb') then
      nmotot=0
      nfopt=0
      nflag=0
      maxlex=0
      totlex=0
      lammax=0
      uperl=0
      endif
c     allocate space for integer arrays
      bndm=1
      bndspn=bndm+nsts
      bndpar=bndspn+nsts
      enbnd=bndpar+nsts
      ncst=enbnd+nsts
      cntm=ncst+nsts
      cntspn=cntm+nsts
      cntpar=cntspn+nsts
      l0cnt=cntpar+nsts
      chnloc=l0cnt+nsts
      jind=chnloc+nsts*lplsmx
c     end integer allocation begin real
      eng=jind+nsts*ncmax
      en=eng+neng
      rk=en+nsts*ncmax
      pt=rk+nsts*ncmax
      wt=pt+nptmx
c     end real begin integer
      nopt=wt+nptmx
      lstopt=nopt+nsts
      kepvnl=lstopt+nfopt*nsts
      kpdiag=kepvnl+nfopt*nsts
      nlag=kpdiag+nsts*nsts
      lstlag=nlag+nsts
c     total number optical and laguerre functions
      nftot=nfopt+nflag
*     ----- memory for local potential -----
*
c     some of this storage is also used for temporaries in setval
c     this is not a problem since we outlay locpot as real
      locpot=lstlag+nflag*nsts
      indxc=locpot
      lch=indxc+ntchn
      ech=lch+ntchn
*     ----- memory for greens functions -----
      gr1=locpot+ntchn*ntchn*nptmx
      gr2=gr1+ntchn*nptmx
*
*     ----- memory for weight vector -----
      wts=gr2+ntchn*nptmx
*     ----- memory for right hand sides -----
      rhs=wts+msize
*     ----- memory for the optical orbitals  -----
      orbtls=rhs+nsol*msize
*     ----- memory for laguerre orbitals -----
      orblag=orbtls+nfopt*uperl*nptmx
*     ----- memory for corresponding inhomogeneities -----
      inhopt=orblag+nflag*uperl*nptmx
      inhlag=inhopt+nfopt*uperl*nptmx
*
*     ----- memory for the solution of the la equations -----
*
*     ----- for the expansion vectors -----
      exvec=inhlag+nflag*uperl*nptmx
*     ----- for the iterates -----
      itvec=exvec+iter*msize
*     ----- for the two small sets of linear equations and the -----
*     ----- corresponding right hand sides -----
      smleq1=itvec+iter*msize
      smleq2=smleq1+iter*iter
      smrhs1=smleq2+iter*iter
      smrhs2=smrhs1+nsol*iter
*     ----- for two scratch arrays -----
*     ----- and a small vector -----
      scr1=smrhs2+nsol*iter
      scr2=scr1+msize*nsol
      iscrch=scr2+nfopt*nsol
*     ----- for the spline arrays needed to fit wavefunction -----
      spln=iscrch+iter
      nbk=spln
      nbkpls=nbk+1
      brk=nbkpls+1
      c=brk+301
      sc=c+900
c     we use these as real arrays in subroutine
      nds1=301
      nds2=900
      nds3=3605
      nds4=4806
c     we allocate 9614 integer words for spln. two integers plus the sum
c     of nds1----nds3
c     nds4=9614
c     nds5=2404
*     ----- memory for the optical potential matrix -----
*     ----- and a scratch array -----
      optmat=sc+3605
      scr3=optmat+nfopt*nfopt
*
*     ----- memory for the overlap matrix and the r matrix -----
*     ----- since they are not used at same time storage is identical --
      ovlpa=scr3+nfopt
      rmat=ovlpa
      nwords=rmat+nsol*max(nflag,nsol)
*     ----- find out how much core available -----
      call getscm (0,z,canget,'m1104: how much core',0)
      if (memmax.eq.0) then
      if (canget.lt.nwords) go to 20
      top=min(nwords,canget)
      call getscm (top,z,maxcor,'m1104',0)
      else
      if (memmax.gt.canget) go to 20
      if (nwords.gt.memmax) go to 20
      call getscm (nwords,z,maxcor,'m1104',0)
      endif
      write (iout,90) nwords
*     ----- set up the various channel arrays, energies and point -----
*     -----           information needed in calculation           -----
      call setval (a(bndm),a(bndspn),a(bndpar),a(enbnd),a(ncst),a(cntm)
     1 ,a(cntspn),a(cntpar),a(l0cnt),ldel,a(chnloc),a(jind),a(indxc),a
     2 (lch),z(ech),z(eng),z(en),z(rk),nsts,ncmax,ntchn,lplsmx,neng)
      call setwpt (z(pt),z(wt),z(wts),a(ncst),nsts,nptmx,msize)
*
*     ----- read in the local potential -----
      call etimes (tim)
      call rdlocp (z(locpot),vloc,ntchn,nptmx)
      call etimes (tim)
      write (iout,30) (tim(i),i=1,3)
*
*     ----- read in the bound state orbitals if needed -----
      call etimes (tim)
      if (orbctl.eq.'orb') then
      if (optctl.eq.'opt') then
      if (nfopt.ne.0) then
      call optinf (a(nopt),a(lstopt),a(kepvnl),a(kpdiag),nsts,nfopt)
      call orbin (z(orbtls),z(pt),a(nopt),a(lstopt),nsts,uperl,nptmx
     1 ,nfopt,ops)
      call etimes (tim)
      write (iout,40) (tim(i),i=1,3)
      endif
      endif
      if (nflag.ne.0) then
      call etimes (tim)
      call laginf (a(nlag),a(lstlag),nsts,nflag)
      call orbin (z(orblag),z(pt),a(nlag),a(lstlag),nsts,uperl,nptmx
     1 ,nflag,ops)
      call etimes (tim)
      write (iout,50) (tim(i),i=1,3)
      endif
      call iosys ('rewind all on mofile read-and-write',0,0,0,0)
      call iosys ('close mofile',0,0,0,0)
      endif
*
*     ----- all basic information is now in place -----
*     ----- begin energy dependent loop -----
*
*     ----- open file to hold wavefunction information -----
*
      call iosys ('open wavefn as new',262144,0,0,cntpsi)
      enadd=eng
      do 10 ieng=1,neng
*
*     ----- generate channel energies -----
*
c     the address for the energy changes with the loop index
c     use enadd as the integer address and add 2
c                   ^
c                   ^
c                   ^
      call geneng (z(enadd),z(enbnd),z(en),z(rk),a(ncst),nsts,ncmax)
*
*     ----- calculate greens function -----
*
      call etimes (tim)
      call grnset (z(gr1),z(gr2),z(en),z(rk),charge,z(pt),rmtrad,a(ncst)
     1 ,a(jind),a(chnloc),nsts,ntchn,nptmx,lplsmx,ncmax,ops)
*
*     ----- calculate the right hand sides -----
*
      call rhscal (z(gr1),z(gr2),z(rhs),a(ncst),a(jind),a(chnloc),nsts
     1 ,nptmx,msize,ntchn,ncmax,lplsmx,ops)
*
*     ----- calculate the inhomogeneities -----
*
      if (nfopt.ne.0) then
      call inhomo (z(gr1),z(gr2),z(orbtls),z(inhopt),z(pt),z(wt),a(ncst)
     1 ,a(nopt),a(chnloc),a(jind),nsts,ntchn,nptmx,uperl,nfopt,ncmax
     2 ,lplsmx,ops)
      endif
      if (nflag.ne.0) then
      call inhomo (z(gr1),z(gr2),z(orblag),z(inhlag),z(pt),z(wt),a(ncst)
     1 ,a(nlag),a(chnloc),a(jind),nsts,ntchn,nptmx,uperl,nflag,ncmax
     2 ,lplsmx,ops)
      endif
*
*     ----- read in the optical potential if needed -----
*
      if (optctl.eq.'opt') then
      call rdopt (z(optmat),optfle,z(scr3),a(nopt),a(kpdiag),a(kepvnl)
     1 ,nsts,nfopt,ieng,ops)
      endif
      call etimes (tim)
      write (iout,60) (tim(i),i=1,3)
*
*     ----- ok. all the easy work is done. -----
*     ----- lets get down to business and do the scattering calculation
*     ----- solve equations using iteration-variation method -----
*
      call etimes (tim)
      call slvit (z(gr1),z(gr2),z(orbtls),z(orblag),z(inhopt),z(inhlag)
     1 ,z(locpot),z(rhs),z(exvec),z(itvec),z(optmat),z(wts),z(smleq1),z
     2 (smleq2),z(smrhs1),z(smrhs2),a(iscrch),z(scr1),z(scr2),a(nbk),a
     3 (nbkpls),z(brk),z(c),z(sc),z(spln),a(nopt),a(nlag),a(jind),nsts,a
     4 (ncst),ncmax,msize,ntchn,nptmx,iter,maxvec,nfopt,nflag,nsol,
     5 numdo,uperl,nds1,nds2,nds3,nds4,convg,ovtol,ieng,ops,soldir,
     6 nowgt)
      call etimes (tim)
      write (iout,70) (tim(i),i=1,3)
*
*     ----- calculate overlaps of solutions and laguerre orbitals -----
*
      call etimes (tim)
      if (nflag.ne.0) then
      call overlp (z(ovlpa),a(scr1),z(rhs),z(orblag),z(wts),a(ncst),a
     1 (nlag),a(jind),nsts,a(nbk),a(nbkpls),z(brk),z(c),z(sc),z(spln)
     2 ,nflag,ncmax,ntchn,uperl,msize,nptmx,nsol,nds1,nds2,nds3,nds4
     3 ,ieng,ops)
      endif
*
*     ----- calculate the r- matrix and we are done -----
*
      call rmatrx (z(rmat),z(rhs),a(ncst),nsts,nsol,ntchn,nptmx,msize
     1 ,ieng,ops)
      call etimes (tim)
      write (iout,80) (tim(i),i=1,3)
*
*     -----               finished           -----
   10 enadd=enadd+1
      call iosys ('rewind all on wavefn read-and-write',0,0,0,0)
      call iosys ('close wavefn',0,0,0,0)
      return
*
*
*
*
   20 write (iout,100) nwords, canget
      stop 'memory'
*
c
   30 format (/,5x,'time to read local potential:',/,10x,'cpu:',f8.3,1x,
     1 'sys:',f8.3,1x,'aio:',f8.3)
   40 format (/,5x,'time to read optical orbs:',/,10x,'cpu:',f8.3,1x,'sy
     1s:',f8.3,1x,'aio:',f8.3)
   50 format (/,5x,'time to read in laguerre orbs:',/,10x,'cpu:',f8.3,1x
     1 'sys:',f8.3,1x,'aio:',f8.3)
   60 format (/,5x,'time to calculate rhs and inhomo:',/,10x,'cpu:',f8.3
     1 ,1x,'sys:',f8.3,1x,'aio:',f8.3)
   70 format (/,5x,'time to solve eqns:',/,10x,'cpu:',f8.3,1x,'sys:',f8.
     1 3,1x,'aio:',f8.3)
   80 format (/,5x,'time to extract final r-matrix:  ',/,10x,'cpu:',f8.3
     1 ,1x,'sys:',f8.3,1x,'aio:',f8.3)
   90 format (//,5x,'memory expanded by',1x,i8,1x,'words')
  100 format (//,5x,'need',1x,i8,1x,'words',2x,i8,1x,'words available',
     1        //,20x,'will quit')
      end
