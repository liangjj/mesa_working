c       10        20        30        40        50        60        70
c23456789012345678901234567890123456789012345678901234567890123456789012
       PROGRAM cknew
       IMPLICIT REAL*8 (a-h,o-z)
c
c*****************modified april1995********************************************
c     includes pseudopotential option
c
c*****************modified september1994*****************************************
c modified to handle closed channels explicitly
c for closed channels, ONLY L-2 functions appear and free functions are dropped
c
c*****************modified august 1992*****************************************
c incorporates either symmetry in denom or standard no symmetry runs in mesa
c also has an off-shell option and a static option
c
c   set istat = exchange   for an exact exchange run
c             = static     for a static run
c             = offshell   for an off-shell run
c             = other      for a standard run
c
c   set ihow  = mesa       for a mesa setup, WITH no solution of linear equations
c             = other      for a full run
c
c   set isym  = symmetry   IF mesa was run in symmetry
c             = other      IF mesa was not run in symmetry
c   
c*************modified to handle explicit inclusion of exchange*****************
c COMPLEX kohn scattering code started  8-8-88
c
c  reads bound-bound, bound-free, and free-free
c  matrices from the SEQUENCE of codes:
c newgrid, ass1, bespline, basis, statnu, ffbf, bndnew2, optpot
c
c changed 5/1/90 to READ nmotot BHL
c
c changed to READ in channel-dependent bound-bound matrices 12/21/88
c
c changed 12/22/88 to handle multi-channel bound-bound matrices
c correctly
c
c  addition of option to orthogonalize free functions to bound fcns
c  8-10-88
c
c  addition of eigenphase computation 10-3-88
c  note this uses routine cgeev for COMPLEX diagonalization from
c  SLATEC library
c
c optical potential added 1/19/89
c
c debugging for 3 channel CASE:  fixed fforth1 indexing 3/18/89
c                                fixed fforth2 indexing 3/18/89
c
c added total cross section  for multichannel CASE 3/18/89
c
c  ### MESA INTERFACE implemented 7/12/89 ####
c  reads file bndmat created my MESA code m950
c
       PARAMETER (nblok=1000,mtop=#maxltop,ltop=#maxltop,nchnl=#maxchan)
c      PARAMETER (lmtop=(2*mtop+1)*(ltop+1)-mtop*(mtop+1))
       PARAMETER (lmtop=#maxlmtop)
       PARAMETER (nstate=nchnl*(nchnl+1)/2)
       PARAMETER (nbfmax=#maxnbfkohn)
c      PARAMETER (nbig=nchnl*(lmtop+nbfmax))
c      PARAMETER (nbig=700)
c      PARAMETER (nsmall=nchnl*lmtop)
      PARAMETER (nbig=#maxbig)
      PARAMETER (nsmall=#maxsmall)
      INTEGER fmt(10),iclosed(nchnl)
      INTEGER nlm(nchnl),ngauss(nchnl),ngch(nbfmax,nchnl)
      INTEGER nscat(nchnl),nsch(nbfmax,nchnl)
      INTEGER lch(lmtop,nchnl),mch(lmtop,nchnl)
      INTEGER ntarg(nchnl),ntch(nbfmax,nchnl)
      INTEGER nbas(nchnl),nbch(nbfmax,nchnl),kpvt(nbig)
      CHARACTER*1 jobvl, jobvr
      REAL*8 kchan(nchnl),kchanp(nchnl)
      REAL*8 scratch(nbfmax,nbfmax,nstate)
      COMPLEX*16 unity(nsmall,nsmall)
      REAL*8 tr(nbfmax,nbfmax)
      REAL*8 hbb(nbfmax**2),hbbe(nbfmax**2)
      REAL*8 hbbtrnd(nbfmax,nbfmax,nstate),hbbtrne(nbfmax,nbfmax,nstate)
      REAL*8 obbtrn(nbfmax,nbfmax,nstate),hbbtrn(nbfmax,nbfmax,nstate)
      REAL*8 vopt(nbig,nbig)
      REAL*8 xsecmat(nchnl,nchnl)
      REAL*8 hpptemp(nbig,nbig)
      COMPLEX*16 ovbf(lmtop,nbfmax,nchnl),ovbfp(lmtop,nbfmax,nchnl)
      COMPLEX*16 hpvb(lmtop,nbfmax,nchnl**2)
      COMPLEX*16 hpvbp(lmtop,nbfmax,nchnl**2)
      COMPLEX*16 hpvbtrn(lmtop,nbfmax,nchnl**2)
      COMPLEX*16 hpvbtrnp(lmtop,nbfmax,nchnl**2)
      COMPLEX*16 ovbftrn(lmtop,nbfmax,nchnl)
      COMPLEX*16 ovbftrnp(lmtop,nbfmax,nchnl)
      COMPLEX*16 hpvhp(lmtop,lmtop,nstate)
      COMPLEX*16 hpvhm(lmtop,lmtop,nchnl**2)
      COMPLEX*16 hpvhmp(lmtop,lmtop,nchnl**2)
      COMPLEX*16 hdenom(nbig,nbig),cov(nbig,nsmall)
      COMPLEX*16 htop(nbig,nsmall)
      COMPLEX*16 bforth(lmtop,nbfmax,nchnl**2)
      COMPLEX*16 bforthp(lmtop,nbfmax,nchnl**2)
      COMPLEX*16 hpvhmps(lmtop,lmtop,nchnl)
      COMPLEX*16 zdotc,zdotu,ai, det, htop2(nbig,nsmall)
      COMPLEX*16 smat(nsmall,nsmall),work(nbig)
      COMPLEX*16 seig(nsmall),svec(nsmall,nsmall)
      COMPLEX*16 swork(2*nsmall)
      REAL*8 srwork(2*nsmall)
c      REAL*8 swork(3*nsmall)
       EQUIVALENCE (vopt,hbb)
       EQUIVALENCE (unity, hdenom),(hbbtrne,hbbe)
       EQUIVALENCE (hpptemp,hdenom)
       EQUIVALENCE (ovbf, hpvb, hpvhm)
       EQUIVALENCE (ovbfp, hpvbp, hpvhmp)
       EQUIVALENCE (hpvhp,scratch),(hpvbtrnp,cov)
c      EQUIVALENCE (fpfmorth,hdenom),(fpfporth,smat),(bforth,htop)
       EQUIVALENCE (bforth,htop),(smat,ovbftrn)
     1 ,(bforthp,htop2)
       CHARACTER*8 ioft,iex,itell,jsym,istatic,ihow,istat,isym
c **       DATA ioft/8hoffshell/
c **       DATA iex/8hexchange/,itell/8hmesa    /,jsym/8hsymmetry/
c **       DATA istatic/8hstatic  /

       ioft = 'offshell'
       iex = 'exchange'
       itell = 'mesa    '
       jsym = 'symmetry'
       istatic = 'static  '

       ai=CMPLX(0.,1.)
       OPEN(5,file='inkohn',status='unknown')
       OPEN(6,file='outkohn',status='unknown')
       OPEN(63,file='eigkohn',status='unknown')
       READ(5,*)iprint,iorth,iopt,ipseudo
       READ(5,734) istat
       IF(istat.EQ.ioft)THEN
       OPEN(77,file='intsffof',form='unformatted',status='old') 
       OPEN(66,file='intsbfof',form='unformatted',status='old') 
       ENDIF
 734   FORMAT(a8)
       WRITE(6,734) istat
       IF(istat.EQ.iex.OR.istat.EQ.ioft)THEN
        OPEN(15,file='exchmat',form='unformatted',status='old')
       ENDIF
       IF(iopt.EQ.1)THEN
       WRITE(6,595)
  595  FORMAT(/" this is an optical potential calculation"/)
       ENDIF
       IF(ipseudo.NE.0)THEN
       WRITE(6,596)
 596   FORMAT(/" this is an effective core potential calculation"/)
       ELSE
       WRITE(6,597)
  597 FORMAT(/" this is not an effective core potential calculation"/)
       ENDIF
       READ(5,734)ihow
       WRITE (6, 598) ihow
 598   FORMAT("ihow is ", a)
       IF(itell.EQ.ihow)THEN
          WRITE (6, 599) "key word is mesa"
 599      FORMAT(a)
          OPEN(7,file='bmesa',form='unformatted',status='unknown')
          OPEN(78,file='cmesa',form='unformatted',status='unknown')
          IF(istat.EQ.istatic)THEN
             OPEN(99,file='bstatic',form='unformatted',status='unknown')
          ENDIF
       ELSE
          WRITE (6, 599) "key word is not mesa"
          OPEN(88,file='tmats',status='unknown')
          OPEN(7,file='punkohn',form='unformatted',status='unknown')
       ENDIF
       READ(5,734)isym
       WRITE(6,734)isym
       OPEN(8,file='intsbf',form='unformatted',status='old') 
       OPEN(9,file='bndmat',form='unformatted',status='old')
       OPEN(10,file='intsff',form='unformatted',status='old') 
       IF(ipseudo.NE.0)THEN
          OPEN(11,file='intsbfps',form='unformatted',status='old') 
          OPEN(12,file='intsffps',form='unformatted',status='old') 
       ENDIF
c
       READ(10) ibcondx
       IF(ibcondx.EQ.0) WRITE(6,731)
 731   FORMAT(//,' S-matrix boundary conditions')
       IF(ibcondx.EQ.1) WRITE(6,732)
 732   FORMAT(//,' T-matrix boundary conditions')
       IF(ibcondx.NE.0.AND.ibcondx.NE.1) THEN
       WRITE(6,733)
 733   FORMAT(//,' Error  boundary cond. on file intsff is nonsense')
       STOP
       ENDIF
       READ(8) nener,nchan,(nlm(ic),ic=1,nchan)
c
       IF(nchan.GT.nchnl) THEN
       WRITE(6,*)' Error nchan gt nchnl ',nchan,nchnl
       STOP 11
       END IF
       READ(8) ((lch(j,ic),mch(j,ic),j=1,nlm(ic)),ic=1,nchan)
       READ(8) eground
       READ(8) (ngauss(ic),ic=1,nchan)
       READ(8) ((ngch(j,ic),j=1,ngauss(ic)),ic=1,nchan)
       nbmax=0
       nlmax=0
       DO 851 i=1,nchan
       IF(nlm(i).GT.nlmax)nlmax=nlm(i)
       IF(ngauss(i).GT.nbmax)nbmax=ngauss(i)
  851     CONTINUE
       IF(nlmax.GT.lmtop)THEN
       WRITE(6,*)' Error nlmax gt lmtop', nlmax,lmtop
       STOP
       ENDIF
       IF(nbmax.GT.nbfmax)THEN
       WRITE(6,*)' Error nbmax gt nbfmax', nbmax,nbfmax
       STOP
       ENDIF
       READ(10)
       READ(10)
       READ(10)
       IF(ipseudo.NE.0)THEN
          READ(11)
          READ(11)
          READ(11)
          READ(11)
          READ(11)
          READ(12)
          READ(12)
          READ(12)
          READ(12)
       ENDIF
c
c  READ ASSIGNMENT of scattering basis to channels
c
       READ(5,*) (nscat(ic),ic=1,nchan)
       WRITE(6,138) (nscat(ic),ic=1,nchan)
138   FORMAT(' number of scattering orbitals in each channel:', 6i3)
       READ(5,*) ((nsch(j,ic),j=1,nscat(ic)),ic=1,nchan)
       DO 875 ic=1,nchan
       WRITE(6,139) ic, (nsch(j,ic),j=1,nscat(ic))
 139   FORMAT(' channel',i3,' orbs:',20i3,/,(17x,20i3))
  875   CONTINUE
c
c READ ASSIGNMENT of TARGET orbitals to channels
c
       IF(iorth.NE.0) THEN
       READ(5,*) (ntarg(ic),ic=1,nchan)
       itarg=0
       WRITE(6,137) (ntarg(ic),ic=1,nchan)
 137   FORMAT(' number of target orbitals in each channel:', 6i3)
       DO 5 i=1,nchan
 5     itarg = itarg + ntarg(i)
       IF(itarg.NE.0) THEN
       READ(5,*) ((ntch(j,ic),j=1,ntarg(ic)),ic=1,nchan)
       DO 876 ic=1,nchan
       WRITE(6,139) ic, (ntch(j,ic),j=1,ntarg(ic))
 876   CONTINUE
       ENDIF
       ELSE
       DO 237 ic=1,nchan
 237   ntarg(ic) = 0
       ENDIF
c
c construct total orbital list for orthogonalization of free fcns.
c IF iorth was zero this list is the same as nscat list
c
       nntot=0
       nbtot=0
       mfree=0
       DO 4 ic=1,nchan
       nntot=nntot+nlm(ic)+nscat(ic)    
       nbtot=nbtot+nscat(ic)
       mfree=mfree+nlm(ic)
       ntic=ntarg(ic)
       nsic=nscat(ic)
       nbas(ic)=ntic+nsic
       DO 2 i=1,nsic
 2     nbch(i,ic)=nsch(i,ic)
       IF (ntic.EQ.0.OR.itarg.EQ.0) go to 4
       DO 3 i=1,ntic
 3     nbch(nsic+i,ic) = ntch(i,ic)
 4     CONTINUE
       IF(nntot.GT.nbig)THEN
       WRITE(6,*)' Error stopping because ntot gt nbig', nntot,nbig
       STOP
       ELSE
       WRITE(6,*)' total dimension of hamiltonian is ',nntot
       ENDIF
c
       IF(ihow.EQ.itell)THEN
       WRITE(7)nchan,nbtot,mfree,(nlm(i),i=1,nchan)
       ENDIF
       IF(iorth.NE.0) THEN
       WRITE(6,787)
 787   FORMAT(//,' *****orthogonalization is on*****',//,
     1 ' free functions will be orthogonalized to following orbitals')
       DO 788 ic=1,nchan
 788   WRITE(6,789) ic,(nbch(i,ic),i=1,nbas(ic))
 789   FORMAT(' in channel',i3,' :  ',10i3/(18x,10i3))
       WRITE(6,791)
 791   FORMAT(//)
       ENDIF
c 
c   skip some preliminary info on file bndmat made by m950
c
       READ(9)
       READ(9)
       READ(9)
       READ(9)
       IF(istat.EQ.iex.OR.istat.EQ.ioft)THEN
       READ(15)
       ENDIF
c
c READ in transformation matrix
c
c..bhl 5/1/90
      READ(9) nbftot,nmotot
c
      IF(nbftot.GT.nbfmax) THEN
      WRITE(6,*)' Error nbftot gt nbfmax ',nbftot,nbfmax
      STOP 11
      END IF
c..bhl      nmotot=nbftot
      nsy = nbftot
      CALL rdbinsqr(tr,nbfmax,nbftot,9)
      WRITE(6,99)
99    FORMAT(" transformation matrix")
      DO 1 i=1,nmotot
      WRITE(6,100)i,(tr(j,i),j=1,nbftot)
1     CONTINUE
100   FORMAT(i5,6f12.7/(5x,6f12.7))
c
c READ in bound-bound overlap matrix
c
      CALL rdbintri(hbb,nsy,9)
c
c  transform overlap matrix as a test
c
      CALL ovtrans(hbb,obbtrn,nsy,nchan,nbas,nbch,tr,scratch,nbfmax
     1 ,nstate,iprint)
c
c READ in energy-INDEPENDENT bound-bound direct hamiltonian
c  and transform hbb direct
c
      CALL bbtrans(hbb,hbbtrnd,nsy,nchan,nbas,nbch,tr,scratch,nbfmax
     1 ,nstate,iprint,hbbe,hbbtrne,istat)
c
c OPEN a loop on energies
c
      IF(itell.EQ.ihow)THEN
      WRITE(7)nener
      ENDIF
      DO 1000 iene=1,nener
      READ(8) (kchan(ic),ic=1,nchan)
      READ(8) (iclosed(ic),ic=1,nchan)
      IF(ihow.EQ.itell)WRITE(78) (iclosed(ic),ic=1,nchan)
      READ(10)
      IF(ipseudo.NE.0)THEN
         READ(11)
         READ(12)
      ENDIF
      WRITE(6,111) (kchan(i),i=1,nchan)
c
c compute the number OPEN OPEN channels at this energy
c
      nopen=nchan
      DO 1001 iopen=1,nchan
 1001    nopen=nopen-iclosed(iopen)
      WRITE(6,*)'nopen=',nopen   
      IF(itell.EQ.ihow)THEN
      WRITE(7)(kchan(ic),ic=1,nchan)
      ENDIF
111   FORMAT(///,' ********************************',
     # ' new energy ********************************',
     # //,' channel momenta:',6f12.5)
      IF(istat.EQ.ioft)THEN
      READ(66)(kchanp(i),i=1,nchan)
      READ(77)(kchanp(i),i=1,nchan)
      WRITE(6,112)(kchanp(i),i=1,nchan)
 112  FORMAT(' off-shell momenta:',6f12.4)
      ENDIF
      DO 507 ic=1,nchan
      nlmic=nlm(ic)
      ngic=ngauss(ic)
      READ(8)((ovbf(ilm,ig,ic),ilm=1,nlmic),ig=1,ngic)
      IF(istat.EQ.ioft)THEN
      READ(66)((ovbfp(ilm,ig,ic),ilm=1,nlmic),ig=1,ngic)
      ENDIF
507   CONTINUE
c
c transform the bound-free overlap matrices
c
      CALL obftrans(ovbf,ovbftrn,tr,nsy,nchan,nbas,nbch,nlm,ngauss,
     1 ngch,nbfmax,lmtop,nchnl,iprint,ovbftrnp,ovbfp,istat,ioft)
      DO 504 ic=1,nchan
      DO 504 jc=1,nchan
      ist=nchan*(ic-1)+jc
      nlmic=nlm(ic)
      ngjc=ngauss(jc)
      READ(8)((hpvb(ilm,ig,ist),ilm=1,nlmic),ig=1,ngjc)
      IF(istat.EQ.ioft)THEN
      READ(66)((hpvbp(ilm,ig,ist),ilm=1,nlmic),ig=1,ngjc)
      ENDIF
504   CONTINUE
c     
c add pseudopotential terms to bound-free integrals
c     
      IF(ipseudo.NE.0)THEN
         DO 505 ic=1,nchan
            nlmic=nlm(ic)
            ngic=ngauss(ic)
            READ(11)((hpvbtrn(ilm,ig,ic),ilm=1,nlmic),ig=1,ngic)
            ist=nchan*(ic-1)+ic
            DO 505 ilm=1,nlmic
               DO 505 ig=1,ngic
 505           hpvb(ilm,ig,ist)=hpvb(ilm,ig,ist)+hpvbtrn(ilm,ig,ic)
      ENDIF
c
c transform the bound-free (h-e) matrices
c
      CALL hbftrans(hpvb,hpvbtrn,tr,nsy,nchan,nbas,nbch,nlm,ngauss,
     1 ngch,nbfmax,lmtop,nchnl,iprint,hpvbp,hpvbtrnp,istat,ioft)
c
c READ the free free integrals file
c
      DO 604 ic=1,nchan
      DO 604 jc=1,ic
      ist = ic*(ic-1)/2 +jc
      nlmic=nlm(ic)
      nlmjc=nlm(jc)
      READ(10) ((hpvhp(ilm,jlm,ist),ilm=1,nlmic),jlm=1,nlmjc)
604   CONTINUE
c     
c add pseudopotential terms to free-free integrals
c
      IF(ipseudo.NE.0)THEN
      DO 665 ic=1,nchan
         nlmic=nlm(ic)
         READ(12) ((hpvhm(ilm,jlm,ic),ilm=1,nlmic),jlm=1,nlmic)
         ist=ic*(ic+1)/2
         DO 665 ilm=1,nlmic
            DO 665 jlm=1,nlmic
 665        hpvhp(ilm,jlm,ist)=hpvhp(ilm,jlm,ist)+
     $          hpvhm(ilm,jlm,ic)
      ENDIF
      DO 608 ic=1,nchan
      DO 608 jc=1,nchan
      ist=nchan*(ic-1)+jc
      nlmic=nlm(ic)
      nlmjc=nlm(jc)
      READ(10) ((hpvhm(ilm,jlm,ist),ilm=1,nlmic),jlm=1,nlmjc)
      IF(istat.EQ.ioft)THEN
      READ(77) ((hpvhmp(ilm,jlm,ist),ilm=1,nlmic),jlm=1,nlmjc)
      ENDIF
 608  CONTINUE
c     
c add pseudopotential terms to free-free integrals
c
      IF(ipseudo.NE.0)THEN
      DO 668 ic=1,nchan
         nlmic=nlm(ic)
         READ(12) ((hpvhmps(ilm,jlm,ic),ilm=1,nlmic),jlm=1,nlmic)
         ist=nchan*(ic-1)+ic
         DO 668 ilm=1,nlmic
            DO 668 jlm=1,nlmic
 668        hpvhm(ilm,jlm,ist)=hpvhm(ilm,jlm,ist)+
     $          hpvhmps(ilm,jlm,ic)
      ENDIF
c
c ***orthogonalization section***
c ***must USE for static exchange***
c
c construct bound-free (h-e) matrices for free fcns orthogonalized
c bound functions
c
      IF(iorth.NE.0) THEN
      CALL bforthog(ovbftrn,lmtop,nbfmax,nchnl,nstate,hpvbtrn,
     1 hbbtrnd,nchan,nlm,nbas,nbas,bforth,iprint,bforthp,
     2 obbtrn,kchan,eground,ovbftrnp,hbbtrne,hpvbtrnp,istat,ioft,
     $ nscat)
c
c construct free-free (h-e) matrices for free functions orthogonalized
c to bound functions
c
      CALL fforth1(ovbftrn,lmtop,nbfmax,nchnl,nstate,
     1 hpvbtrn,bforth,nchan,nlm,nbas,hpvhp,iprint)
      CALL fforth2(ovbftrn,lmtop,nbfmax,nchnl,nstate,
     1 hpvbtrn,bforth,nchan,nlm,nbas,hpvhm,iprint,
     2 ibcondx,hpvbtrnp,hpvhmp,istat,ioft,nscat,ovbftrnp)
c
c load matrices WITH orthogonalized free functions back into
c the arrays which held the same matrices WITH unorthogonalized free functions
c so that remainder of the code is INDEPENDENT of orthogonalization PROCEDURE
c
      IF(istat.NE.ioft)THEN
      DO 452 ic=1,nchan
      nlmic=nlm(ic)
      DO 452 jc=1,nchan
      nlmjc=nlm(jc)
      icc=nchan*(ic-1)+jc
      DO 452 ilm=1,nlmic
      DO 452 jlm=1,nlmjc
       hpvhmp(ilm,jlm,icc)=hpvhm(ilm,jlm,icc)
452   CONTINUE
      ENDIF
      DO 453 ic=1,nchan
      nlmic=nlm(ic)
      DO 453 jc=1,nchan
      nsjc=nscat(jc)
      icc=nchan*(ic-1)+jc
      DO 453 ilm=1,nlmic
      DO 453 jsc=1,nsjc
      IF(istat.NE.ioft)THEN
      hpvbtrnp(ilm,jsc,icc)=bforth(ilm,jsc,icc)
       ELSE
       hpvbtrnp(ilm,jsc,icc)=bforthp(ilm,jsc,icc)
      ENDIF
453   hpvbtrn(ilm,jsc,icc)=bforth(ilm,jsc,icc)
c     
      ELSE
         IF(istat.NE.ioft)THEN
      DO 943 ic=1,nchan
      nlmic=nlm(ic)
      DO 943 jc=1,nchan
      nlmjc=nlm(jc)
      icc=nchan*(ic-1)+jc
      DO 943 ilm=1,nlmic
      DO 943 jlm=1,nlmjc
 943     hpvhmp(ilm,jlm,icc)=hpvhm(ilm,jlm,icc)
      DO 944 ic=1,nchan
      nlmic=nlm(ic)
      DO 944 jc=1,nchan
      nsjc=nscat(jc)
      icc=nchan*(ic-1)+jc
      DO 944 ilm=1,nlmic
      DO 944 jsc=1,nsjc
 944  hpvbtrnp(ilm,jsc,icc)=hpvbtrn(ilm,jsc,icc)
      ENDIF
      ENDIF
c
c  *** END of orthogonalization section ***
c
c construct denominator
c
      CALL denom(hdenom,nbig,hpvhp,lmtop,nstate,hpvbtrn,nbfmax,nchnl
     1 ,hbbtrn,nchan,nscat,nlm,kchan,eground,obbtrn,iprint,ntot,vopt,
     1 iopt,nsch,hpptemp,nmotot,hbbtrnd,istat,itell,ihow,isym,jsym,
     $ hbbtrne,nopen)
c
c cnstruct numerator
c
      CALL top(htop,nbig,nsmall,hpvhm,lmtop,nstate,hpvbtrn,
     1 nbfmax,nchnl,nchan,nscat,nlm,iprint,ntot,nfree,ibcondx
     2 ,htop2,hpvhmp,hpvbtrnp,itell,ihow,ovbftrnp,istat,ioft,cov,
     $ nbas,nopen)
      IF(mfree.NE.nfree)THEN
      WRITE(6,*)' note that nfree and mfree are different', mfree,nfree
      ENDIF
c
c IF this is a mesa run, WRITE out born terms and quit
c
      IF(itell.EQ.ihow)THEN
         iprev = 0
         DO 651 ic=1,nopen
            nlmic=nlm(ic)
            jprev = 0
            DO 650 jc=1,nopen
               nlmjc = nlm(jc)
               ii=max0(ic,jc)
               jj=ic+jc-ii
               ist=ii*(ii-1)/2 + jj
               icc=nchan*(ic-1) + jc
               DO 640 ilm=1,nlmic
                  isub=ilm+iprev
                  IF(ibcondx.EQ.0) THEN
                     DO 638 jlm=1,nlmjc
                        jsub = jlm+jprev
c******note**************
c this is incorrect for off-shell s-matrix because hpvhpp was never calculated
c
  638                 smat(isub,jsub) =CONJG(hpvhp(ilm,jlm,ist))
                   ELSE
                      DO 639 jlm=1,nlmjc
                         jsub = jlm+jprev
 639                  smat(isub,jsub) =imag(hpvhmp(ilm,jlm,icc))
                   ENDIF
 640            CONTINUE
 650         jprev=jprev+nlmjc
 651      iprev=iprev+nlmic
          WRITE(7)((smat(i,j),i=1,nfree),j=1,nfree)
          WRITE(6,*)' Born matrix'
          DO 59 i=1,nfree
 59          WRITE(6,108) i,(-2.*smat(i,j),j=1,nfree)

      ENDIF
      IF(itell.EQ.ihow)go to 1000
c
c SAVE numerator matrix
c
      nterm=nbig*nsmall
c      CALL ccopy(nterm,htop,1,htop2,1)
c
c solve linear equations
c
c     CALL matinv(hdenom,ntot,htop2,nfree,det,nbig)

c      WRITE (6, 669) ntot, nfree
c 669  FORMAT("hdenom  ntot=", i10, "  nfree =", i10)
c      DO 667 i = 1, ntot
c         WRITE (6, 671) i, (hdenom(j, i), j = 1, ntot)
c 667  CONTINUE
c 671  FORMAT("Row", i5, /, (4e17.8))
c      WRITE (6, 673)
c 673  FORMAT("htop")
c      DO 672, i = 1, nfree
c         WRITE (6, 674) i, (htop2(j, i), j = 1, ntot)
c 674     FORMAT("Row", i5, /, (4e17.8))
c 672  CONTINUE

      CALL cgefs(hdenom,nbig,ntot,htop2(1,1),1,ind,work,kpvt)
      DO 666 i=2,nfree
 666     CALL cgefs(hdenom,nbig,ntot,htop2(1,i),2,ind,work,kpvt)
c
      IF(iprint.NE.0) THEN
      WRITE(6,107)
107   FORMAT(//' solution matrix ')
      DO 60 i=1,ntot
60    WRITE(6,108) i,(htop2(i,j),j=1,nfree)
108   FORMAT(1x,i3,6f12.5,/,(4x,6f12.5))
      ENDIF
      DO 500 ilm=1,nfree
      DO 500 jlm=1,nfree
      smat(ilm,jlm)=zdotu(ntot,htop2(1,jlm),1,htop(1,ilm),1)
500   CONTINUE
      IF(iprint.NE.0) THEN
      WRITE(6,109)
109   FORMAT(//' scattered wave part  ')
      DO 61 i=1,nfree
61    WRITE(6,108) i,(smat(i,j),j=1,nfree)
      ENDIF
c
c add born terms
c
      iprev = 0
      DO 551 ic=1,nopen
      nlmic=nlm(ic)
      jprev = 0
      DO 550 jc=1,nopen
      nlmjc = nlm(jc)
      ii=max0(ic,jc)
      jj=ic+jc-ii
      ist=ii*(ii-1)/2 + jj
      icc=nchan*(ic-1) + jc
      DO 540 ilm=1,nlmic
      isub=ilm+iprev
      IF(ibcondx.EQ.0) THEN
      DO 538 jlm=1,nlmjc
      jsub = jlm+jprev
c******note**************
c this is incorrect for off-shell s-matrix because hpvhpp was never calculated
c
538   smat(isub,jsub) =ai*(-smat(isub,jsub)+CONJG(hpvhp(ilm,jlm,ist)))
      ELSE
      DO 539 jlm=1,nlmjc
      jsub = jlm+jprev
 539   smat(isub,jsub) =-2.0*(-smat(isub,jsub)+imag(hpvhmp(ilm,jlm
     1 ,icc)))
      ENDIF
540   CONTINUE
550   jprev=jprev+nlmjc
  551 iprev=iprev+nlmic
c
c WRITE out s matrix
c
      IF(ibcondx.EQ.0) THEN
      WRITE(6,776) (kchan(i),i=1,nopen)
776   FORMAT(//,' S-matrix for channel momenta:',3f12.6)
      ELSE
c
c  Here we have the T-matrix.  Compute the cross sections from it.
c
      istart = 0
      DO 879 ic=1,nopen
      jstart = 0
      ni = nlm(ic)
      DO 878 jc=1,nopen
      nj=nlm(jc)
      ist=istart+1
      ifin=istart+ni
      jst=jstart+1
      jfin=jstart+nj
      WRITE(7)ic,jc,ni,nj,kchan(ic),kchan(jc)
      WRITE(88,177)ic,jc,ni,nj,kchan(ic),kchan(jc)
 177  FORMAT(4i5,2f20.10)
      WRITE(7)((smat(ii,jj),ii=ist,ifin),jj=jst,jfin)
      WRITE(88,277)((smat(ii,jj),ii=ist,ifin),jj=jst,jfin)
 277  FORMAT(4e20.10)
      summod = 0.
      DO 877 ilm=1,ni
      DO 877 jlm=1,nj
      isub = istart + ilm
      jsub = jstart + jlm
877   summod = summod + ABS(smat(isub,jsub))**2
      xsecmat(ic,jc) = 4.0*3.141592654*summod/kchan(ic)**2
      jstart = jstart + nlm(jc)
878   CONTINUE
      istart = istart + nlm(ic)
879   CONTINUE
       WRITE(6,874)
874   FORMAT(//,' total cross sections: row index = initial chnl,',
     # ' column index = final chnl')
      WRITE(6,871) (i,i=1,nopen)
871   FORMAT(4x,6(6x,i2,4x))
      crpt=xsecmat(1,1)
      eept=(kchan(1)**2)/2.
      DO 873 i=1,nopen
      WRITE(6,872) i, (xsecmat(i,j),j=1,nopen)
873   CONTINUE
872   FORMAT(1x,i2,1x,6f12.5/(4x,6f12.5))
c
      WRITE(6,778) (kchan(i),i=1,nopen)
778   FORMAT(//,' T-matrix for channel momenta:',3f12.6)
      ENDIF
      DO 570 i=1,nfree
570   WRITE(6,777) i,(smat(j,i),j=1,nfree)
777   FORMAT(1x,i3,6f12.5,/,(4x,6f12.5))
c
c  test unitarity
c
c watch out; unity(i,j) is equivalenced to denom(i,j)
c
      IF(ibcondx.NE.0) THEN
      DO 795 i=1,nfree
      DO 795 j=1,nfree
795   smat(i,j) =  (0.0,2.0)*smat(i,j)
      DO 796 i=1,nfree
  796 smat(i,i) = smat(i,i) + 1.0
      ENDIF
c
c  we now own the S-matrix -- no matter what TYPE of boundary
c  we used in the calculation
c  first test unitarity and THEN compute eigenphases
c
      DO 792 i=1,nfree
      DO 792 j=1,nfree
  792 unity(i,j)=zdotc(nfree,smat(1,i),1,smat(1,j),1)
      WRITE(6,793)
793   FORMAT(//,' unitarity check: Hermitian conjg(S) * S')
      DO 794 i=1,nfree
794   WRITE(6,777) i,(unity(j,i),j=1,nfree)
c      job = 0
c      CALL cgeev(smat,nsmall,nfree,seig,svec,nsmall,swork,job,inform)
	WRITE(6,*)' S-matrix'
      DO 670 i=1,nfree
670   WRITE(6,777) i,(smat(j,i),j=1,nfree)
      jobvl='n'
      jobvr='v'
      lwork=2*nsmall
      CALL zgeev(jobvl,jobvr,nfree,smat,nsmall,seig,svec,nsmall,
     &           svec,nsmall,swork,lwork,srwork,inform)
      IF(inform.NE.0) WRITE(6,1299) inform
1299  FORMAT(//,' *** S-matrix diagonalization failed, inform =',i5)
      phasesum = 0.0
      DO 573 i = 1,nfree
         arg1=imag(seig(i))
         arg2=DBLE(seig(i))
      phase = ATAN2(arg1,arg2)
      phase=phase*0.5
      rtest = ABS(seig(i))
      phasesum=phasesum + phase
      WRITE(6,574) i, phase,rtest
574   FORMAT(' eigenphase #',i3,' =',f12.5,' modulus =',f12.5)
573   CONTINUE
      WRITE(63,587) eept,crpt,phasesum
 587  FORMAT(3f15.7)
      WRITE(6,576)phasesum,(kchan(i),i=1,nopen)
  576 FORMAT(' eigenphase sum =',f12.5,
     # ' for channel momenta:', (6f12.5))
      IF (ibcondx.NE.0) THEN
      eincidnt = kchan(1)**2/2.
884   FORMAT(' incident energy =',f15.10)
      ENDIF
1000  CONTINUE
c
      STOP
      END
