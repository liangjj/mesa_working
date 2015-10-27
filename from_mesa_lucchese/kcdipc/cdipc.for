         PROGRAM dipc
         IMPLICIT REAL*8 (a-h,o-z)
         CHARACTER*8 iword
c
c computes bound-free dipole matrix elements
c
c created by modifying ffbfc 5-9-89 aeo
c
c tnr** added closed-channel modifications October, 2008
c
c  reads file of r grid (transformed)
c  reads file of ylms
c  reads file of spline coefficients for bessel functions computation
c  reads file of L**2 basis functions on grid
c
c  changed to READ in lch, mch channel quantum numbers 8-2-88
c
c COMPLEX io changed to REAL io because of system error 11-30-88
c
c
c parameters:
c
      PARAMETER (nblok=500,mtop=#maxltop,ltop=#maxltop,nchnl=#maxchan)
      PARAMETER (nsblok=25000)
c      PARAMETER (lmtop=(2*mtop+1)*(ltop+1)-mtop*(mtop+1))
      PARAMETER (lmtop=#maxlmtop)
      PARAMETER (lm2=2*lmtop)
      PARAMETER (nbfmax=#maxnbfkohn)
      PARAMETER (maxene=200)
c
c  nblok = number of points in BLOCK of grid etc. values
c  mtop = largest m VALUE
c  ltop = largest l VALUE
c  nchnl = largest number of channels allowed
c  nbfmax  =  largest number of l**2 basis functions allowed
c      (transformation to MO's occurs in scattering code)
c  nsblok = max number of spline points for bessel fcns
c  maxene = maximum # of incident energies
c
      complex*16 hp(nblok,0:ltop,nchnl)
      complex*16 cj(nsblok,0:ltop,nchnl)
      complex*16 hs(nsblok,0:ltop,nchnl)
      real*8 cjr(2*nsblok,0:ltop,nchnl), cy(nsblok,0:ltop,nchnl)
      real*8 hsr(2*nsblok,0:ltop,nchnl), hsder(nsblok,0:ltop,nchnl)
      real*8 rvec(nblok),krvec(nblok), x(nsblok)
      integer ipow(0:ltop),nscat(nchnl),ntarg(nchnl),nbas(nchnl)
      integer nsch(nbfmax,nchnl),ntch(nbfmax,nchnl),nbch(nbfmax,nchnl)
      real*8 grid(nblok,3),tr(nbfmax,nbfmax)
      real*8 wt(nblok),cvec(nblok),dvec(nblok)
      real*8 echan(nchnl), energy(maxene)
      real*8 kchan(nchnl)
      real*8 ylm(nblok,0:ltop,0:2*ltop)
      real*8 buff(nblok*4),biff(2*(nblok*(ltop+1)))
      real*8 basis(nbfmax*nblok)
      real*8 rit(3)
      integer it(3),iclosed(nchnl)
      dimension nlm(nchnl),lch(lmtop,nchnl),mch(lmtop,nchnl)
      complex*16 xvec(nblok),yvec(nblok),zvec(nblok)
        dimension hbx(nbfmax,nbfmax),hbtx(nbfmax,nbfmax)
        dimension hby(nbfmax,nbfmax)
        dimension hbz(nbfmax,nbfmax)
        complex*16 hpvb(lmtop,nbfmax,nchnl)
        complex*16 hpvbx(lmtop,nbfmax,nchnl)
        complex*16 hpvby(lmtop,nbfmax,nchnl)
        complex*16 hpvbz(lmtop,nbfmax,nchnl)
        complex*16 hpvbt(lmtop,nbfmax,nchnl)
        complex*16 hpvbtx(lmtop,nbfmax,nchnl)
        complex*16 hpvbty(lmtop,nbfmax,nchnl)
        complex*16 hpvbtz(lmtop,nbfmax,nchnl)
      common/parms/mgrid,npt,nbf
      common/parmbes/ngrid2,nchunk,nener,nchan,lmax2
      common/parm/nbig,mpt,lmax,mumax
      data pi/3.14159265358979/
      equivalence (cj(1,0,1),cjr(1,0,1)),(hs(1,0,1),hsr(1,0,1))
c
c all io units defined here (no calls to "create")
c
c      call link("unit6=(outdipc,create,hc),unit9=(grid,abs,open)
c     1 ,unit10=(ylms,abs,open),unit13=(vbas,abs,open)
c     3 ,unit8=dipbf,unit5=indipc,unit88=besspl//")
c
      open(unit=5,file='indipc',form='formatted')
c tnr**********
      read(5,*)npt
      irecl=4*npt*8 
c tnr**********
      open(unit=6,file='outdipc',form='formatted')
      call openabs(9,'grid',irecl)
c tnr**********
      call rdabs(9,rit,3,0)
      it(1)=rit(1)
      it(2)=rit(2)
      it(3)=rit(3)
      write(6,*)it
      ngrid=it(3)
      if(it(1).ne.npt)then
         write(6,*)" stopping because grid buffer length is wrong"
         write(6,*)it(1),it(2),it(3)
         write(6,*)npt,irecl
         stop
      endif
c tnr**********
c	call openabs(10,'ylms')
c	call openabs(13,'vbas')
      open(10,file='ylms',form='unformatted',status='old')
      open(13,file='vbas',form='unformatted',status='old')
      open(unit=14,file='bndmat',form='unformatted')
      open(unit=8,file='dipbf',form='unformatted')
c
c      call second(time)
      write(6,119)
119   format(' bound-free dipole matrix element code',//)
      ipow(0)=0
      do 444 i=1,ltop
  444 ipow(i)=1
c
c read in spline information for free functions
c bessel case only!!!
c
      read(5,888) iword
 888  format(a8)
      if(iword.eq."bessel") then
	open(unit=88,file='bessplr',form='unformatted')
       icp=1
       read(88)llmax,nr,rmin,rdel,alpha
       if(nr.gt.nsblok)then
       write(6,*)' nsblok is too small',nsblok,nr
       call exit
       endif
       if(llmax.gt.ltop)then
       write(6,*)' ltop is too small',ltop,llmax
       call exit
       endif
c complex io equivalenced to real io
      nrr=2*nr
       read(88)(x(i),i=1,nr)
       read(88)((hsr(i,k,icp),i=1,nrr),k=0,llmax)
       read(88)((hsder(i,k,icp),i=1,nr),k=0,llmax)
       read(88)((cjr(i,k,icp),i=1,nrr),k=0,llmax)
       read(88)((cy(i,k,icp),i=1,nr),k=0,llmax)
       rd26= rdel* rdel/6.
       rmax=(nr-1)*rdel+rmin
       write(6,666)rmin,rmax,rdel,llmax,alpha
666   format(" the spline points run from ",f10.5," to ",f10.5,
     1 " in steps of ",f10.5/" max l must be less than ",i3/
     2 " test function parameter is :",f6.3)
c
	else
	open(unit=88,file='cesspl',form='unformatted')	
        read(88)
c
      endif
c
c read grid parameters from basis set file
c
c tnr ***********
c      call rdabs(13,mgrid,3,0)
      read(13)mgrid,npt,nbf
c
c read parameters from ylm file and check
c
      read(10)nbig,mpt,lmax,mumax
      if(mgrid.ne.nbig.or.mpt.ne.npt)then
       write(6,701)mgrid,nbig,npt,mpt
701    format("stopping because of mismatch between ylm and basis files"
     1/" mgrid,nbig,npt,mpt :",4i5)
      endif
c      call rdabs(10,nbig,4,0)
c      write(6,*)' nbig,mpt,lmax,mumax from ylms:',
c     $ nbig,mpt,lmax,mumax
      if(lmax.gt.ltop)then
      write(6,*)' ltop is too small',ltop,lmax
      stop
      endif
c tnr ***********
c
c
c read in information about energies and channels
c the scattering energies are incident electron energies
c
      read(5,*) iprint,iorth
      read(14)nchan
      nchan2=nchan
      read(14) (echan(i),i=1,nchan)
      ignd=ismin(nchan,echan,1)
      eground=echan(ignd)
      read(14) nener
      read(14) (energy(i),i=1,nener)
      write(6,202)nchan,(echan(i),i=1,nchan)
202   format(' TARGET energies for ',i3,' channels:',/,(2x,5e15.8))
      write(6,203) nener,(energy(i),i=1,nener)
203   format(1x,i3,' incident energies:',/,(2x,5e15.8))
      read(14) nbftot,nmotot
c
      if(nbftot.ne.nbf)then
         write(6,*)' basis number inconsistency', nbf,nbftot
         stop
      endif
      if(nbftot.gt.nbfmax) then
      write(6,*)' nbftot gt nbfmax ',nbftot,nbfmax
      stop 11
      end if
c..bhl      nmotot=nbftot
      nsy = nbftot
      call rdbinsqr(tr,nbfmax,nbftot,14)
      write(6,99)
99    format(" transformation matrix")
      do 1 i=1,nmotot
      write(6,101)i,(tr(j,i),j=1,nbftot)
1     continue
 101  format(i5,8e14.7/(5x,8e14.7))
c
c for absolute files:
c  jwhere is location in ylm file
c  iset is location in r grid file
c
c tnr***********
c      call rdabs(9,ngrid,1,0)
c      if(mgrid.ne.ngrid)then
c        write(6,700)mgrid,ngrid
c       stop
c700    format(" stopping because of grid mismatch, ngrid,mgrid:",2i5)
c      endif
c tnr***********
c
c read channel quantum numbers
c
      do 22 ic=1,nchan2
      read(5,*) nlm(ic)
       if(nlm(ic).gt.lmtop)then
       write(6,*)' lmtop is too small',lmtop,nlm(ic)
       call exit
       endif
      nlmic=nlm(ic)
      read(5,*) (lch(j,ic),mch(j,ic),j=1,nlmic)
      write(6,752) ic
  752 format(' l s and m s for channel:',i5)
      write(6,'(5x,2i3)') (lch(j,ic),mch(j,ic),j=1,nlmic)
      do 23 j=1,nlmic
      m=2*iabs(mch(j,ic))
      if(mch(j,ic).lt.0) m = m-1
      mch(j,ic)=m
23    continue
      read(5,*)nscat(ic)
      write(6,138) (nscat(ic))
138   format(' number of scattering orbitals in this channel:',i3)
      read(5,*) (nsch(j,ic),j=1,nscat(ic))
      write(6,139) ic, (nsch(j,ic),j=1,nscat(ic))
139   format(' channel',i3,' orbs:',20i3,/,(17x,20i3))
      read(5,*)ntarg(ic)
      write(6,137) ntarg(ic)
137   format(' number of TARGET orbitals in this channel:',6i3)
      if(ntarg(ic).ne.0)then
      read(5,*) (ntch(j,ic),j=1,ntarg(ic)) 
      write(6,139) ic, (ntch(j,ic),j=1,ntarg(ic))
      endif
c
22      continue
      nbtot=0
      mfree=0
      do 44 ic=1,nchan
      nbtot=nbtot+nscat(ic)
      mfree=mfree+nlm(ic)
      ntic=ntarg(ic)
      nsic=nscat(ic)
      nbas(ic)=ntic+nsic
      do 223 i=1,nsic
 223     nbch(i,ic)=nsch(i,ic)
      if (ntic.eq.0) go to 44
      do 333 i=1,ntic
333     nbch(nsic+i,ic) = ntch(i,ic)
44     continue
      ngauss=nbf
      write(6,753)nbf
  753 format('number of ao basis functions=  ',i5)
c
c write headers for output integrals files
c
c  bound-free:
      write(8) nener,nchan,(nlm(ic),ic=1,nchan),
     >         (nscat(ic),ic=1,nchan)
      write(8) ((lch(j,ic),mch(j,ic),j=1,nlm(ic)),ic=1,nchan)
      write(8) ((nsch(j,ic),j=1,nscat(ic)),ic=1,nchan)
      write(8) eground
      write(8) nmotot
        do 29 j=1,nbfmax
	do 29 k=1,nbfmax
        hbx(j,k)=0.0 
	hby(j,k)=0.0 
	hbz(j,k)=0.0
        hbtx(j,k)=0.0 
 29   continue
c
c open a loop on energies
c
      do 1000 iene=1,nener
c
c construct channel momenta
c
      do 17 ichan=1,nchan
      ec = energy(iene) - (echan(ichan)-echan(ignd))
      if(ec.le.0.0) then
         iclosed(ichan)=1
         kchan(ichan) = sqrt(-2.0*ec)
      else
         iclosed(ichan)=0
         kchan(ichan) = sqrt(2.0*ec)
      endif
17    kchan(ichan) = sqrt(2.0*ec)
      write(6,*)'iclosed',(iclosed(i),i=1,nchan)
      write(6,717) energy(iene), (kchan(i),i=1,nchan)
717   format(//,' incident E = ',f10.6,/,' channel momenta = ',(6e12.5))
c
c reset all energy-independent file pointers
c
c tnr************
c      jset=3
c      jwhere=4
c tnr************
      iset=1
      iread=0
      iquit=0
c
c initialize bound-free matrix elements
c
      do 25 j=1,nbfmax
      do 25 i=1,lmtop
      do 26 k=1,nchnl
        hpvbx(i,j,k)=0.0
        hpvby(i,j,k)=0.0
        hpvbz(i,j,k)=0.0
        hpvb(i,j,k)=0.
        hpvbtx(i,j,k)=0.0
        hpvbty(i,j,k)=0.0
        hpvbtz(i,j,k)=0.0
26    hpvbt(i,j,k)=0.
25    continue
c
c read in a block of grid points and transfer to a temporary location
c
      marg=min0(ngrid,npt)
      nread=4*marg
      call rdabs(9,buff(1),nread,iset)
c tnr***********
c      iset=iset+nread
      iset=iset+1
c tnr***********
      iread=iread+marg
      ibes=0
 32   continue
33    continue
      narg=marg
      call dcopy(narg,buff(1),4,grid(1,1),1)
      call dcopy(narg,buff(2),4,grid(1,2),1)
      call dcopy(narg,buff(3),4,grid(1,3),1)
       call dcopy(narg,buff(4),4,wt(1),1)
      iremn=ngrid-iread
      if(iremn.eq.0)then
        iquit=1
        go to 34
      endif
      marg=min0(iremn,npt)
      nread=4*marg
      call rdabs(9,buff(1),nread,iset)
c tnr***********
c      iset=iset+nread
      iset=iset+1
c tnr***********
      iread=iread+marg
34    continue
c
c read in a block of gaussians
c
c tnr************
c      call rdabs(13,basis,nbf*narg,jset)
c      jset=jset+nbf*narg
      read(13)(basis(ib),ib=1,nbf*narg)
c tnr************
c
c skip over basis function second derivatives
c
c tnr***********
c      jset=jset+nbf*narg
      read(13)
c tnr***********
c
c read in a block of ylm's
c
      DO 100 m=0,mumax
      IF(m.EQ.0)THEN
      jbuf=narg*(lmax+1)
      ELSE
      jbuf=2*narg*(lmax-m+1)
      ENDIF
c tnr************
      READ(10)(biff(ib),ib=1,jbuf)
c      CALL rdabs(10,biff(1),jbuf,jwhere)
c      jwhere=jwhere+jbuf
c tnr************
      IF(m.EQ.0)THEN
      DO 10 l=0,lmax
      DO 10 i=1,narg
10    ylm(i,l,0)=biff(i+narg*l)
      ELSE
      DO 11 l=m,lmax
      DO 11 i=1,narg
      ylm(i,l,2*m-1)=biff(i+(l-m)*2*narg)
      ylm(i,l,2*m) = biff(i+(l-m)*2*narg +narg)
11    CONTINUE
      ENDIF
100   CONTINUE
c
c compute bessel functions for this BLOCK
c
      DO 12 i = 1,narg
  12  rvec(i)=SQRT(grid(i,1)**2 + grid(i,2)**2 +grid(i,3)**2)
      DO 800 ic=1,nchan
      icp=1
      IF(iword.EQ."coulomb") THEN
      icp=ic
      IF(ibes.EQ.0) THEN
         READ(88) ak,znuc
         IF(iclosed(ic).EQ.0) THEN
            READ(88)llmax,nr,rmin,rdel,alpha
            IF(nr.GT.nsblok)THEN
               WRITE(6,*)' nsblok is too small',nsblok,nr
               CALL EXIT
            ENDIF
            IF(llmax.GT.ltop)THEN
               WRITE(6,*)' ltop is too small',ltop,llmax
               CALL EXIT
            ENDIF
            WRITE(6,342) ak,kchan(ic)
 342        FORMAT(' momenta from tape and computed',2e12.3)
c COMPLEX io equivalenced to REAL io
            nrr=2*nr
            READ(88)(x(i),i=1,nr)
            READ(88)((hsr(i,k,icp),i=1,nrr),k=0,llmax)
            READ(88)((hsder(i,k,icp),i=1,nr),k=0,llmax)
            READ(88)((cjr(i,k,icp),i=1,nrr),k=0,llmax)
            READ(88)((cy(i,k,icp),i=1,nr),k=0,llmax)
            rd26= rdel* rdel/6.
            rmax=(nr-1)*rdel+rmin
            WRITE(6,666)rmin,rmax,rdel,llmax,alpha
         ENDIF
         IF(lmax.GT.llmax) THEN
            WRITE(6,702)
 702  FORMAT(' stopping because lmax gt maximum l in bessel splines',//)
            STOP
         ENDIF
      ENDIF
      ENDIF
      DO 121 i=1,narg
121   krvec(i)=kchan(ic)*rvec(i)
      aksq = SQRT(kchan(ic))
      ak52 = aksq*kchan(ic)**2
      IF(iclosed(ic).EQ.0) THEN
      DO 750 l=0,lmax
      DO 13 i=1,narg
c********************
c the expression for klo is split into two fortran
c statements so that klo is never less than 1
c*****************
      klo=(krvec(i)-rmin)/rdel
      klo=klo+1
      a=(x(klo+1)-krvec(i))/rdel
      b=(krvec(i)-x(klo))/rdel
      hp(i,l,ic)=a*hs(klo,l,icp)+b*hs(klo+1,l,icp)+(a*(a*a-1.)
     1 *cj(klo,l,icp)+b*(b*b-1.)*cj(klo+1,l,icp))*rd26
      hp(i,l,ic)=hp(i,l,ic)*aksq
13    CONTINUE
750   CONTINUE
      ENDIF
c CLOSE loop on channels
  800 CONTINUE
c
c-------------------------------------------------------------
c  at this point a BLOCK of everything on the grid has been
c  READ in.  now CALL routines which accumulate the bound-free
c  dipole matrix elements from contributions at
c  each grid point
c-------------------------------------------------------------
c
      CALL cdipint(hpvb,hpvbx,hpvby,hpvbz,nchan2,nlm,lch,mch
     1  ,nchnl,lmtop,ngauss,nbfmax,cvec,dvec,basis,grid,wt,
     2  ylm,hp,nblok,narg,ltop,hbx,hby,hbz,iene,iclosed)
c
c RETURN to start a fetch another BLOCK of points
c
      ibes=1
      IF(iquit.EQ.0)go to 32
c
c transform to mo basis
c
      IF(iprint.NE.0) THEN
      IF(iene.EQ.1)THEN
      WRITE(6,*)'x'
      DO i=1,ngauss
         WRITE(6,200)i,(hbx(i,j),j=1,ngauss)
      ENDDO
      WRITE(6,*)'y'
      DO i=1,ngauss
         WRITE(6,200)i,(hby(i,j),j=1,ngauss)
      ENDDO
      WRITE(6,*)'z'
      DO i=1,ngauss
         WRITE(6,200)i,(hbz(i,j),j=1,ngauss)
      ENDDO
 200  FORMAT(i5,10e12.4/(5x,10e12.4))
      ENDIF
      ENDIF
      CALL trnmat(hpvbx,hpvby,hpvbz,hpvbtx,hpvbty,hpvbtz,
     $hpvb,hpvbt,hbx,hby,hbz,hbtx,tr,nchan2,nlm,
c     $ngauss,nmotot,nchnl,nbfmax,lm2,iene)
     $ngauss,nmotot,nchnl,nbfmax,lm2,iene)
c-----------------------------------------------------------
c output section
c--------------------------------------------------------
c
c WRITE bound-free matrices to printer IF requested
c
      IF(iprint.NE.0) THEN
      DO 603 ic=1,nchan2
      WRITE(6,757) ic
757   FORMAT(//,' bound-free dipole matrixs for channel no.:',i3)
      nlmic = nlm(ic)
      ngjc = nmotot
        WRITE(6,576)
576     FORMAT( ' x ')
      DO 503 ilm = 1,nlmic
      WRITE(6,756) ilm,(hpvbtx(ilm,ig,ic),ig=1,ngjc)
756   FORMAT(1x,i3,8e12.5,/,(4x,8e12.5))
503   CONTINUE
        WRITE(6,577)
577     FORMAT( ' y ')
      DO 513 ilm = 1,nlmic
      WRITE(6,756) ilm,(hpvbty(ilm,ig,ic),ig=1,ngjc)
513   CONTINUE
        WRITE(6,578)
578     FORMAT( ' z ')
      DO 523 ilm = 1,nlmic
      WRITE(6,756) ilm,(hpvbtz(ilm,ig,ic),ig=1,ngjc)
523   CONTINUE
        WRITE(6,579)
579     FORMAT( ' bound-free ')
      DO 533 ilm = 1,nlmic
      WRITE(6,756) ilm,(hpvbt(ilm,ig,ic),ig=1,ngjc)
533   CONTINUE
  603 CONTINUE
      ENDIF
c
c orthogonalize free functions to bound functions
c     
      IF(iorth.EQ.1)THEN
         WRITE(6,*)' calling bforthog'
      CALL bforthog(hpvbt,lmtop,nbfmax,nchnl,nstate,hpvbtx,
     1hbx, nchan,nlm,nbas,nmotot,iprint,nbch)
      CALL bforthog(hpvbt,lmtop,nbfmax,nchnl,nstate,hpvbty,
     1hby, nchan,nlm,nbas,nmotot,iprint,nbch)
      CALL bforthog(hpvbt,lmtop,nbfmax,nchnl,nstate,hpvbtz,
     1     hbz, nchan,nlm,nbas,nmotot,iprint,nbch)
      ENDIF
      IF(iprint.NE.0) THEN
      IF(iene.EQ.1)THEN
        WRITE(6,576)
      DO 633 jg = 1,nmotot
      WRITE(6,756) jg,(hbx(jg,ig),ig=1,nmotot)
633   CONTINUE
        WRITE(6,577)
      DO 613 jg = 1,nmotot
      WRITE(6,756) jg,(hby(jg,ig),ig=1,nmotot)
613   CONTINUE
        WRITE(6,578)
      DO 623 jg = 1,nmotot
      WRITE(6,756) jg,(hbz(jg,ig),ig=1,nmotot)
623   CONTINUE
c
      ENDIF
      ENDIF
c
c WRITE bound-bound matrices to disk
c
      IF(iene.EQ.1) THEN
      WRITE(8) ((hbx(ig,jg),ig=1,nmotot),jg=1,nmotot)
      WRITE(8)((hby(ig,jg),ig=1,nmotot),jg=1,nmotot)
      WRITE(8)((hbz(ig,jg),ig=1,nmotot),jg=1,nmotot)
      ENDIF
c
c WRITE OPEN channel bound-free matrices to disk
c
      WRITE(8) (kchan(ic),ic=1,nchan2)
      WRITE(8) (iclosed(ic),ic=1,nchan2)
      nopen=0
      DO i=1,nchan2
         IF(iclosed(i).EQ.0)nopen=nopen+1
      ENDDO
      DO 504 ic=1,nopen
      nlmic=nlm(ic)
      ngjc=nmotot
      WRITE(8)((hpvbtx(il,ig,ic),il=1,nlmic),ig=1,ngjc)
      WRITE(8)((hpvbty(il,ig,ic),il=1,nlmic),ig=1,ngjc)
      WRITE(8)((hpvbtz(il,ig,ic),il=1,nlmic),ig=1,ngjc)
504   CONTINUE
c  tnr*************
c
c reset pointers
c
      REWIND 13
      READ(13)
      REWIND 10
      READ(10)
c  tnr*************
c
c CLOSE big loop on incident energies
c
 1000 CONTINUE
      CALL closeabs(9)
c      CALL second(tnow)
c      elapse=tnow-time
c      WRITE(6,777)elapse
c777   FORMAT(" elapsed time is ",e12.4)
      STOP
      END
