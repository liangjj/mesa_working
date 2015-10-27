      program ffbf
      implicit real*8 (a-h,o-z)
c
c computes free-free and  bound-free matrix elements
c
c created by combining the original bound-free and free codes 11-30-88
c
c  added read of splines and onboard bessel fcn computation 11-30-88
c
c  reads file of r grid (transformed)
c  reads file of ylms
c  reads file of spline coefficients for bessel functions computation
c  reads file of L**2 basis functions on grid
c  reads file of static potential on grid
c
c
c  free-free integrals:
c    computes  M = (H - 0.5*kchan**2) matrix elements of
c    two classes: <h+ M h+> and <h+ M h->
c     the definition of h- changes depending on
c    whether T-matrix or S-matrix boundary conditions are specified
c
c  bound-free integrals:
c  computes <h+, (H-E) b>
c
c   note that codes allow for there to be potentials for several states
c
c  changed to read in lch, mch channel quantum numbers 8-2-88
c
c complex io changed to real io because of system error 11-30-88
c
c
c parameters:
c
      parameter (nblok=500,mtop=#maxltop,ltop=#maxltop,nchnl=#maxchan,nchnlc=#maxchan)
      parameter (nsblok=30000)
c     parameter (lmtop=(2*mtop+1)*(ltop+1)-mtop*(mtop+1))
c     parameter (lmtop=(ltop+1)**2)
      parameter (lmtop=#maxlmtop)
      parameter (nstate=nchnl**2)
      parameter (nsymst = nchnl*(nchnl+1)/2)
      parameter (nbfmax=#maxnbfkohn)
      parameter (maxene=200)
c
c  nblok = number of points in block of grid etc. values
c  mtop = largest m value
c  ltop = largest l value
c  nchnl = largest number of channels allowed
c  nchnlc = largest number of channels allowed coulomb case
c  nbfmax  =  largest number of l**2 basis functions allowed
c      (transformation to MO's occurs in scattering code)
c  nsblok = max number of spline points for bessel fcns
c  maxene = maximum # of incident energies
c
      complex*16 hp,hd,cj,hs
      common hp(nblok,0:ltop,nchnl), hd(nblok,0:ltop,nchnl)
     1 , cj(nsblok,0:ltop,nchnlc), cy(nsblok,0:ltop,nchnlc)
     1 , hs(nsblok,0:ltop,nchnlc), hsder(nsblok,0:ltop,nchnlc)
      real*8 cjr(2*nsblok,0:ltop,nchnlc), hsr(2*nsblok,0:ltop,nchnlc)
      real*8 rit(3)
      real*8 rvec(nblok),krvec(nblok), x(nsblok)
      integer ipow(0:ltop),iclosed(nchnl),it(3)
      character*8 smatrix, tmatrix, bcondx,iword,bessel,coulomb
      real*8 grid(nblok,3)
      real*8 vbuf(nblok*nchnl*(nchnl+1)/2)
      real*8 vpot(nblok,nsymst)
      real*8 wt(nblok)
      real*8 echan(nchnl), energy(maxene)
      real*8 kchan(nchnl)
      real*8 ylm(nblok,0:ltop,0:2*ltop)
      real*8 buff(nblok*4),biff(2*(nblok*(ltop+1)))
      real*8 basis(nbfmax*nblok)
      dimension nlm(nchnl),lch(lmtop,nchnl),mch(lmtop,nchnl)
      dimension ngauss(nchnl),ngch(nbfmax,nchnl)
      complex*16 cvec(nblok),hpvb(lmtop,nbfmax,nstate)
      complex*16 ovbf(lmtop,nbfmax,nchnl)
      complex*16 hpvhp(lmtop,lmtop,nsymst),hpvhm(lmtop,lmtop,nchnl**2)
      complex*16 evec(nblok), dvec(nblok)
      common/ivparms/ istate,ngrid3
      common/parms/mgrid,npt,nbf
      common/parmbes/ngrid2,nchunk,nener,nchan,lmax2
      common/parm/nbig,mpt,lmax,mumax
      data pi/3.14159265358979/
      data smatrix/"smatrix"  /,tmatrix /"tmatrix"  /,bessel/"bessel"/
     $ ,coulomb/"coulomb"/
      equivalence (cj(1,0,1),cjr(1,0,1)),(hs(1,0,1),hsr(1,0,1))
c
c all io units defined here (no calls to "create")
c
c      call link("unit6=(outffbf,create,hc),unit9=(grid,abs,open)
c     1 ,unit12=(vstat,abs,open),unit10=(ylms,abs,open),
c     2 unit13=(vbas,abs,open)
c     # ,unit14=(bndmat)
c     3 ,unit8=intsbf,unit5=inffbf,unit88=cesspl,unit7=intsff//")
c
      open(5,file='inffbf',status='unknown')
      open(6,file='outffbf',status='unknown')
      read(5,*)npt
      irecl=4*npt*8 
      call openabs(9,'grid',irecl)
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

      open(7,file='intsff',form='unformatted',status='unknown')
      open(8,file='intsbf',form='unformatted',status='unknown')
      open(88,file='cesspl',form='unformatted',status='old')
      open(13,file='vbas',form='unformatted',status='old')
      open(10,file='ylms',form='unformatted',status='old')
      open(12,file='vstat',form='unformatted',status='old')
      open(14,file='bndmat',form='unformatted',status='old')
      write(6,119)
119   format(' bound-free and free-free integrals code',//)
      ipow(0)=0
      do 444 i=1,ltop
  444 ipow(i)=1
c
c read in spline information for free functions
c bessel case only!!!
c
      read(88) iword
      if(iword.eq.bessel) then
       icp=1
       read(88)llmax,nr,rmin,rdel,alpha
c complex io equivalenced to real io
      nrr=2*nr
       read(88)(x(i),i=1,nr)
       read(88)((hsr(i,k,icp),i=1,nrr),k=0,llmax)
       read(88)((hsder(i,k,icp),i=1,nr),k=0,llmax)
       read(88)((cjr(i,k,icp),i=1,nrr),k=0,llmax)
       read(88)((cy(i,k,icp),i=1,nr),k=0,llmax)
       rd26= rdel* rdel/6.
       rmax=(nr-1)*rdel+rmin
       write(6,666)rmin,rmax,rdel,llmax,alpha,nmax
666   format(" the spline points run from ",f10.5," to ",f10.5,
     1 " in steps of ",f10.5/" max l must be less than ",i3/
     2 " cutoff parameters are alpha,nmax:",f6.3,i3)
      endif
c
c read grid parameters from basis set file
c
      read(13)mgrid,npt,nbf
c
c read parameters from ylm file and check
c
      read(10)nbig,mpt,lmax,mumax
      if(mgrid.ne.nbig.or.mpt.ne.npt)then
       write(6,701)mgrid,nbig,npt,mpt
701    format("stopping because of mismatch between ylm and basis files"
     1/" mgrid,nbig,npt,mpt :",4i5)
      stop
      endif
c
c
c read in information about energies and channels
c the scattering energies are incident electron energies
c
      read(5,*) iprint
      read(14) nchan
      nchan2=nchan
      read(14) (echan(i),i=1,nchan)
      ignd=ismin(nchan,echan,1)
      read(14) nener
      read(14) (energy(i),i=1,nener)
      write(6,202)nchan,(echan(i),i=1,nchan)
202   format(' target energies for ',i3,' channels:',/,(2x,5e15.8))
      write(6,203) nener,(energy(i),i=1,nener)
203   format(1x,i3,' incident energies:',/,(2x,5e15.8))
      if(mgrid.ne.ngrid)then
        write(6,700)mgrid,ngrid
       stop
700    format(" stopping because of grid mismatch, ngrid,mgrid:",2i5)
      endif
c
c read parameters from static potential file and check
c
      read(12)istate,ngrid3
      if(ngrid3.ne.ngrid) then
        write(6,732)
732   format(' stopping because of mismatch between grid and potential')
        call exit
      endif
      write(6,771) istate
771   format(//,' number of coupling blocks',i5)
c
c read boundary condition specification
c
      read(5,734) bcondx
734   format(a8)
      iflag=-1
      if(bcondx.eq.smatrix) then
      iflag=0
      write(6,731)
731   format(//,' S-matrix boundary conditions')
      endif
      if(bcondx.eq.tmatrix) then
      iflag=1
      write(6,733)
733   format(//,' T-matrix boundary conditions')
      endif
      if(iflag.lt.0) then
      write(6,736)
736   format(//,' stopping because of unrecognizable boundary cond.')
      stop
      endif
c read channel quantum numbers and assignment of
c  l**2 basis functions to channels
c
      if(istate.ne.nchan*(nchan+1)/2) then
      write(6,587)
587   format(' stopping because: istate ne nchan*(nchan+1)/2')
      stop
      endif
      do 22 ic=1,nchan2
      read(5,*) nlm(ic)
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
c
      read(5,*) ngauss(ic)
      ngic=ngauss(ic)
      read(5,*) (ngch(j,ic),j=1,ngic)
      write(6,753) ic
  753 format(' basis functions for channel:',i3)
      write(6,'(5x,10i4)') (ngch(j,ic),j=1,ngic)
22    continue
c
c write headers for output integrals files
c
c  bound-free:
      write(8) nener,nchan,(nlm(ic),ic=1,nchan)
      write(8) ((lch(j,ic),mch(j,ic),j=1,nlm(ic)),ic=1,nchan)
      eground=echan(ignd)
      write(8) eground
      write(8) (ngauss(ic),ic=1,nchan)
      write(8) ((ngch(j,ic),j=1,ngauss(ic)),ic=1,nchan)
c  free-free:
      write(7) iflag
      write(7) nener,nchan,(nlm(ic),ic=1,nchan)
      write(7) ((lch(j,ic),mch(j,ic),j=1,nlm(ic)),ic=1,nchan)
      write(7) eground
c
c open a loop on energies
c
      do 1000 iene=1,nener
         ibes=0
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
17    continue
      write(6,717) energy(iene), (kchan(i),i=1,nchan)
717   format(//,' incident E = ',f10.6,/,' channel momenta = ',(6e12.5))
c
c reset all energy-independent file pointers
c
c      jset=3
c      jwhere=4
c      lwhere = 2
      iset=1
      iread=0
      iquit=0
c initialize free-free matrix elements
c
      do 1 i=1,lmtop
      do 1 ist=1,nsymst
      do 1 j=1,lmtop
1     hpvhp(j,i,ist) = 0.0
      nn2=nchnl**2
      do 2 i=1,lmtop
      do 2 ist=1,nn2
      do 2 j=1,lmtop
      hpvhm(j,i,ist)=0.
2     continue
c
c initialize bound-free matrix elements
c
      do 25 i=1,lmtop
      do 25 j=1,nbfmax
      do 26 k=1,nstate
26    hpvb(i,j,k)=0.
      do 27 k=1,nchnl
27    ovbf(i,j,k)=0.
25    continue
c
c read in a block of grid points and transfer to a temporary location
c
      marg=min0(ngrid,npt)
      nread=4*marg
      call rdabs(9,buff(1),nread,iset)
      iset=iset+1
      iread=iread+marg
 32    continue  
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
      iset=iset+1
      iread=iread+marg
34    continue
c
c read in a block of gaussians
c
      read(13)(basis(ib),ib=1,nbf*narg)
c
c skip over basis function second derivatives
c
      read(13)
c
c read in a block of ylm's
c
      do 100 m=0,mumax
      if(m.eq.0)then
      jbuf=narg*(lmax+1)
      else
      jbuf=2*narg*(lmax-m+1)
      endif
      read(10)(biff(ib),ib=1,jbuf)
      if(m.eq.0)then
      do 10 l=0,lmax
      do 10 i=1,narg
10    ylm(i,l,0)=biff(i+narg*l)
      else
      do 11 l=m,lmax
      do 11 i=1,narg
      ylm(i,l,2*m-1)=biff(i+(l-m)*2*narg)
      ylm(i,l,2*m) = biff(i+(l-m)*2*narg +narg)
11    continue
      endif
100   continue
c
c compute bessel functions for this block
c
      do 12 i = 1,narg
  12  rvec(i)=sqrt(grid(i,1)**2 + grid(i,2)**2 +grid(i,3)**2)
      do 800 ic=1,nchan
      icp=1
      if(iword.eq.coulomb) then
      icp=ic
      if(ibes.eq.0) then
         read(88) ak,znuc
      if(iclosed(ic).eq.0) then
       read(88)llmax,nr,rmin,rdel,alpha
      write(6,342) ak,kchan(ic)
342   format(' momenta from tape and computed',2e12.3)
c complex io equivalenced to real io
      nrr=2*nr
       read(88)(x(i),i=1,nr)
       read(88)((hsr(i,k,icp),i=1,nrr),k=0,llmax)
       read(88)((hsder(i,k,icp),i=1,nr),k=0,llmax)
       read(88)((cjr(i,k,icp),i=1,nrr),k=0,llmax)
       read(88)((cy(i,k,icp),i=1,nr),k=0,llmax)
       rd26= rdel* rdel/6.
       rmax=(nr-1)*rdel+rmin
       write(6,666)rmin,rmax,rdel,llmax,alpha,nmax
      endif
      endif
      if(lmax.gt.llmax) then
      write(6,702)
702   format(' stopping because lmax gt maximum l in bessel splines',//)
      stop
      endif
      endif
      do 121 i=1,narg
121   krvec(i)=kchan(ic)*rvec(i)
      aksq = sqrt(kchan(ic))
      ak52 = aksq*kchan(ic)**2
      if(iclosed(ic).eq.0) then
      do 750 l=0,lmax
      do 13 i=1,narg
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
      hd(i,l,ic)=a*hsder(klo,l,icp)+b*hsder(klo+1,l,icp)+(a*(a*a-1.)
     1 *cy(klo,l,icp)+b*(b*b-1.)*cy(klo+1,l,icp))*rd26
      hd(i,l,ic)=hd(i,l,ic)*ak52/krvec(i)**(1-ipow(l))
13    continue
 750  continue
      endif
c close loop on channels
  800 continue
c
c  read static potentials for this block
c  note ******  weights are already multiplied into potentials *****
c
      read(12)(vbuf(ib),ib=1,narg*istate)
      do 607 iii=1,istate
      idum=(iii-1)*narg+1
      call dcopy(narg,vbuf(idum),1,vpot(1,iii),1)
607   continue
c
c-------------------------------------------------------------
c  at this point a block of everything on the grid has been
c  read in.  now call routines which accumulate the bound-free
c  and free-free matrix elements from contributions at
c  each grid point
c-------------------------------------------------------------
c
      call ffints(hpvhp,hpvhm,nchan2,nlm,lch,mch,nchnl,lmtop
     1  ,cvec,dvec,evec,vpot,basis,wt,ylm,hp,hd,nblok
     2  ,iflag,narg,ltop,iclosed)
      call bfints(hpvb,ovbf,nchan2,nlm,lch,mch,nchnl,lmtop
     1  ,ngauss, ngch,nbfmax,cvec,basis,vbuf,wt,ylm,hp,hd,nblok,
     2  nstate,narg,ltop,iclosed)
c
c return to start a fetch another block of points
c
      ibes=1
      if(iquit.eq.0)go to 32
c-----------------------------------------------------------
c output section
c--------------------------------------------------------
c
c write bound-free matrices to printer if requested
c
      if(iprint.ne.0) then
      do 503 ic=1,nchan2
      do 503 jc=1,nchan2
      ist=nchan2*(ic-1)+jc
      write(6,757) ic,jc
757   format(//,' bound-free (h-e) matrix for channels:',2i3)
      nlmic = nlm(ic)
      ngjc = ngauss(jc)
      do 503 ilm = 1,nlmic
      write(6,756) ilm,(hpvb(ilm,ig,ist),ig=1,ngjc)
756   format(1x,i3,8e12.5,/,(4x,8e12.5))
503   continue
      do 505 ic=1,nchan2
      write(6,506)ic
506   format(//,' bound-free overlap matrix for channel',i3)
      nlmic=nlm(ic)
      ngic=ngauss(ic)
      do 505 ilm=1,nlmic
505   write(6,756)ilm,(ovbf(ilm,ig,ic),ig=1,ngic)
      endif
c
c write bound-free matrices to disk
c
      write(8) (kchan(ic),ic=1,nchan2)
      write(8) (iclosed(ic),ic=1,nchan2)
      write(6,*) 'iclosed',(iclosed(ic),ic=1,nchan2)
      do 507 ic=1,nchan2
      nlmic=nlm(ic)
      ngic=ngauss(ic)
      write(8)((ovbf(ilm,ig,ic),ilm=1,nlmic),ig=1,ngic)
507   continue
      do 504 ic=1,nchan2
      do 504 jc=1,nchan2
      ist=nchan2*(ic-1)+jc
      nlmic=nlm(ic)
      ngjc=ngauss(jc)
      write(8)((hpvb(ilm,ig,ist),ilm=1,nlmic),ig=1,ngjc)
504   continue
c
c write free-free matrices to printer if requested
c
      if(iprint.ne.0) then
      do 603 ic=1,nchan2
      do 603 jc=1,ic
      ist=ic*(ic-1)/2 + jc
      nlmic = nlm(ic)
      nlmjc = nlm(jc)
      write(6,663) ic,jc
663   format(//,'h+ h+ (h-e) matrix (ml,ml) for channels:',2i3)
      do 600 ilm=1,nlmic
600   write(6,66) ilm,(hpvhp(ilm,jlm,ist),jlm=1,nlmjc)
66    format(i5,(8e12.4),/,(5x,8e12.4))
603   continue
      do 609 ic=1,nchan2
      do 609 jc=1,nchan2
      ist=nchan2*(ic-1)+jc
      nlmic = nlm(ic)
      nlmjc = nlm(jc)
      write(6,664) ic,jc
664   format(//,'h+ h- (h-e) matrix (ml,ml) for channels:',2i3)
      do 601 ilm=1,nlmic
601   write(6,66) ilm,(hpvhm(ilm,jlm,ist),jlm=1,nlmjc)
609   continue
      endif
c
c write free-free (h-e) matrices to disk
c
      write(7) (kchan(ic),ic=1,nchan2)
      do 604 ic=1,nchan2
      do 604 jc=1,ic
      ist = ic*(ic-1)/2 +jc
      nlmic=nlm(ic)
      nlmjc=nlm(jc)
      write(7) ((hpvhp(ilm,jlm,ist),ilm=1,nlmic),jlm=1,nlmjc)
604   continue
      do 608 ic=1,nchan2
      do 608 jc=1,nchan2
      ist=nchan2*(ic-1)+jc
      nlmic=nlm(ic)
      nlmjc=nlm(jc)
      write(7) ((hpvhm(ilm,jlm,ist),ilm=1,nlmic),jlm=1,nlmjc)
608   continue
c
c reset pointers
c
      rewind 13
      read(13)
      rewind 10
      read(10)
      rewind 12
      read(12)
c
c close big loop on incident energies
c
 1000 continue
      call closeabs(9)
      stop
      end

