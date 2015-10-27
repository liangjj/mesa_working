      program ffbf
c
c computes free-free and  bound-free matrix elements
c
c created by combining the original bound-free and free codes 11-30-88
c
c  added read of splines and onboard bessel fcn computation 11-30-88
c
c  changed to handle green's function type continuum functions
c     also changed by TNR on 8/8/91 to eliminate unnecessary complex
c      arithmetic in ffint and bfint
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
c changed 5-8-90 nbfmax is now only used to specify the array basis,
c nbfcmx is the maximum no. of basis functions in any one channel.
c also lmtop, the max. no. of l,m pairs, can be adjusted to suit need
c
c parameters:
c
      parameter (nblok=500,mtop=6,ltop=6,nchnl=7)
      parameter (nsblok=3000)
c     parameter (lmtop=(2*mtop+1)*(ltop+1)-mtop*(mtop+1))
      parameter (lmtop=16)
      parameter (nstate=nchnl**2)
      parameter (nsymst = nchnl*(nchnl+1)/2)
      parameter (nbfmax=500)
      parameter (nbfcmx =74)
      parameter (maxene=200)
      complex cvec,dvec,evec,hp,ovbf,hpvb,hpvhp,hpvhm
      common ovbf(lmtop,nbfcmx,nchnl),hpvb(lmtop,nbfcmx,nstate)
     1 , hpvhp(lmtop,lmtop,nsymst),hpvhm(lmtop,lmtop,nchnl**2)
     2 , basis(nbfmax*nblok),cvec(nblok),dvec(nblok),evec(nblok)
     3, hp(nblok,0:ltop,nchnl), hd(nblok,0:ltop,nchnl)
     4, vbuf(nblok*nsymst),vpot(nblok,nsymst),wt(nblok)
     5 ,ylm(nblok,0:ltop,0:2*ltop), ngauss(nchnl),ngch(nbfcmx,nchnl)
     6 , nlm(nchnl),lch(lmtop,nchnl),mch(lmtop,nchnl)
c
c  nblok = number of points in block of grid etc. values
c  mtop = largest m value
c  ltop = largest l value
c  nchnl = largest number of channels allowed
c  nbfmax  =  largest number of l**2 basis functions allowed
c      (transformation to MO's occurs in scattering code)
c  nsblok = max number of spline points for bessel fcns
c  maxene = maximum # of incident energies
c
      complex cj(nsblok,0:ltop)
      complex hs(nsblok,0:ltop)
      real cjr(2*nsblok,0:ltop), cy(nsblok,0:ltop)
      real hsr(2*nsblok,0:ltop), hsder(nsblok,0:ltop)
      real rvec(nblok),krvec(nblok), x(nsblok)
      integer ipow(0:ltop)
      integer smatrix, tmatrix, bcondx
      real grid(nblok,3)
      real echan(nchnl), energy(maxene)
      real kchan(nchnl)
      real buff(nblok*4),biff(2*(nblok*(ltop+1)))
      complex cdotu
      common/ivparms/ istate,ngrid3
      common/parms/mgrid,npt,nbf
      common/parmbes/ngrid2,nchunk,nener,nchan,lmax2
      common/parm/nbig,mpt,lmax,mumax
      data pi/3.14159265358979/
      data smatrix/8hsmatrix  /,tmatrix /8htmatrix  /
      equivalence (cj(1),cjr(1)),(hs(1),hsr(1))
c
c all io units defined here (no calls to "create")
c
c      call link("unit6=(outffbf,create,hc),unit9=(grid,abs,open)
c     1 ,unit12=(vstat,abs,open),unit10=(ylms,abs,open),
c     # unit14=bndmat,
c     2 unit13=(vbas,abs,open)
c     3 ,unit8=intsbf,unit5=inffbf,unit88=bessplr,unit7=intsff//")
c
      open(5,file='inffbf')
      open(6,file='outffbf')
      open(7,file='intsff',form='unformatted')
      open(8,file='intsbf',form='unformatted')
      open(88,file='bessplr',form='unformatted',status='old')
      call openabs(9,'grid')
      call openabs(10,'ylms')
      call openabs(12,'vstat')
      call openabs(13,'vbas')
      open(14,file='bndmat',form='unformatted',status='old')
c
      call second(time)
      write(6,119)
119   format(' bound-free and free-free integrals code',//)
      ipow(0)=0
      do 444 i=1,ltop
  444 ipow(i)=1
c
c read in spline information for free functions
c
       read(88)llmax,nr,rmin,rdel,alpha
c complex io equivalenced to real io
      nrr=2*nr
       read(88)(x(i),i=1,nr)
       read(88)((hsr(i,k),i=1,nrr),k=0,llmax)
       read(88)((hsder(i,k),i=1,nr),k=0,llmax)
       read(88)((cjr(i,k),i=1,nrr),k=0,llmax)
       read(88)((cy(i,k),i=1,nr),k=0,llmax)
c********************
c       write(6,840)(cjr(k,0),k=1,nrr)
c 840   format(10e12.4)
c*******************
       rd26= rdel* rdel/6.
       rmax=(nr-1)*rdel+rmin
       write(6,666)rmin,rmax,rdel,llmax,alpha
666   format(" the spline points run from ",f10.5," to ",f10.5,
     1 " in steps of ",f10.5/" max l must be less than ",i3/
     2 " cutoff parameter is:",f6.3)
c
c read grid parameters from basis set file
c
      call rdabs(13,mgrid,3,0)
c..unicos
c      call waitfor(13)
c..unicos
c
c read parameters from ylm file and check
c
      call rdabs(10,nbig,4,0)
c..unicos
c      call waitfor(10)
c..unicos
      if(mgrid.ne.nbig.or.mpt.ne.npt)then
       write(6,701)mgrid,nbig,npt,mpt
701    format("stopping because of mismatch between ylm and basis files"
     1/" mgrid,nbig,npt,mpt :",4i5)
      endif
      if(lmax.gt.llmax) then
      write(6,702)
702   format(' stopping because lmax gt maximum l in bessel splines',//)
      stop
      endif
c
c
c read in information about energies and channels
c the scattering energies are incident electron energies
c
      read(5,*) iprint
      read(14) nchan
      read(14) (echan(i),i=1,nchan)
      nchan2=nchan
      ignd=ismin(nchan,echan,1)
      read(14) nener
      read(14) (energy(i),i=1,nener)
      write(6,202)nchan,(echan(i),i=1,nchan)
202   format(' target energies for ',i3,' channels:',/,(2x,5e15.8))
      write(6,203) nener,(energy(i),i=1,nener)
203   format(1x,i3,' incident energies:',/,(2x,5e15.8))
c
c for absolute files:
c  jwhere is location in ylm file
c  iset is location in r grid file
c
      call rdabs(9,ngrid,1,0)
c..unicos
c      call waitfor(9)
c..unicos
      if(mgrid.ne.ngrid)then
        write(6,700)mgrid,ngrid
       stop
700    format(" stopping because of grid mismatch, ngrid,mgrid:",2i5)
      endif
c
c read parameters from static potential file and check
c
      call rdabs(12,istate,2,0)
c..unicos
c      call waitfor(12)
c..unicos
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
c
c construct channel momenta
c
      do 17 ichan=1,nchan
      ec = energy(iene) - (echan(ichan)-echan(ignd))
      if(ec.le.0.0) then
          write(6,207) ec
  207     format(' *** no closed channels allowed (yet) ****',e15.8)
      stop
      endif
17    kchan(ichan) = sqrt(2.0*ec)
      write(6,717) energy(iene), (kchan(i),i=1,nchan)
717   format(//,' incident E = ',f10.6,/,' channel momenta = ',(6e12.5))
c
c reset all energy-independent file pointers
c
      jset=3
      jwhere=4
      lwhere = 2
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
      do 25 j=1,nbfcmx
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
      iset=iset+nread
      iread=iread+marg
c..unicos
 32    continue  
c32    if(iowait(9))32,33
c33    continue
c..unicos
      narg=marg
      call scopy(narg,buff(1),4,grid(1,1),1)
      call scopy(narg,buff(2),4,grid(1,2),1)
      call scopy(narg,buff(3),4,grid(1,3),1)
       call scopy(narg,buff(4),4,wt(1),1)
      iremn=ngrid-iread
      if(iremn.eq.0)then
        iquit=1
        go to 34
      endif
      marg=min0(iremn,npt)
      nread=4*marg
      call rdabs(9,buff(1),nread,iset)
      iset=iset+nread
      iread=iread+marg
34    continue
c
c read in a block of gaussians
c
      call rdabs(13,basis,nbf*narg,jset)
      jset=jset+nbf*narg
c
c skip over basis function second derivatives
c
      jset=jset+nbf*narg
c..unicos
c      call waitfor(13)
c..unicos
c
c read in a block of ylm's
c
      do 100 m=0,mumax
      if(m.eq.0)then
      jbuf=narg*(lmax+1)
      else
      jbuf=2*narg*(lmax-m+1)
      endif
      call rdabs(10,biff(1),jbuf,jwhere)
      jwhere=jwhere+jbuf
c..unicos
c      call waitfor(10)
c..unicos
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
      do 121 i=1,narg
121   krvec(i)=kchan(ic)*rvec(i)
      aksq = sqrt(kchan(ic))
      ak52 = aksq*kchan(ic)**2
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
      hp(i,l,ic)=a*hs(klo,l)+b*hs(klo+1,l)+(a*(a*a-1.)*cj(klo,l)+b*
     1 (b*b-1.)*cj(klo+1,l))*rd26
      hp(i,l,ic)=hp(i,l,ic)*aksq
      hd(i,l,ic)=a*hsder(klo,l)+b*hsder(klo+1,l)+(a*(a*a-1.)*cy(klo,l)
     1 +b*(b*b-1.)*cy(klo+1,l))*rd26
      hd(i,l,ic)=hd(i,l,ic)*ak52/krvec(i)**(1-ipow(l))
13    continue
750   continue
c close loop on channels
  800 continue
c
c  read static potentials for this block
c  note ******  weights are already multiplied into potentials *****
c
      call rdabs(12,vbuf(1),narg*istate,lwhere)
      lwhere=lwhere+narg*istate
c..unicos
c      call waitfor(12)
c..unicos
      do 607 iii=1,istate
      idum=(iii-1)*narg+1
      call scopy(narg,vbuf(idum),1,vpot(1,iii),1)
607   continue
c
c-------------------------------------------------------------
c  at this point a block of everything on the grid has been
c  read in.  now call routines which accumulate the bound-free
c  and free-free matrix elements from contributions at
c  each grid point
c-------------------------------------------------------------
c
      call ffints(nchan2,iflag,narg)
      call bfints(nchan2,narg)
c
c return to start a fetch another block of points
c
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
c close big loop on incident energies
c
 1000 continue
      call second(tnow)
      elapse=tnow-time
      write(6,777)elapse
777   format(" elapsed time is ",e12.4)
      call exit
      end
      subroutine ffints(nchan2,iflag,narg)
      parameter (nblok=500,mtop=6,ltop=6,nchnl=7)
      parameter (nsblok=3000)
c     parameter (lmtop=(2*mtop+1)*(ltop+1)-mtop*(mtop+1))
      parameter (lmtop=16)
      parameter (nstate=nchnl**2)
      parameter (nsymst = nchnl*(nchnl+1)/2)
      parameter (nbfmax=500)
      parameter (nbfcmx =74)
      parameter (maxene=200)
      complex cvec,dvec,evec,hp,ovbf,hpvb,hpvhp,hpvhm,ai,toboth
      real crvec(2*nblok),rvec(nblok)
      common ovbf(lmtop,nbfcmx,nchnl),hpvb(lmtop,nbfcmx,nstate)
     1 , hpvhp(lmtop,lmtop,nsymst),hpvhm(lmtop,lmtop,nchnl**2)
     2 , basis(nbfmax*nblok),cvec(nblok),dvec(nblok),evec(nblok)
     3, hp(nblok,0:ltop,nchnl), hd(nblok,0:ltop,nchnl)
     4, vbuf(nblok*nsymst),vpot(nblok,nsymst),wt(nblok)
     5 ,ylm(nblok,0:ltop,0:2*ltop), ngauss(nchnl),ngch(nbfcmx,nchnl)
     6 , nlm(nchnl),lch(lmtop,nchnl),mch(lmtop,nchnl)
      equivalence (cvec(1),crvec(1)),(rvec(1),dvec(1))
c
c  free-free matrices
c
      complex cdotu,cdotc
c
c compute potential matrices
c
      ai=cmplx(0.,1.)
      do 500 ic=1,nchan2
      nlmic=nlm(ic)
      do 500 ilm=1,nlmic
      l1=lch(ilm,ic)
      m1=mch(ilm,ic)
      do 400 jc=1,ic
      fac=1.
      if(ic.eq.jc)fac=.5
      ist = ic*(ic-1)/2 + jc
      icjc=nchan2*(ic-1)+jc
      jcic=nchan2*(jc-1)+ic
      nlmjc=nlm(jc)
      do 400 jlm=1,nlmjc
      l2=lch(jlm,jc)
      m2=mch(jlm,jc)
      if(iflag.eq.0) then
      do 300 i=1,narg
      potval=vpot(i,ist)*ylm(i,l1,m1)*ylm(i,l2,m2)
      cvec(i)=potval*hp(i,l1,ic)
      dvec(i)=fac*potval*conjg(hp(i,l2,jc))
      evec(i)=fac*potval*conjg(hp(i,l1,ic))
300   continue
      else
      do 299 i=1,narg
      potval=vpot(i,ist)*ylm(i,l1,m1)*ylm(i,l2,m2)
      cvec(i)=potval*hp(i,l1,ic)
      dvec(i)=fac*potval*aimag(hp(i,l2,jc))
      evec(i)=fac*potval*aimag(hp(i,l1,ic))
299   continue
      endif
      hpvhp(ilm,jlm,ist) = hpvhp(ilm,jlm,ist) +
     1 cdotu(narg,cvec,1,hp(1,l2,jc),1)
      hpvhm(ilm,jlm,icjc) = hpvhm(ilm,jlm,icjc) +
     1 cdotu(narg,hp(1,l1,ic),1,dvec,1)
      hpvhm(jlm,ilm,jcic) = hpvhm(jlm,ilm,jcic) +
     1 cdotu(narg,hp(1,l2,jc),1,evec,1)
400   continue
500   continue
c
c compute matrix of (kchan**2/2 - T)
c
      do 501 ic=1,nchan2
      nlmic=nlm(ic)
      icc=ic*(ic+1)/2
      iic=nchan2*(ic-1)+ic
      do 501 ilm=1,nlmic
      l1=lch(ilm,ic)
      m1=mch(ilm,ic)
      if(iflag.eq.0) then
      do 301 i=1,narg
      cvec(i)=wt(i)*ylm(i,l1,m1)**2*hp(i,l1,ic)
      rvec(i) = hd(i,l1,ic)
301   continue
      toboth= sdot(narg,rvec,1,crvec,1)+
     2 sdot(narg,rvec,1,crvec(2),2)*ai
      hpvhp(ilm,ilm,icc)=hpvhp(ilm,ilm,icc)+
     1    toboth
      hpvhm(ilm,ilm,iic)=hpvhm(ilm,ilm,iic)+
     1   toboth
      else
      do 302 i=1,narg
      cvec(i)=wt(i)*ylm(i,l1,m1)**2*hp(i,l1,ic)
302   continue
      hpvhp(ilm,ilm,icc)=hpvhp(ilm,ilm,icc)+
     1     sdot(narg,crvec,2,hd(1,l1,ic),1)+
     2 sdot(narg,crvec(2),2,hd(1,l1,ic),1)*ai
      endif
501   continue
      return
      end
      subroutine bfints(nchan2,narg)
      parameter (nblok=500,mtop=6,ltop=6,nchnl=7)
      parameter (nsblok=3000)
c     parameter (lmtop=(2*mtop+1)*(ltop+1)-mtop*(mtop+1))
      parameter (lmtop=16)
      parameter (nstate=nchnl**2)
      parameter (nsymst = nchnl*(nchnl+1)/2)
      parameter (nbfmax=500)
      parameter (nbfcmx =74)
      parameter (maxene=200)
      complex cvec,dvec,evec,hp,ovbf,hpvb,hpvhp,hpvhm,ai
      real crvec(2*nblok)
      common ovbf(lmtop,nbfcmx,nchnl),hpvb(lmtop,nbfcmx,nstate)
     1 , hpvhp(lmtop,lmtop,nsymst),hpvhm(lmtop,lmtop,nchnl**2)
     2 , basis(nbfmax*nblok),cvec(nblok),dvec(nblok),evec(nblok)
     3, hp(nblok,0:ltop,nchnl), hd(nblok,0:ltop,nchnl)
     4, vbuf(nblok*nsymst),vpot(nblok,nsymst),wt(nblok)
     5 ,ylm(nblok,0:ltop,0:2*ltop), ngauss(nchnl),ngch(nbfcmx,nchnl)
     6 , nlm(nchnl),lch(lmtop,nchnl),mch(lmtop,nchnl)
      equivalence (cvec(1),crvec(1))
c
c  bound-free matrices
c
      complex cdotu
c
c compute potential matrix
c
      ai=cmplx(0.,1.)
      do 500 ic=1,nchan2
      nlmic=nlm(ic)
      do 500 ilm=1,nlmic
      l1=lch(ilm,ic)
      m1=mch(ilm,ic)
      do 400 jc=1,nchan2
      ist=nchan2*(ic-1)+jc
      ii=max0(ic,jc)
      jj=ic+jc-ii
      ipot=narg*(ii*(ii-1)/2+jj-1)
      ngjc = ngauss(jc)
      do 400 ig=1,ngjc
      nbch = ngch(ig,jc)
      do 300 i=1,narg
      cvec(i)=ylm(i,l1,m1)*basis(narg*(nbch-1)+i)*vbuf(i+ipot)
300   continue
      hpvb(ilm,ig,ist)=cdotu(narg,cvec,1,hp(1,l1,ic),1)
     1 + hpvb(ilm,ig,ist)
400   continue
500   continue
c
c compute matrix of (E - T)
c
      do 501 ic=1,nchan2
      icc=nchan2*(ic-1)+ic
      nlmic=nlm(ic)
      ngic = ngauss(ic)
      do 501 ilm=1,nlmic
      l1=lch(ilm,ic)
      m1=mch(ilm,ic)
      do 401 ig=1,ngic
      nbch = ngch(ig,ic)
      do 301 i=1,narg
      cvec(i)=wt(i)*ylm(i,l1,m1)*basis(narg*(nbch-1)+i)
301   continue
      ovbf(ilm,ig,ic)=cdotu(narg,hp(1,l1,ic),1,cvec,1)
     1 + ovbf(ilm,ig,ic)
      hpvb(ilm,ig,icc)=hpvb(ilm,ig,icc)+sdot(narg,hd(1,l1,ic),1,crvec
     1 ,2)+ai*sdot(narg,hd(1,l1,ic),1,crvec(2),2)
401   continue
501   continue
      return
      end
      subroutine waitfor(godot)
c waits for godot until he is no longer busy
      integer godot
10    if(iowait(godot)) 10,20
20    continue
      return
      end
