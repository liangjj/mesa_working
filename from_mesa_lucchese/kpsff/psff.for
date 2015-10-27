      parameter (mxcen=20,npseudmx=6,lmax=3,ntmx=20,ltop=#maxltop,mtop=#maxltop)
      parameter (nchnl=#maxchan,nsblok=3000,nqmx=500,nqmx2=2*nqmx)
      parameter (limit=300,len=4*limit)
c     parameter (lmtop=(2*mtop+1)*(ltop+1)-mtop*(mtop+1))
      parameter (lmtop=#maxlmtop,maxene=200,nblok=500,nsplmx=500)
      parameter (nbfcmx =130)
      parameter (nbfmax=#maxnbfkohn)
      implicit real*8(a-h,o-z)
      real time, tnow
      dimension iwk(limit),wk(len)
      integer jbasis(npseudmx),icum(0:(lmax-1)),nclm(0:lmax,npseudmx)
     $     ,it(3)
      real*8 xplm(nblok),cphi(nblok),sphi(nblok),cmphi(nblok),
     $     smphi(nblok),p(nblok,0:lmax),pp(nblok,0:lmax)
     $     ,pm(nblok,0:lmax),pscr(7*nblok)
      real*8 wl(10),xl(10),uuu(nblok,lmax*npseudmx)
      real*8 uspl(nsplmx,npseudmx,lmax+1),cspl(nsplmx,npseudmx,lmax+1)
      real*8 rspl(nsplmx),scr(nsplmx),bvec(nblok,ntmx*lmax**2)
      real*8 xx(nblok),yy(nblok),zz(nblok),buff(4*nblok),ww(nblok)
      real*8 ppvr(nqmx2),pmvr(nqmx2),biff(2*(nblok*(ltop+1))),
c     $     barf(2*(lmax-1)*nblok,0:(lmax-1))
     $ barf(2*(lmax-1)*nblok),vec(nblok,lmax**2),gauss(nbfmax*nblok)
      real*8 echan(nchnl), energy(maxene),rrvec(nblok)
      real*8 kchan(nchnl),rr(nqmx),wt(nqmx),u(nqmx,lmax+1)
      real*8 ylm(nblok,0:ltop,0:2*ltop),rvec(nblok),uu(nblok)
      integer nlm(nchnl),lch(lmtop,nchnl),mch(lmtop,nchnl)
     $,iclosed(nchnl),njlm(npseudmx),ngauss(nchnl),ngch(nbfcmx,nchnl)
      real*8 cjr(2*nsblok,0:ltop),hsr(2*nsblok,0:ltop),x(nsblok)
      real*8 ppmatr(2*nblok,lmtop)
      real*8 vecschr(2*npseudmx*ntmx*lmax**2,lmtop,nchnl)
      real*8 vccschr(2*npseudmx*ntmx*lmax**2,lmtop),
     $ bdbd(ntmx*npseudmx*lmax**2,nbfmax)
      real*8 rit(3)
      complex*16 hpvb(lmtop,nbfcmx,nchnl)
      complex*16 vecsch(npseudmx*ntmx*lmax**2,lmtop,nchnl)
      complex*16 vccsch(npseudmx*ntmx*lmax**2,lmtop)
      complex*16 cterm,cmp(nblok),ppmat(nblok,lmtop),zdotu,dterm
      complex*16 cj(nsblok,0:ltop),hpp,hvec(nblok,0:ltop)
      complex*16 hs(nsblok,0:ltop),ppvec(nqmx),pmvec(nqmx)
      complex*16 hpvhp(lmtop,lmtop,nchnl),hpvhm(lmtop,lmtop,nchnl)
      external fun
      dimension a(3,mxcen),icenps(npseudmx),lps(npseudmx),
     $exps(ntmx,lmax+1,npseudmx),cf(ntmx,lmax+1,npseudmx),
     $nn(ntmx,lmax+1,npseudmx),nt(npseudmx,lmax+1)
     $,expsf(npseudmx,lmax+1),nnf(npseudmx,lmax+1),h(npseudmx,lmax+1)
     $,ntf(npseudmx,lmax+1),fac(0:6),v(ntmx,ntmx),
     $vi(ntmx,ntmx),vsch((ntmx*(ntmx+1))/2,lmax*npseudmx),
     $basis(nblok,lmax*ntmx),lskip(npseudmx),
     $lcum(npseudmx)
      common/intg/alpa,beta,aa,bb,mm
      common/parm/nbig,mpt,lbig,mumax
      common/parms/mbig,mmpt,nbf
      equivalence (cj(1,1),cjr(1,1)),(hs(1,1),hsr(1,1))
      equivalence (ppvec(1),ppvr(1)),(pmvec(1),pmvr(1))
      equivalence (ppmat(1,1),ppmatr(1,1))
     $ ,(vecsch(1,1,1),vecschr(1,1,1)),(vccsch(1,1),vccschr(1,1))
      parameter (fourpi=4.*3.14159265358979)
      data fac/1.,1.,2.,6.,24.,120.,720./
      data xl/-0.9739065285,-0.8650633667,-0.6794095683,-0.4333953941,
     $-0.1488743390, 0.1488743390, 0.4333953941, 0.6794095683,
     $ 0.8650633667, 0.9739065285/
      data wl/ 0.0666713443, 0.1494513492, 0.2190863625, 0.2692667193,
     $ 0.2955242247, 0.2955242247, 0.2692667193, 0.2190863625,
     $ 0.1494513492, 0.0666713443/
      data tol/1.e-8/
      bd=0.
      inf=1
      eabs=1.e-8
      erel=1.e-5
      pi=3.14159265358
      sqpi=sqrt(pi)
      open(5,file='inpsnew')
      open(6,file='outpsnew')
      open(7,file='intsffps',form="unformatted")
      open(8,file='intsbfps',form="unformatted")
      open(14,file='bndmat',form='unformatted',status='old')
c      call openabs(9,'grid')
c
c read in centers
c     
      read(5,*)ncent,npt
      if (ncent .gt. mxcen) then
         write(6,*) 'Error ncent > mxcen'
         stop 'bas ncent'
      end if
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
c      call openabs(10,'ylms')
c      call openabs(13,'vbas')
      open(13,file='vbas',form='unformatted',status='old')
      open(10,file='ylms',form='unformatted',status='old')
c      call openabs(11,'scr11')
c      call openabs(20,'scr20')
      open(11,file='scr11',form='unformatted',status='unknown')
      open(20,file='scr20',form='unformatted',status='unknown')
      open(88,file='bessplr',form='unformatted',status='old')
c
c read parameters from ylm file and check
c
      read(10)nbig,mpt,lbig,mumax
      if(nbig.ne.ngrid.or.mpt.ne.npt)then
       write(6,701)nbig,ngrid,npt,mpt
701    format("stopping because of mismatch between ylm and grid files"
     1/" nbig,ngrid,npt,mpt :",4i5)
      stop
      endif
c
c read parameters from basis file and check
c     
      read(13)mbig,mmpt,nbf
      if(mbig.ne.ngrid.or.mmpt.ne.npt)then
         write(6,801)mbig,ngrid,npt,mmpt
 801  format("stopping because of mismatch between basis and grid files"
     1/" mbig,ngrid,npt,mmpt :",4i5)
         stop
         endif
      write(6,445)lbig,mumax
 445  format(/' info from ylm tape: lbig=',i3,'  mumax=',i3/)
      write(6,446)nbf
 446  format(/' info from vbas tape: nbf=',i3)
      do 4 i=1,ncent
 4    read(5,*)(a(iii,i),iii=1,3)
 100  format(3f10.5)
      write(6,*)' atomic centers'
      write(6,100)((a(iii,i),iii=1,3),i=1,ncent)
c     
c     which center is origin?
c     
      read(5,*)iorig
c     
c     read in pseudopotential parameters
c
      read(5,*)npseudo
      write(6,101)npseudo
 101  format(' number of effective core potentials =',i3)
      kl=0
      do 1 i=1,npseudo
         read(5,*)icenps(i),lps(i)
         write(6,102)i,icenps(i),lps(i)
 102  format(/'  ECP number ',i3,/' associated with center',i3,
     $ '   highest L-component is',i3)
         nl=lps(i)+1
         do 2 j=1,nl
          if(j.eq.1)then  
             write(6,103)lps(i)
          else 
             jj=j-2
             write(6,104)jj,lps(i)
          endif
 103      format(/' L=',i3,' component')
 104      format(/i3,'-',i3,' component')
          write(6,105)
 105      format('       exponent     coefficient   power')
          read(5,*)nt(i,j)  
          do 3 k=1,nt(i,j)
            read(5,*)exps(k,j,i),cf(k,j,i),nn(k,j,i)
            write(6,106)exps(k,j,i),cf(k,j,i),nn(k,j,i)
 106        format(2f15.6,i5)
 3        continue
 2       continue
      if(nl.gt.1)then
       write(6,107)
 107   format(/' separable representation of l-dependent terms')
       do 5 j=2,nl
       kl=kl+1
       jj=j-2
       write(6,104)jj,lps(i)
       write(6,109)
 109   format('       exponent     power    center')
       read(5,*)ntf(i,j),range,nnf(i,j)
       h(i,j)=range/(ntf(i,j)-1.)
       expsf(i,j)=1./h(i,j)/h(i,j)
       do 6 k=1,ntf(i,j)
            centr=(k-1.)*h(i,j)
            write(6,108)expsf(i,j),nnf(i,j),centr
 108        format(f15.6,i5,f15.6)
 6       continue
       nfit=ntf(i,j)
       n=nt(i,j)
       beta=expsf(i,j)
       do 7 ii=1,nfit
          aa=(ii-1.)*h(i,j)
          do 7 jj=1,ii
             bb=(jj-1.)*h(i,j)
            term=0.
            do 8 k=1,n
               mm=nn(k,j,i)+2*nnf(i,j)
               alpa=exps(k,j,i)
               call qagi(fun,bd,inf,eabs,erel,tab,aberr,neval,
     $             ier,limit,len,last,iwk,wk)
 8             term=term+cf(k,j,i)*tab
         v(jj,ii)=term   
         vi(ii,jj)=0.
         vi(jj,ii)=0.
 7    v(ii,jj)=v(jj,ii)
      do 77 ii=1,nfit
 77   vi(ii,ii)=1.
c      write(6,*)' potential matrix'
c      do 856 ii=1,nfit
c 856  write(6,303)(v(ii,jj),jj=1,nfit)
 303  format(10e12.4)
      call matinvr(v,nfit,vi,nfit,det,ntmx)
      ij=0
c      write(6,*)'potential inverse',kl
      do 9 ii=1,nfit
         do 10 jj=1,ii
            ij=ij+1
 10         vsch(ij,kl)=vi(ii,jj)
c            write(6,110)(vsch(jj,kl),jj=1,ii)
 110        format(10e12.4)
 9       continue
 5    continue
      endif
 1    continue
c     
c     read in free-function info
c
c read in spline information for free functions
c
       read(88)llmax,nr,rmin,rdel,alpha
c complex io equivalenced to real io
      nrr=2*nr
       read(88)(x(i),i=1,nr)
       read(88)((hsr(i,k),i=1,nrr),k=0,llmax)
c       read(88)((hsder(i,k),i=1,nr),k=0,llmax)
       read(88)
       read(88)((cjr(i,k),i=1,nrr),k=0,llmax)
c     read(88)((cy(i,k),i=1,nr),k=0,llmax)
       read(88)
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
      read(5,*)nchan
      do 22 ic=1,nchan
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
c     read basis function assignments
c
      read(5,*) ngauss(ic)
      ngic=ngauss(ic)
      if (ngic .gt. nbfcmx) then
         write(6,*) 'Error ngic larger than nbfcmx', ngic, nbfcmx
         stop 'stopping with bad ngic in psff'
      end if
      read(5,*) (ngch(j,ic),j=1,ngic)
      write(6,753) ic
  753 format(' basis functions for channel:',i3)
      write(6,'(5x,10i4)') (ngch(j,ic),j=1,ngic)
22    continue

      read(14) mchan
	if(nchan.ne.mchan)stop
      read(14) (echan(i),i=1,nchan)
      nchan2=nchan
      ignd=ismin(nchan,echan,1)
      read(14) nener
      read(14) (energy(i),i=1,nener)
      write(6,202)nchan,(echan(i),i=1,nchan)
202   format(' target energies for ',i3,' channels:',/,(2x,5e15.8))
      write(6,203) nener,(energy(i),i=1,nener)
 203  format(1x,i3,' incident energies:',/,(2x,5e15.8))
c     
c     load origin-centered ECP's, if present, into an array
c     otherwise, spline the ECP's for further use.
      ibingo=0
      jjread=0
      do 16 ips=1,npseudo
         if(icenps(ips).ne.iorig)go to 29
         ibingo=ips
         nts=nt(ips,1)
         go to 18
 29   continue
      jjread=jjread+1
      if(jjread.eq.1)then
         read(5,*)rsplmx,nspline
         rdspl=rsplmx/dfloat(nspline-1)
         rspmin=.00001
         rspl(1)=rspmin
         do 19 i=2,nspline
 19      rspl(i)=rdspl+rspl(i-1)
         rdsp26=(rdspl**2)/6.
         write(6,*)' spline info for non-origin centered ecps:'
         write(6,40)nspline,rsplmx
 40      format(' no. of points=',i4,'   max. value of r=',f10.2)
      endif
      nll=lps(ips)+1
      do 52 il=1,nll
      do 30 i=1,nspline
 30      uspl(i,ips,il)=0.
      nspl=nt(ips,il)
      do 31 i=1,nspline
         do 31 j=1,nspl
            arg=exps(j,il,ips)*rspl(i)**2
 31      uspl(i,ips,il)=uspl(i,ips,il)+cf(j,il,ips)*exp(-arg)*rspl(i)
     $    **nn(j,il,ips)
         yp1=(uspl(2,ips,il)-uspl(1,ips,il))/rdspl
         yp2=(uspl(nspline,ips,il)-uspl(nspline-1,ips,il))/rdspl
         call spline(rspl,uspl(1,ips,il),nspline,yp1,yp2,scr,
     $      cspl(1,ips,il))
c         write(6,*)'  rspl     uspl'
c         do 44 i=1,nspline
c 44      write(6,43)rspl(i),uspl(i,ips)
 52   continue
      go to 16
 18   continue
c     
c     read in quadrature info
c     
       write(6,*)' quadrature info for origin centered ecps:'
      read(5,*)rqmx,ninter
      nquad=ninter*10
      write(6,41)rqmx,ninter,nquad
 41   format(' the distance from the origin to',f6.2,
     $     ' will be split into',i3,' intervals with 10 pts per '
     $ ,'interval for a total of ',i4,'points')
      if(nquad.gt.nqmx)then
         write(6,*)' stopping because nquad exceeds nqmx'
         stop
      endif
      sint=rqmx/dfloat(ninter)
      aa=0.
      b=sint
      j=0
      do 38 i=1,ninter
         amb=(b-aa)/2.
         apb=(b+aa)/2.
         do 39 k=1,10
            j=j+1
            rr(j)=amb*xl(k)+apb
 39      wt(j)=amb*wl(k)
         aa=aa+sint
 38   b=b+sint
      nl=lps(ibingo)+1
      do 20 j=1,nl
      do 20 i=1,nquad
 20   u(i,j)=0.
      do 21 k=1,nl
         nts=nt(ibingo,k)
      do 21 i=1,nquad
      do 21 j=1,nts
         arg=exps(j,k,ibingo)*rr(i)**2
         u(i,k)=u(i,k)+cf(j,k,ibingo)*wt(i)*exp(-arg)
     $      *rr(i)**nn(j,k,ibingo)
 21   continue
 16   continue
c      write(6,*)'  rr       wt          wt*u'
c      do 42 i=1,nquad
c 42   write(6,43)rr(i),wt(i),u(i)
 43   format(f9.5,2e12.4)
c     
c make sure there are l-dependent, off-center terms present     
c
      nldep=0
      ii=0
      lcum(1)=0
      do 48 i=1,npseudo
         nclm(0,i)=0
         do 449 j=1,lps(i)
 449     nclm(j,i)=nclm(j-1,i)+ntf(i,j+1)*(2*j-1)
         if(i.eq.ibingo)go to 48
         ii=ii+1
         lskip(ii)=lps(i)
         nldep=nldep+lps(i)
 48   continue
      npsoff=ii
      if(nldep.eq.0)go to 49
      do 448 i=2,npsoff
 448  lcum(i)=lcum(i-1)+lskip(i-1)
c     
c     generate off-center ylm's and basis functions 
c     on grid 
c
      iset=1
      iread=0
      iquit=0
      jset=1
      jwhere=1
c
c read in a block of grid points and transfer to a temporary location
c
      marg=min0(ngrid,npt)
      nread=4*marg
      call rdabs(9,buff(1),nread,iset)
c************tnr
c      iset=iset+nread
      iset=iset+1
      iread=iread+marg
 46   continue  
      narg=marg
      call dcopy(narg,buff(1),4,xx(1),1)
      call dcopy(narg,buff(2),4,yy(1),1)
      call dcopy(narg,buff(3),4,zz(1),1)
      call dcopy(narg,buff(4),4,ww(1),1)
      iremn=ngrid-iread
      if(iremn.eq.0)then
        iquit=1
        go to 47
      endif
      marg=min0(iremn,npt)
      nread=4*marg
      call rdabs(9,buff(1),nread,iset)
c************tnr
c      iset=iset+nread
      iset=iset+1
      iread=iread+marg
 47   continue
      do 50 ips=1,npseudo
         jbasis(ips)=0
         jj=0
         if(ips.eq.ibingo.or.lps(ips).eq.0)go to 50
         ldep=1+lps(ips)
c     
c   basis functions
c     
         do 51 i=1,narg
         xarg=(xx(i)-a(1,icenps(ips)))   
         yarg=(yy(i)-a(2,icenps(ips)))
         zarg=(zz(i)-a(3,icenps(ips)))
         arg2=(xarg**2+yarg**2 +zarg**2) 
         rvec(i)=sqrt(arg2)
         xplm(i)=zarg/rvec(i)
         xsq=xplm(i)**2
         xqq=abs(1.-xsq)
         if (xqq.lt.tol) then
            cphi(i)=1.
            sphi(i)=0.
            if(xplm(i)) 4466,4467,4468
 4466          xplm(i)=-1.
             go to 51
 4467          write(6,*) 'error '
             stop
 4468          xplm(i)=1.
             else
            cphi(i)=xarg/rvec(i)/sqrt(1.-xplm(i)**2)
         sphi(i)=yarg/rvec(i)/sqrt(1.-xplm(i)**2)
         endif
 51   continue

         do 53 il=2,ldep
            nbas=ntf(ips,il)
            do 53 j=1,nbas
            aa=h(ips,il)*(j-1.)
            jj=jj+1
            do 53 i=1,narg
               if(rvec(i).gt.rsplmx)then
                  basis(i,jj)=0.
               else
                  arg=expsf(ips,il)*(rvec(i)-aa)**2
                  basis(i,jj)=exp(-arg)*rvec(i)**nnf(ips,il)
               endif
 53         continue
c
         jwrite=nblok*jj
c         call wrabs(20,basis,jwrite,jset)
         write(20)((basis(i,j),i=1,nblok),j=1,jj)
c         jset=jset+jwrite
c         jset=jset+1
         jbasis(ips)=jj
c
c ylm's
c     
       mubig=lps(ips)-1           
       do 54 mu=0,mubig
          call plm(xplm,narg,mu,mubig,p
     $      ,pscr(1+nblok),pscr(1+2*nblok),pscr(1+5*nblok)
     $      ,pscr,pscr(1+3*nblok),pscr(1+4*nblok),pscr(1+6*nblok)
     $      ,nblok)
         if(mu.eq.0)then
            do 55 i=0,mubig
               const=sqrt((2*i+1)/fourpi)
               do 555 j=1,narg
                  p(j,i)=p(j,i)*const
 555           continue
 55            continue
               do 56 i=0,mubig
                  call scopy(narg,p(1,i),1,biff(1+i*narg),1)
 56            continue
               jbuf=narg*(mubig+1)
c               call wrabs(11,biff(1),jbuf,jwhere)
               write(11) (biff(i),i=1,jbuf)
c               jwhere=jwhere+jbuf
               jwhere=jwhere+1
            else
               do 57 i=1,narg
                  cmp(i)=(cmplx(cphi(i),sphi(i)))**mu
                  cmphi(i)=dble(cmp(i))
                  smphi(i)=imag(cmp(i))
 57            continue
               do 58 i=mu,mubig
                  const=sqrt((2*i+1)/fourpi*fac(i-mu)/fac(i+mu))
c     sqrt(2) factor added to normalize "real valued" ylms
                  const=const*sqrt(2.0)
                  do 588 j=1,narg
                     pp(j,i)=p(j,i)*const*cmphi(j)
                     pm(j,i)=p(j,i)*const*smphi(j)
 588              continue
 58         continue
                  do 59 i=mu,mubig
                     call scopy(narg,pp(1,i),1,biff(1+(i-mu)*2*narg),1)
                     call scopy(narg,pm(1,i),1,biff(1+narg+(i-mu)*
     $                  2*narg),1)
 59               continue
                  jbuf=2*narg*(1+mubig-mu)
c                  call wrabs(11,biff(1),jbuf,jwhere)
               write(11) (biff(i),i=1,jbuf)
c                  jwhere=jwhere+jbuf
                  jwhere=jwhere+1
               endif
 54    continue
 50   continue
c
c return to start a fetch another block of points
c
      if(iquit.eq.0)go to 46
 49   continue
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
      iflag=1
      write(7) iflag
      write(7) nener,nchan,(nlm(ic),ic=1,nchan)
      write(7) ((lch(j,ic),mch(j,ic),j=1,nlm(ic)),ic=1,nchan)
      write(7) eground
c
c open a loop on energies
c
      call second(time)
      izero=ntmx*npseudmx*lmax**2
      do 83 i=1,izero
         do 83 j=1,nbfmax
 83      bdbd(i,j)=0.
      do 1000 iene=1,nener
c
c reset all energy-independent file pointers
c
         rewind 20
         rewind 11
      iset=1
      iread=0
      iquit=0
c      jwhere=4
      iwrite=0
c      kset=3
         jset=1
         kwhere=1
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
 17   continue
      write(6,717) energy(iene), (kchan(i),i=1,nchan)
 717  format(//,' incident E = ',f10.6,/,' channel momenta = ',(6e12.5))

c initialize free-free matrix elements
c
      nzero=npseudo*ntmx*(lmax)**2
      do 11 i=1,lmtop
      do 11 ist=1,nchan
         do 111 j=1,nzero
 111     vecsch(j,i,ist)=0.
         do 112 k=1,nbfcmx
 112     hpvb(i,k,ist)=0.   
      do 11 j=1,lmtop
      hpvhm(j,i,ist) = 0.0
11    hpvhp(j,i,ist) = 0.0
c     
c     compute contribution to free-free integrals from L-independent terms
      if(ibingo.eq.0)go to 12
      do 14 ic=1,nchan
         if(iclosed(ic).eq.1)go to 14
c     
c     do origin-centered contribution first
c
      nlmic=nlm(ic)
      lprior=-1
            do 15 ilm=1,nlmic
               l=lch(ilm,ic)
               if(l.eq.lprior)go to 244
c     
c L-independent part
c
            do 24 i=1,nquad
               arg=rr(i)*kchan(ic)
               klo=(arg-rmin)/rdel
               klo=klo+1
               aa=(x(klo+1)-arg)/rdel
               b=(arg-x(klo))/rdel
               hpp=aa*hs(klo,l)+b*hs(klo+1,l)+(aa*(aa*aa-1.)*
     $            cj(klo,l)+b*(b*b-1.)*cj(klo+1,l))*rd26
c               hpp=sin(arg)/arg
               ppvec(i)=hpp*hpp
 24            pmvec(i)=hpp*imag(hpp)
            ppr=ddot(nquad,u,1,ppvr(1),2)
            ppi=ddot(nquad,u,1,ppvr(2),2)
            pmr=ddot(nquad,u,1,pmvr(1),2)
            pmi=ddot(nquad,u,1,pmvr(2),2)
 244        hpvhp(ilm,ilm,ic)=hpvhp(ilm,ilm,ic)+cmplx(ppr,ppi)
     $        *kchan(ic)
            hpvhm(ilm,ilm,ic)=hpvhm(ilm,ilm,ic)+cmplx(pmr,pmi)   
     $        *kchan(ic)
c     
c L-dependent parts
c
            if(nl.eq.1)go to 15
            if(l.eq.lprior)go to 246
            qqr=0.
            qqi=0.
            qmr=0.
            qmi=0.
            do 245 i=2,nl
               ll=i-2
               if(ll.eq.l)then
            qqr=ddot(nquad,u(1,i),1,ppvr(1),2)
            qqi=ddot(nquad,u(1,i),1,ppvr(2),2)
            qmr=ddot(nquad,u(1,i),1,pmvr(1),2)
            qmi=ddot(nquad,u(1,i),1,pmvr(2),2)
            go to 246
            endif
 245     continue
 246     hpvhp(ilm,ilm,ic)=hpvhp(ilm,ilm,ic)+cmplx(qqr,qqi)
     $        *kchan(ic)
         hpvhm(ilm,ilm,ic)=hpvhm(ilm,ilm,ic)+cmplx(qmr,qmi)   
     $           *kchan(ic)
 15         continue
            write(6,*)' hpvhp'
            call cmwrite(hpvhp(1,1,ic),nlmic,nlmic,lmtop)
            write(6,*)' hpvhm'
            call cmwrite(hpvhm(1,1,ic),nlmic,nlmic,lmtop)
 14      continue
 12      continue
c     
c     do contributions from other centers
c     
      if(ibingo.ne.0.and.npseudo.eq.1)go to 1000

c
c read in a block of grid points and transfer to a temporary location
c
      marg=min0(ngrid,npt)
      nread=4*marg
      call rdabs(9,buff(1),nread,iset)
c      iset=iset+nread
      iset=iset+1
      iread=iread+marg
 32    continue  
      iwrite=iwrite+1
      narg=marg
      call scopy(narg,buff(1),4,xx(1),1)
      call scopy(narg,buff(2),4,yy(1),1)
      call scopy(narg,buff(3),4,zz(1),1)
      call scopy(narg,buff(4),4,ww(1),1)
      iremn=ngrid-iread
      if(iremn.eq.0)then
        iquit=1
        go to 34
      endif
      marg=min0(iremn,npt)
      nread=4*marg
      call rdabs(9,buff(1),nread,iset)
c      iset=iset+nread
      iset=iset+1
      iread=iread+marg
34    continue
      do 28 i=1,narg
  28  rvec(i)=sqrt(xx(i)**2 + yy(i)**2 +zz(i)**2)
c
c read in a block of gaussians
c
c      call rdabs(13,gauss,nbf*narg,kset)
      read(13)(gauss(ib),ib=1,nbf*narg)
c      kset=kset+nbf*narg
c
c skip over basis function second derivatives
c
c      kset=kset+nbf*narg
      read(13)
c     
c     calculate LMAX(uu) and L-dependent(uuu) componentss of non-origen-centered ECP's at these points
c
      do 35 i=1,narg
         uu(i)=0.   
         do 335 j=1,nldep
 335     uuu(i,j)=0.
      jj=0   
      do 33 j=1,npseudo
         if(j.eq.ibingo)go to 33
         jj=jj+1
         arg2=((xx(i)-a(1,icenps(j)))**2 
     $            +(yy(i)-a(2,icenps(j)))**2 
     $            +(zz(i)-a(3,icenps(j)))**2) 
         arg=sqrt(arg2)
         if(arg.gt.rsplmx)go to 33
          klo=(arg-rspmin)/rdspl
          klo=klo+1
          aa=(rspl(klo+1)-arg)/rdspl
          b=(arg-rspl(klo))/rdspl
          hpr=aa*uspl(klo,j,1)+b*uspl(klo+1,j,1)+(aa*(aa*aa-1.)*
     $      cspl(klo,j,1)+b*(b*b-1.)*cspl(klo+1,j,1))*rdsp26
          uu(i)=uu(i)+hpr/arg2
          do 333 k=1,lps(j)
             hpr=aa*uspl(klo,j,k+1)+b*uspl(klo+1,j,k+1)+(aa*(aa*aa-1.)*
     $            cspl(klo,j,k+1)+b*(b*b-1.)*cspl(klo+1,j,k+1))*rdsp26
             kk=k+lcum(jj)
 333      uuu(i,kk)=hpr/arg2
 33   continue
 35   continue
c      if(iwrite.eq.iwww)then
c               write(6,*)iwww,narg,kchan(ic)
c               write(6,654)(rvec(i),uu(i),ww(i),i=1,narg)
 654           format(4e12.4)
c               endif
c
c read in a block of ylm's associated with free functions
c
      do 266 m=0,mumax
      if(m.eq.0)then
      jbuf=narg*(lbig+1)
      else
      jbuf=2*narg*(lbig-m+1)
      endif
c      call rdabs(10,biff(1),jbuf,jwhere)
      read(10)(biff(ib),ib=1,jbuf)
c      jwhere=jwhere+jbuf
      if(m.eq.0)then
      do 26 l=0,lbig
      do 26 i=1,narg
26    ylm(i,l,0)=biff(i+narg*l)
      else
      do 27 l=m,lbig
      do 27 i=1,narg
      ylm(i,l,2*m-1)=biff(i+(l-m)*2*narg)
      ylm(i,l,2*m) = biff(i+(l-m)*2*narg +narg)
27    continue
      endif
266   continue
c     
c loop over pseudopotentials
c
      kk=0
      jlmpr=0
      do 13 ips=1,npseudo
         if(ips.eq.ibingo)go to 13
         kk=kk+1
c*****************!!!!!!!*******************
c         jj=nblok*jbasis(ips)
         jj=jbasis(ips)
         ldep=1+lps(ips)
         icum(0)=0
         do 79 i=3,ldep
 79      icum(i-2)=icum(i-3)+ntf(ips,i-1)
         if(jj.gt.0)then
c     
c     get fitting functions
c     
c            call rdabs(20,basis,jj,jset)
         read(20)((basis(i,j),i=1,nblok),j=1,jj)
c            jset=jset+jj
            jset=jset+1
c
c read in a block of ylm's associated with ECP's
c
            mubig=lskip(kk)-1
            do 61 m=0,mubig
               if(m.eq.0)then
                  mtag=1
                  jbuf=narg*(mubig+1)
               else
                  mtag=2*m
                  jbuf=2*narg*(mubig-m+1)
               endif
c               call rdabs(11,barf(1),jbuf,kwhere)
               read(11) (barf(i),i=1,jbuf)
c               kwhere=kwhere+jbuf
               kwhere=kwhere+1
               jstart=0
               do 62 l=m,mubig
                  lm=l**2+mtag
                  do 663 i=1,narg
  663              vec(i,lm)=barf(i+jstart)*uuu(i,lcum(kk)+l+1)
                   if(iene.eq.1)then
                      nbas=ntf(ips,l+2)
                      do 763 ii=1,nbas
                         ibtag=ii+nbas*(mtag-1)+nclm(l,ips)+jlmpr
                         do 773 jj=1,narg
 773                     xx(jj)=basis(jj,ii+icum(l))*ww(jj)*vec(jj,lm)
                         do 763 j=1,nbf
                         bdbd(ibtag,j)=bdbd(ibtag,j)+
     $                           ddot(narg,xx,1,gauss(narg*(j-1)+1),1)
 763                  continue   
                      endif
                  jstart=jstart+narg
                  if(m.ne.0)then
                     do 64 i=1,narg
c 64                  vec(i,lm+1)=barf(i+jstart+narg)*uuu(i,lcum(kk)+l+1)
  64                  vec(i,lm+1)=barf(i+jstart)*uuu(i,lcum(kk)+l+1)
                     if(iene.eq.1)then
                        do 764 ii=1,nbas
                         ibtag=ii+nbas*(mtag)+nclm(l,ips)+jlmpr
                         do 774 jj=1,narg
 774                     xx(jj)=basis(jj,ii+icum(l))*ww(jj)*vec(jj,lm+1)
                         do 764 j=1,nbf
                         bdbd(ibtag,j)=bdbd(ibtag,j)+
     $                           ddot(narg,xx,1,gauss(narg*(j-1)+1),1)
 764                     continue
                     endif
                     jstart=jstart+narg
                  endif
 62            continue
 61         continue
            lm=0
            jlm=0
            jprev=0
            do 65 l=0,mubig
            mm=2*l+1
            nbas=ntf(ips,l+2)
            do 66 m=1,mm
               lm=lm+1
               do 67 j=1,nbas
                  jlm=jlm+1
                  jj=j+jprev
                  do 68 i=1,narg
 68                 bvec(i,jlm)=vec(i,lm)*basis(i,jj)
 67               continue
 66            continue
 65         continue
            njlm(ips)=jlm
      endif
c     
c     loop on channels
c     
         do 25 ic=1,nchan
         if(iclosed(ic).eq.1)go to 25
            nlmic=nlm(ic)
c     
c     calculate and store free functions at points
c     
c     do 121 i=1,narg
           do 36 i=1,narg
            arg=kchan(ic)*rvec(i)
               klo=(arg-rmin)/rdel
               klo=klo+1
               aa=(x(klo+1)-arg)/rdel
               b=(arg-x(klo))/rdel
            do 36 l=0,llmax
                hvec(i,l)=aa*hs(klo,l)+b*hs(klo+1,l)+(aa*(aa*aa-1.)*
     $          cj(klo,l)+b*(b*b-1.)*cj(klo+1,l))*rd26
c                hvec(i,l)=sin(arg)/arg
  36        continue
c     
c do l-independent terms
c     
c     
c     note that uu is already summed over all ECPs, so only do this part once
c
            if(ips.eq.1)then
               ngic=ngauss(ic)
            do 37 j=1,nlmic
               l=lch(j,ic)
               m=mch(j,ic)
               do 81 i=1,narg
 81            ppvec(i)=ww(i)*hvec(i,l)*uu(i)*ylm(i,l,m)
               do 82 i=1,ngic
                ppr=ddot(narg,ppvr(1),2,gauss(1+narg*(ngch(i,ic)-1)),1)
                ppi=ddot(narg,ppvr(2),2,gauss(1+narg*(ngch(i,ic)-1)),1)
 82          hpvb(j,i,ic)=hpvb(j,i,ic)+sqrt(kchan(ic))*cmplx(ppr,ppi)
             do 37 k=1,nlmic
               lp=lch(k,ic)
               mp=mch(k,ic)
               do 121 i=1,narg
               term=kchan(ic)*ww(i)*uu(i)*ylm(i,l,m)*ylm(i,lp,mp)
               cterm=hvec(i,l)*hvec(i,lp)
               rrvec(i)=term
               ppvec(i)=cterm
 121           pmvec(i)=hvec(i,l)*imag(hvec(i,lp))
               ppr=ddot(narg,rrvec,1,ppvr(1),2)
               ppi=ddot(narg,rrvec,1,ppvr(2),2)
               pmr=ddot(narg,rrvec,1,pmvr(1),2)
               pmi=ddot(narg,rrvec,1,pmvr(2),2)
               hpvhp(j,k,ic)=hpvhp(j,k,ic)+cmplx(ppr,ppi)
               hpvhm(j,k,ic)=hpvhm(j,k,ic)+cmplx(pmr,pmi)
 37         continue
            endif
c     
c do l-dependent terms
c     
            do 63 ilm=1,nlmic
               l=lch(ilm,ic)
               m=mch(ilm,ic)
               do 63 i=1,narg
 63            ppmat(i,ilm)=hvec(i,l)*ylm(i,l,m)*ww(i)
               do 69 ilm=1,nlmic
                  do 69 jlm=1,njlm(ips)
                     ppr=ddot(narg,ppmatr(1,ilm),2,bvec(1,jlm),1) 
                     ppi=ddot(narg,ppmatr(2,ilm),2,bvec(1,jlm),1) 
                     vecsch(jlm+jlmpr,ilm,ic)=
     $                     vecsch(jlm+jlmpr,ilm,ic)+cmplx(ppr,ppi)
 69               continue
 25      continue
         jlmpr=jlmpr+njlm(ips)
 13      continue
c
c return to start a fetch another block of points
c
      if(iquit.eq.0)go to 32
c     
c assemble l-dependent off-center terms
c     
c      write(6,*)'  ic  ilm  ips   l    m     ppr         ppi'
      do 70 ic=1,nchan
         if(iclosed(ic).eq.1)go to 70
         nlmic=nlm(ic)
c         write(6,*)' vecsch'
c         llll=npseudmx*ntmx*lmax**2
c         call cmwrite(vecsch(1,1,ic),jlmpr,nlmic,llll)
         do 71 ilm=1,nlmic
            jlmpr=0
            kk=0
            kl=0
            do 72 ips=1,npseudo
               if(ips.eq.ibingo)go to 72
               kk=kk+1
               mubig=lskip(kk)-1
               lm=0
               do 73 l=0,mubig
                  nbas=ntf(ips,l+2)
                  kl=kl+1
                  ii=0
                  do 74 i=1,nbas
                     do 74 j=1,i
                        ii=ii+1
                        vi(i,j)=vsch(ii,kl)
 74                  vi(j,i)=vi(i,j)
                     mm=2*l+1
                     do 75 m=1,mm
                        marker=2*(lm+jlmpr)
                        do 76 ibas=1,nbas
                           ppr=ddot(nbas,vi(1,ibas),1,
     $                       vecschr(1+marker,ilm,ic),2)
                           ppi=ddot(nbas,vi(1,ibas),1,
     $                       vecschr(2+marker,ilm,ic),2)
c                           write(6,776)ic,ilm,ips,l,m,ppr,ppi
 776                       format(5i5,2e12.4)
 76                     ppvec(ibas)=cmplx(ppr,ppi)
                        do 677 ibas=1,nbas
 677                    vccsch(ibas+lm+jlmpr,ilm)=ppvec(ibas)
 75                  lm=lm+nbas
 73               continue
 72            jlmpr=jlmpr+njlm(ips)
 71         continue
c     
c     complete free-free terms
c     
            write(6,*)' hpvb before l-dependent part'
            ngic=ngauss(ic)
            call cmwrite(hpvb(1,1,ic),nlmic,ngic,lmtop)
            write(6,*)' hpvhp before l-dependent part'
            call cmwrite(hpvhp(1,1,ic),nlmic,nlmic,lmtop)
            write(6,*)' hpvhm before l-dependent part'
            call cmwrite(hpvhm(1,1,ic),nlmic,nlmic,lmtop)
            do 78 i=1,nlmic
               do 778 j=1,nlmic
                  cterm=zdotu(jlmpr,vccsch(1,i),1,vecsch(1,j,ic),1)
                  hpvhp(i,j,ic)=cterm*kchan(ic)+hpvhp(i,j,ic)
                  ppr=ddot(jlmpr,vccschr(1,i),2,vecschr(2,j,ic),2)
                  ppi=ddot(jlmpr,vccschr(2,i),2,vecschr(2,j,ic),2)
                  hpvhm(i,j,ic)=cmplx(ppr,ppi)*kchan(ic)+hpvhm(i,j,ic)
 778            continue
c     
c     complete bound-free terms
c
c                write(6,*)'bdbd'
              ngic=ngauss(ic)
               do 80 j=1,ngic
                  ilabel=ngch(j,ic)
c                  write(6,888)ilabel,(bdbd(ii,ilabel),ii=1,jlmpr)
 888              format(i5,10e12.4/(5x,10e12.4))
                  ppr=ddot(jlmpr,vccschr(1,i),2,bdbd(1,ilabel),1)
                  ppi=ddot(jlmpr,vccschr(2,i),2,bdbd(1,ilabel),1)
 80            hpvb(i,j,ic)=hpvb(i,j,ic)
     $            +sqrt(kchan(ic))*cmplx(ppr,ppi)
 78         continue
70      continue
         do 45 ic=1,nchan
            nlmic=nlm(ic)
            write(6,*)' hpvhp'
            call cmwrite(hpvhp(1,1,ic),nlmic,nlmic,lmtop)
            write(6,*)' hpvhm'
            call cmwrite(hpvhm(1,1,ic),nlmic,nlmic,lmtop)
            write(6,*)' hpvb'
            ngic=ngauss(ic)
            call cmwrite(hpvb(1,1,ic),nlmic,ngic,lmtop)
 45      continue
c     
c     write results to disk
c     
         write(7) (kchan(ic),ic=1,nchan)
         do 84 ic=1,nchan
            nlmic=nlm(ic)
 84      write(7) ((hpvhp(ilm,jlm,ic),ilm=1,nlmic),jlm=1,nlmic)
         do 85 ic=1,nchan
            nlmic=nlm(ic)
 85      write(7) ((hpvhm(ilm,jlm,ic),ilm=1,nlmic),jlm=1,nlmic)
         write(8) (kchan(ic),ic=1,nchan)
         do 86 ic=1,nchan
            nlmic=nlm(ic)
            ngic=ngauss(ic)
            write(8)((hpvb(ilm,ig,ic),ilm=1,nlmic),ig=1,ngic)
 86      continue
      rewind 13
      read(13)
      rewind 10
      read(10)
 1000 continue
      call second(tnow)
      tnow=tnow-time
      write(6,*)' time for run is ',tnow
      stop
      end
