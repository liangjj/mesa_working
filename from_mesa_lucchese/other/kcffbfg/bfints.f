      subroutine bfints(nchan2,narg,iclosed)
      implicit real*8 (a-h,o-z)
      parameter (nblok=500,mtop=10,ltop=10,nchnl=10)
      parameter (nsblok=3000)
c     parameter (lmtop=(2*mtop+1)*(ltop+1)-mtop*(mtop+1))
      parameter (lmtop=30)
c      parameter (lmtop=10)
      parameter (nstate=nchnl**2)
      parameter (nsymst = nchnl*(nchnl+1)/2)
      parameter (nbfmax=500)
c      parameter (nbfcmx =100)
      parameter (nbfcmx =80)
      parameter (maxene=200)
      complex*16 cvec,dvec,evec,hp,ovbf,hpvb,hpvhp,hpvhm,ai
      real*8 crvec(2*nblok)
	dimension iclosed(1)
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
      complex*16 zdotu
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
	if(iclosed(ic).eq.1)then
         hpvb(ilm,ig,ist)=0.
         go to 400
        else
      do 300 i=1,narg
      cvec(i)=ylm(i,l1,m1)*basis(narg*(nbch-1)+i)*vbuf(i+ipot)
300   continue
      hpvb(ilm,ig,ist)=zdotu(narg,cvec,1,hp(1,l1,ic),1)
     1 + hpvb(ilm,ig,ist)
	endif
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
	if(iclosed(ic).eq.1)then
	hpvb(ilm,ig,icc)=0.
	ovbf(ilm,ig,ic)=0.
	go to 401
        else
      do 301 i=1,narg
      cvec(i)=wt(i)*ylm(i,l1,m1)*basis(narg*(nbch-1)+i)
301   continue
      ovbf(ilm,ig,ic)=zdotu(narg,hp(1,l1,ic),1,cvec,1)
     1 + ovbf(ilm,ig,ic)
      hpvb(ilm,ig,icc)=hpvb(ilm,ig,icc)+ddot(narg,hd(1,l1,ic),1,crvec
     1 ,2)+ai*ddot(narg,hd(1,l1,ic),1,crvec(2),2)
	endif
401   continue
501   continue
      return
      end
