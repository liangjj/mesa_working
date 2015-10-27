      subroutine bfints(hpvb,ovbf,nchan2,nlm,lch,mch,nchnl,lmtop
     1  ,ngauss, ngch,nbfmax,cvec,basis,vbuf,wt,ylm,hp,hd,nblok,
     2  nstate,narg,ltop,iclosed)
c
c  bound-free matrices
c
      implicit real*8 (a-h, o-z)
      real*8 ylm(nblok,0:ltop,0:2*ltop)
      real*8 vbuf(nblok*nchnl*(nchnl+1)/2)
      real*8 basis(nbfmax*nblok)
      dimension nlm(nchnl),lch(lmtop,nchnl),mch(lmtop,nchnl)
      dimension ngauss(nchnl),ngch(nbfmax,nchnl)
      dimension iclosed(1)
      complex*16 zdotu,cvec(nblok),hpvb(lmtop,nbfmax,nstate)
      complex*16 ovbf(lmtop,nbfmax,nchnl)
      real*8 wt(nblok)
      complex*16 hp(nblok,0:ltop,nchnl), hd(nblok,0:ltop,nchnl)
c
c compute potential matrix
c
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
c     
c     closed channels
c     
      if(iclosed(ic).eq.1) then
      hpvb(ilm,ig,ist)=0.0
      go to 400
      else
c
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
c     
c     closed channels
c     
      if(iclosed(ic).eq.1) then
      hpvb(ilm,ig,icc)=0.0
      ovbf(ilm,ig,ic)=0.0
      go to 401
      else
c
      do 301 i=1,narg
      cvec(i)=wt(i)*ylm(i,l1,m1)*basis(narg*(nbch-1)+i)
301   continue
      ovbf(ilm,ig,ic)=zdotu(narg,hp(1,l1,ic),1,cvec,1)
     1 + ovbf(ilm,ig,ic)
      hpvb(ilm,ig,icc)=hpvb(ilm,ig,icc)+zdotu(narg,hd(1,l1,ic),1,cvec
     1 ,1)
      endif
401   continue
501   continue
      return
      end
