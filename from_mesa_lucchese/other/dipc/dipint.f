      subroutine dipint(hpvb,hpvbx,hpvby,hpvbz,nchan2,nlm,lch,mch
     1  ,nchnl,lmtop,ngauss,nbfmax,cvec,dvec,basis,grid,wt,ylm,
     2  hp,nblok,narg,ltop,hbx,hby,hbz,iene)
c
c  bound-free matrices
c
      implicit real*8 (a-h,o-z)
      real*8 ylm(nblok,0:ltop,0:2*ltop)
      real*8 basis(nbfmax*nblok)
      dimension nlm(nchnl),lch(lmtop,nchnl),mch(lmtop,nchnl)
      dimension grid(nblok,3)
      real*8 cvec(nblok),dvec(nblok)
        dimension hbx(nbfmax,nbfmax)
        dimension hby(nbfmax,nbfmax)
        dimension hbz(nbfmax,nbfmax)
        complex*16 hpvb(lmtop,nbfmax,nchnl)
        complex*16 hpvbx(lmtop,nbfmax,nchnl)
        complex*16 hpvby(lmtop,nbfmax,nchnl)
        complex*16 hpvbz(lmtop,nbfmax,nchnl)
      real*8 wt(nblok)
      complex*16 hp(nblok,0:ltop,nchnl),ai
      ai=(0.,1.)
c
c compute dipole matrix
c
      do 500 ic=1,nchan2
      nlmic=nlm(ic)
      do 500 ilm=1,nlmic
      l1=lch(ilm,ic)
      m1=mch(ilm,ic)
      ngjc = ngauss
      do 400 ig=1,ngjc
      do 300 i=1,narg
      term=ylm(i,l1,m1)*basis(narg*(ig-1)+i)*wt(i)
      cvec(i)=term*dble(hp(i,l1,ic))
      dvec(i)=term*imag(hp(i,l1,ic))
300   continue
      hpvbx(ilm,ig,ic)=ddot(narg,cvec,1,grid(1,1),1)
     1 + hpvbx(ilm,ig,ic)+ai*ddot(narg,dvec,1,grid(1,1),1)
      hpvby(ilm,ig,ic)=ddot(narg,cvec,1,grid(1,2),1)
     1 + hpvby(ilm,ig,ic)+ai*ddot(narg,dvec,1,grid(1,2),1)
      hpvbz(ilm,ig,ic)=ddot(narg,cvec,1,grid(1,3),1)
     1 + hpvbz(ilm,ig,ic)+ai*ddot(narg,dvec,1,grid(1,3),1)
        do 405 i=1,narg
      hpvb(ilm,ig,ic)=cvec(i)+ai*dvec(i) + hpvb(ilm,ig,ic)
405     continue
400   continue
500   continue
c
       if(iene.eq.1)then
	do 600 ig=1,ngauss
	do 600 jg=1,ngauss
	do 601 i=1,narg
	fac= wt(i)*basis(narg*(ig-1)+i)*basis(narg*(jg-1)+i)
	hbx(ig,jg)=hbx(ig,jg)+fac*grid(i,1)
	hby(ig,jg)=hby(ig,jg)+fac*grid(i,2)
	hbz(ig,jg)=hbz(ig,jg)+fac*grid(i,3)
c	hbz(ig,jg)=hbz(ig,jg)+fac
601	continue
600	continue
      endif
c
      return
      end
