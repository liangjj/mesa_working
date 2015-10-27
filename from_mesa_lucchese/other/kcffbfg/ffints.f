      subroutine ffints(nchan2,iflag,narg,iclosed)
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
	dimension iclosed(1)
      complex*16 cvec,dvec,evec,hp,ovbf,hpvb,hpvhp,hpvhm,ai,toboth
      real*8 crvec(2*nblok),rvec(nblok)
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
      complex*16 zdotu
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
	if(iclosed(ic).eq.1.or.iclosed(jc).eq.1)then
         hpvhp(ilm,jlm,ist) = 0.
         hpvhm(ilm,jlm,icjc) = 0.
         hpvhm(jlm,ilm,jcic) = 0.
         go to 400
	else
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
      dvec(i)=fac*potval*imag(hp(i,l2,jc))
      evec(i)=fac*potval*imag(hp(i,l1,ic))
299   continue
      endif
      hpvhp(ilm,jlm,ist) = hpvhp(ilm,jlm,ist) +
     1 zdotu(narg,cvec,1,hp(1,l2,jc),1)
      hpvhm(ilm,jlm,icjc) = hpvhm(ilm,jlm,icjc) +
     1 zdotu(narg,hp(1,l1,ic),1,dvec,1)
      hpvhm(jlm,ilm,jcic) = hpvhm(jlm,ilm,jcic) +
     1 zdotu(narg,hp(1,l2,jc),1,evec,1)
	endif
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
	if(iclosed(ic).eq.1)then
         hpvhp(ilm,ilm,icc) = 0.
         hpvhm(ilm,ilm,iic) = 0.
         go to 501
	else
      if(iflag.eq.0) then
      do 301 i=1,narg
      cvec(i)=wt(i)*ylm(i,l1,m1)**2*hp(i,l1,ic)
      rvec(i) = hd(i,l1,ic)
301   continue
      toboth= ddot(narg,rvec,1,crvec,1)+
     2 ddot(narg,rvec,1,crvec(2),2)*ai
      hpvhp(ilm,ilm,icc)=hpvhp(ilm,ilm,icc)+
     1    toboth
      hpvhm(ilm,ilm,iic)=hpvhm(ilm,ilm,iic)+
     1   toboth
      else
      do 302 i=1,narg
      cvec(i)=wt(i)*ylm(i,l1,m1)**2*hp(i,l1,ic)
302   continue
      hpvhp(ilm,ilm,icc)=hpvhp(ilm,ilm,icc)+
     1     ddot(narg,crvec,2,hd(1,l1,ic),1)+
     2 ddot(narg,crvec(2),2,hd(1,l1,ic),1)*ai
      endif
      endif
501   continue
      return
      end
