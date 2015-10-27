      subroutine ffints(hpvhp,hpvhm,nchan2,nlm,lch,mch,nchnl,lmtop
     1  ,cvec,dvec,evec,vpot,basis,wt,ylm,hp,hd,nblok
     2  ,iflag,narg,ltop,iclosed)
c
c  free-free matrices
c
      implicit real*8 (a-h, o-z)
      real*8 ylm(nblok,0:ltop,0:2*ltop)
      real*8 vpot(nblok,nchnl*(nchnl+1)/2)
      dimension nlm(nchnl),lch(lmtop,nchnl),mch(lmtop,nchnl)
      dimension iclosed(1)
      complex*16 zdotu
      complex*16 cvec(nblok),dvec(nblok),evec(nblok)
      complex*16 hpvhp(lmtop,lmtop,nchnl*(nchnl+1)/2)
      complex*16 hpvhm(lmtop,lmtop,nchnl**2)
      real*8 wt(nblok)
      complex*16 hp(nblok,0:ltop,nchnl), hd(nblok,0:ltop,nchnl)
c
c compute potential matrices
c
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
c     
c     closed channels
c     
      if(iclosed(ic).eq.1.or.iclosed(jc).eq.1) then
      hpvhp(ilm,jlm,ist)=0.0
      hpvhm(ilm,jlm,icjc)=0.0
      hpvhm(jlm,ilm,jcic)=0.0
      else
c
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
c     
c     closed channels
c     
      if(iclosed(ic).eq.1) then
      hpvhp(ilm,ilm,icc)=0.0
      hpvhm(ilm,ilm,iic)=0.0
      go to 501
      else
c
      if(iflag.eq.0) then
      do 301 i=1,narg
      cvec(i)=wt(i)*ylm(i,l1,m1)**2*hp(i,l1,ic)
      dvec(i) = conjg(hd(i,l1,ic))
301   continue
      else
      do 302 i=1,narg
      cvec(i)=wt(i)*ylm(i,l1,m1)**2*hp(i,l1,ic)
      dvec(i) = imag(hd(i,l1,ic))
302   continue
      endif
      hpvhp(ilm,ilm,icc)=hpvhp(ilm,ilm,icc)+
     1 zdotu(narg,cvec,1,hd(1,l1,ic),1)
      hpvhm(ilm,ilm,iic)=hpvhm(ilm,ilm,iic)+
     1 zdotu(narg,dvec,1,cvec,1)
      endif
501   continue
      return
      end
