*deck @(#)shells.f	1.1 9/9/91
      subroutine shells(grid,save,igrid,isave)
      parameter (numshl=50, nptdim=100)
      implicit real*8 (a-h,o-z) 
      real*8 nmulphi,nmulth
      character *20 itype, jtype, ktype
      dimension grid(4,igrid), save(isave,4)
      common/spheri/ nshell,nr(numshl),nthell(numshl), nphell(numshl),
     1               ntheta(numshl,numshl), nphi(numshl,numshl)
      common/spherr/ r(numshl),theta(numshl,numshl),phi(numshl,numshl)
      common /gauss/ x(nptdim,3),w(nptdim,3),scr(nptdim),dummy(2)
      common/io/inp,iout
      common/quadtp/ itype, jtype, ktype
      pi=3.14159265358d+00
      fac=pi/180.d+00
      do 14 i=1,numshl
      do 14 j=1,numshl
      phi(i,j)=fac*phi(i,j)
      theta(i,j)=theta(i,j)*fac
14    continue
      n1=nthell(1)+1
      n2=nphell(1)+1
      fac1=theta(n1,1)-theta(1,1)
      fac2=phi(n2,1)-phi(1,1)
      do 15 i=2,nshell
      n1=nthell(i)+1
      n2=nphell(i)+1
      fac3=theta(n1,i)-theta(1,i)
      fac4=phi(n2,i)-phi(1,i)
      if(fac1.ne.fac3.or.fac2.ne.fac4)then
      call lnkerr(' theta/phi boundary inconsistent')
      stop
      endif
15    continue
      nmulphi=2.d+00*pi/fac2
      nmulth=pi/fac1
      write(iout,16)nmulth,nmulphi
16    format(' multiplicative factors for theta and phi  :',2f10.5)
      icount=0
c
c open a loop on shells
c
       vol=0.d0
      do 1 i=1,nshell
c
c get quadrature points for this shell
c
      if(i.gt.1.and.nr(i).eq.nr(i-1))go to 11
      nrr=nr(i)
      call gaussq(itype,nrr,0.d+00,0.d+00,0,dummy,scr,x(1,1),w(1,1))
11    continue
      ampr=(r(i+1)-r(i))/2.d+00
      abpr=(r(i+1)+r(i))/2.d+00
      do 1 j=1,nthell(i)
      if(j.gt.1.and.ntheta(j,i).eq.ntheta(j-1,i))go to 12
      ntt=ntheta(j,i)
      call gaussq(jtype,ntt,0.d+00,0.d+00,0,dummy,scr,x(1,2),w(1,2))
12    continue
      ampt=(cos(theta(j,i))-cos(theta(j+1,i)))/2.d+00
      abpt=(cos(theta(j,i))+cos(theta(j+1,i)))/2.d+00
      do 1 k=1,nphell(i)
      if(k.gt.1.and.nphi(k,i).eq.nphi(k-1,i))go to 13
      npp=nphi(k,i)
      call gaussq(ktype,npp,0.d+00,0.d+00,0,dummy,scr,x(1,3),w(1,3))
13    continue
      ampp=(phi(k+1,i)-phi(k,i))/2.d+00
      abpp=(phi(k+1,i)+phi(k,i))/2.d+00
      do 51 jj=1,npp
      phii=ampp*x(jj,3)+abpp
      wphi=nmulphi*ampp*w(jj,3)
      sphi=sin(phii)
      cphi=cos(phii)
      do 51 kk=1,ntt
      ctheta=ampt*x(kk,2)+abpt
      stheta=sqrt(1-ctheta**2)
      wth=nmulth*ampt*w(kk,2)
      do 52 l=1,nrr
      zz=ampr*x(l,1)+abpr
      wz=zz*zz*ampr*w(l,1)
      icount=icount+1
      grid(1,icount)=zz*stheta*cphi
      grid(2,icount)=zz*stheta*sphi
      grid(3,icount)=zz*ctheta
      grid(4,icount)=wth*wphi*wz
      vol=vol+grid(4,icount)
   52 continue
51    continue
1     continue

      return
      end
