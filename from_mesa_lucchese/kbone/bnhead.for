      program bnhead
c
c  DIFFERENTIAL CROSS SECTIONS
c dipole transition option added 8/31/90
c
c corrected to apply MEAN approximation to scattering amplitude
c
c and handle multiple symmetry components correctly
c
c...changed 12/1/89 by t.n.r. to include phase factors in T-matrix elements
c
c  **********
c changed 6/8/90 to include scale factors for T-matrix elements
c************
c     changed 2/95 to use monder instead of spline1
c     changed to include rottational splittting for polar molecules
c************
      implicit real*8(a-h,o-z)
      parameter (nmx=500,ltop=#maxltop,mutop=#maxltop,ith=50,nspl=200)
      parameter (nsymx=8,nlmmax=72)
      parameter(nenemax=20,nchanmx=10,nrmx=50,npmx=64*nrmx)
      dimension xeu(nmx),xpeu(nmx),yeu(nmx),ypeu(nmx),zeu(nmx),zpeu(nmx)
      dimension totc(nsymx,nenemax),totcb(nsymx,nenemax)
      real*8 totanb(nenemax),totan(nenemax),kchan(nchanmx)
      dimension dummy(2),xx(64),wwt(64),scr(64)
      dimension rr(nrmx),ng(nrmx),r(npmx),wt(npmx)
      dimension aj(0:ltop),ajp(0:ltop),aplus(ltop,nenemax)
      dimension aminus(ltop,nenemax),ichan(nsymx),jchan(nsymx)
      dimension ay(0:ltop)
c..unicos
c      dimension th(ith),name(nsymx),x(64,3),w(64,3)
c..unicos
      dimension th(ith),x(64,3),w(64,3)
      character*8 name(nsymx)
c
      dimension sgni(nchanmx)
      dimension nlm(nsymx),lch(nsymx,nlmmax),mch(nsymx,nlmmax)
      dimension istart(nsymx),istop(nsymx),jstart(nsymx),jstop(nsymx)
      dimension ctk(nmx),ctkp(nmx),cpk(nmx),cpkp(nmx),spk(nmx),spkp(nmx)
      dimension cmphi(nmx),cmphip(nmx),smphi(nmx),smphip(nmx)
      dimension ww(nmx)
      dimension ylm(nmx,0:ltop,-mutop:mutop)
      dimension ylmp(nmx,0:ltop,-mutop:mutop)
      dimension p(nmx,0:ltop),pp(nmx,0:ltop)
      dimension xs(ith,nenemax),egy(nenemax),twoe(nenemax),xsec(nenemax)
      dimension xsecb(nenemax),twoep(nenemax)
      dimension tborn(nlmmax,nlmmax), xsecpol(ith,nenemax)
      dimension xsb(ith,nenemax),fac(0:60)
      dimension cspl(ith+3),dspl(ith+3),ascr(ith+3),tspl(nspl),
     $ aspl(nspl),xss(nspl)
      complex*16 fborn(nmx,nenemax),fb(nmx,nenemax)
      complex*16 ai,cmp, tmat(nlmmax,nlmmax), f(nmx,nenemax)
      ai=(0.,1.)
      data pi / 3.141592653589793  /
c
c  input of quadrature info, scattering angles on unit5
c  input of channel  info on unit8
c  input of T matrices on unit8, this unit is rewound and reread many times.
c
c      call link("unit5=inpdcs,unit6=(outdcs,hc,create),
c     1 unit7=(plthead,hc,create)//")
c
      open(5,file='inpdcs',status='unknown')
      open(6,file='outdcs',status='unknown')
      open(7,file='pltdcs',status='unknown')
c
      fac(0)=1.
      lmmax = 60
      do 61 i=1,lmmax
 61   fac(i)=i*fac(i-1)
c
c get quadratures for Euler angles alpha, beta, and gamma
c
      call pts(nalpha,nbeta,ngamma,x,w)
      write(6,159) nalpha, nbeta, ngamma
159   format(' Orders of alpha, beta, and gamma quadratures ',3i5)
      narg=nbeta*ngamma
      read(5,*) lmax,mumax
      write(6,158) lmax,mumax
158   format( ' maximum l and m are ',2i5)
c
c Is this a dipole transition?
c
      read(5,*) ifpolar
      if(ifpolar.ne.0) then
      write(6,161)
161   format(' ** dipole flag is on **')
      read(5,*) dx,dy,dz,fudge,erot
      mdip=1
      mxdip=1
      mydip=1
      dipole=sqrt(dx**2+dy**2+dz**2)
      if(abs(dz).gt.1.e-4)mdip=0
      if(abs(dx).gt.1.e-4)mxdip=0
      if(abs(dy).gt.1.e-4)mydip=0
      write(6,323) dipole,erot
 323  format(' transition moment = ',e12.5,
     $ " rotational splitting=",e12.4)
      read(5,*)ngrid
      ng1=ngrid+1
      read(5,*)(rr(i),i=1,ng1)
      read(5,*)(ng(i),i=1,ngrid)
      k=0
      nold=0
      do 15 i=1,ngrid
      amb=(rr(i+1)-rr(i))/2.
      apb=(rr(i+1)+rr(i))/2.
      nsh=ng(i)
      if(nsh.ne.nold)call gaussq(1,nsh,0.,0.,0,dummy,scr,xx,wwt)
      do 11 j=1,nsh
      k=k+1
      r(k)=amb*xx(j)+apb
      wt(k)=amb*wwt(j)
   11 continue
15    nold=nsh
      np=k
      endif
c read symmetry and channel information
      read(5,*) nsym,nchan
      write(6,77)nchan,nsym
      do 962 i=1,nsym
      read(5,*)ichan(i),jchan(i)
 962  write(6,963)i,ichan(i),jchan(i)
77    format(//" this is a",i3," channel calculation with",i3," symmetry
     1     components")
 963  format(' for symmetry',i3,' cross sections for channel',i3,
     1 ' to channel',i3)
       do 23 i=1,nsym
      read(5,750)name(i)
750   format(a8)
      write(6,78)name(i)
78    format(" symmetry  ",a8," has the following l,m pairs:")
      read(5,*) nlm(i),sgni(i)
      read(5,*) (lch(i,j),mch(i,j),j=1,nlm(i))
      write(6,79) (lch(i,j),mch(i,j),j=1,nlm(i))
79    format(2i5)
      read(5,*) istart(i),istop(i),jstart(i),jstop(i)
      write(6,76) istart(i),istop(i),jstart(i),jstop(i)
76    format(" starting and stopping ranges for channels i and j are:",
     1 4i4)
23    continue
      read(5,*) nenergy,nskip
      write(6,*)' number of energies=',nenergy,nskip
c
c read scattering angle theta
      itnum=0
 1000 read(5,*,end=999) theta
      itnum=itnum+1
      th(itnum)=theta
      theta = theta*pi/180.
      s=sin(theta)
      c=cos(theta)
      do 70 izero=1,nenemax
      xsecb(izero) = 0.
70    xsec(izero) = 0.0
c
c construct k and k prime unit vectors for quadrature pts
c in Euler angles alpha, beta, and gamma
c
      do 1 ialp=1,nalpha
      cosa=cos(x(ialp,1))
      wa=w(ialp,1)
      sina=sin(x(ialp,1))
      ii=0
      do 2 j=1,nbeta
      cosb=x(j,2)
      sinb=sqrt(1.-cosb**2)
      wb=w(j,2)
      do 3 k=1,ngamma
      ii=ii+1
      cosg=cos(x(k,3))
      sing=sin(x(k,3))
      wg=w(k,3)
      r11=cosa*cosb*cosg-sina*sing
c     r12=-cosa*cosb*sing-sina*cosg
      r13=cosa*sinb
      r21=sina*cosb*cosg+cosa*sing
c     r22=-sina*cosb*sing+cosa*cosg
      r23=sina*sinb
      r31=-sinb*cosg
c     r32=sinb*sing
      r33=cosb
      akx=r13
      aky=r23
      akz=r33
      akpx=r11*s+r13*c
      akpy=r21*s+r23*c
      akpz=r31*s+r33*c
      ctk(ii)=akz
      ctkp(ii)=akpz
      stk = sqrt(1.-ctk(ii)**2)
      stkp = sqrt( 1.-ctkp(ii)**2)
      cpk(ii) = 1.0
      spk(ii) = 0.
      if (stk.gt.1.e-15) then
      cpk(ii) = akx/stk
      spk(ii) = aky/stk
      endif
      cpkp(ii) = 1.0
      spkp(ii) = 0.0
      if(stkp.gt.1.e-15) then
      cpkp(ii) = akpx/stkp
      spkp(ii) = akpy/stkp
      endif
      ww(ii)=wa*wb*wg
      xeu(ii)=akx
      xpeu(ii)=akpx
      yeu(ii)=aky
      ypeu(ii)=akpy
      zeu(ii)=akz
      zpeu(ii)=akpz
3     continue
2     continue
c
c  we now have a block of quadrature points in beta and gamma
c for one value of alpha.
c Now we compute Ylm's for that block
c
      do 13 mu=0,mumax
      call plm(ctk,narg,mu,lmax,p)
      call plm(ctkp,narg,mu,lmax,pp)
      if(mu.eq.0)then
       do 10 i=0,lmax
       const=sqrt((2*i+1)/(4.*pi))
       do 10 j=1,narg
       ylmp(j,i,0) =  pp(j,i)*const
       ylm(j,i,0)=p(j,i)*const
10     continue
      else
       do 4 i=1,narg
       cmp=(cmplx(cpk(i),spk(i)))**mu
       cmphi(i)=dble(cmp)
       smphi(i)=dimag(cmp)
       cmp=(cmplx(cpkp(i),spkp(i)))**mu
       cmphip(i)=dble(cmp)
       smphip(i)=dimag(cmp)
4      continue
       do 12 i=mu,lmax
       const=sqrt((2*i+1)/(4.*pi)*fac(i-mu)/fac(i+mu))
c  sqrt(2) factor added to normalize "real valued" ylms
       const=const*sqrt(2.0)
       do 12 j=1,narg
       ylm(j,i,-mu)=p(j,i)*const*cmphi(j)
       ylm(j,i,mu)=p(j,i)*const*smphi(j)
       ylmp(j,i,-mu)=pp(j,i)*const*cmphip(j)
       ylmp(j,i,mu)=pp(j,i)*const*smphip(j)
12     continue
      endif
 13   continue
c
c we have the ylm's. Now read T matrices
c
c***********note**************
c we are only going to use the tmatrix block corresponding to a
c specified channel pair.
c****************************
      do 49 izero=1,nenergy
      do 49 jzero=1,narg
      fborn(jzero,izero) = 0.0
49    f(jzero,izero) = (0.0,0.0)
      do 25 isym=1,nsym
      if(ifpolar.ne.0.and.itnum.eq.1.and.ialp.eq.1)then
      write(6,"(a8)")name(isym)
      endif
c..unicos
c      call open(8,name(isym),1,len)
c..unicos
      open(8,file=name(isym),form='unformatted',status='unknown')
c
      rewind(8)
      do 24 iene=1,nenergy
66    continue
      if(iene.le.nskip) go to 24
      read(8)icw,jcw,ni,nj,aki,akj
      twoe(iene)=aki*aki
      twoep(iene)=akj*akj
      read(8)((tmat(i,j),i=1,ni),j=1,nj)
      if(ichan(isym).eq.icw.and.jchan(isym).eq.jcw)go to 67
      go to 66
67    continue
222   format(4x,6e12.5)
      degi=sgni(isym)
      do 304 i=1,ni
      do 304 j=1,nj
304   tmat(i,j)=tmat(i,j)*degi
c get total cross sections
      if(itnum.eq.1.and.ialp.eq.1)then
      totc(isym,iene)=0.
      ifirst=istart(isym)
      ilast=istop(isym)
      jfirst=jstart(isym)
      jlast=jstop(isym)
      do 71 i=ifirst,ilast
      ii=i+1-ifirst
      do 71 j=jfirst,jlast
      jj=j+1-jfirst
      totc(isym,iene)=totc(isym,iene)+cdabs(tmat(ii,jj))**2
71    continue
      endif
c
c if polar case compute tborn
c
      if(ifpolar.ne.0) then
      do 97 i=1,nlmmax
      do 97 j=1,nlmmax
97    tborn(i,j) = 0.0
      if(itnum.eq.1.and.ialp.eq.1.and.isym.eq.1)then
      call radial(lmax,np,r,wt,aki,akj,aplus(1,iene),aminus(1,iene))
      endif
      ifirst=istart(isym)
      ilast=istop(isym)
      jfirst=jstart(isym)
      jlast=jstop(isym)
      ii=0
      do 81 i=ifirst,ilast
      jj=0
      ii=ii+1
      do 81 j=jfirst,jlast
      jj=jj+1
      ldiff=lch(isym,i)-lch(isym,j)
      if(iabs(ldiff).ne.1) go to 81
c dz case
      if(mdip.eq.0)then
       if(mch(isym,i).ne.mch(isym,j)) go to 62
       if(ldiff)85,81,86
85     lg=lch(isym,j)
       rad=aplus(lg,iene)
       go to 87
86     lg=lch(isym,i)
       rad=aminus(lg,iene)
87     continue
       am=mch(isym,i)
       tborn(ii,jj) =2.*dz*rad*sqrt(aki*akj)
       tborn(ii,jj)=tborn(ii,jj)*sqrt((lg+am)*(lg-am)/(2.*lg+1.)
     1 /(2.*lg-1.))
      if(itnum.eq.1.and.ialp.eq.1)then
       write(6,867)aki,akj,ii,jj,lch(isym,i),lch(isym,j),tmat(ii,jj)
     1 ,tborn(ii,jj)
867    format(2f10.5,4i3,3e12.4)
      endif
      go to 81
c dx and dy cases
62    if(mxdip.eq.1.and.mydip.eq.1)go to 81
      else
       mi=mch(isym,i)
       mj=mch(isym,j)
       if(iabs(mi-mj).ne.1)go to 81
       li=lch(isym,i)
       lj=lch(isym,j)
       if(ldiff)64,81,65
64     lg=lj
       rad=aplus(lg,iene)
       go to 68
65     lg=li
       rad=aminus(lg,iene)
68     rad=rad*2*sqrt(aki*akj)
c case where mi is zero
        if(mi.eq.0)then
        if(ldiff)40,81,41
40      rad=rad*sqrt((lg+1.)*lg/2./(2.*lg+1.)/(2.*lg-1.))
        go to 42
41      rad=-rad*sqrt(lg*(lg-1.)/2./(2.*lg+1.)/(2.*lg-1.))
42      dip=dx
        if(mj.eq.1)dip=dy
        tborn(ii,jj)=dip*rad
        if(itnum.eq.1.and.ialp.eq.1)then
        write(6,867)aki,akj,ii,jj,lch(isym,i),lch(isym,j),tmat(ii,jj)
     1  ,tborn(ii,jj)
      endif
        go to 81
       endif
c case where mj is zero
        if(mj.eq.0)then
        if(ldiff)43,81,44
43      rad=-rad*sqrt(lg*(lg-1.)/2./(2.*lg+1.)/(2.*lg-1.))
        go to 45
44      rad=rad*sqrt((lg+1.)*lg/2./(2.*lg+1.)/(2.*lg-1.))
45      dip=dx
        if(mi.eq.1)dip=dy
        tborn(ii,jj)=dip*rad
        if(itnum.eq.1.and.ialp.eq.1)then
        write(6,867)aki,akj,ii,jj,lch(isym,i),lch(isym,j),tmat(ii,jj)
     1  ,tborn(ii,jj)
      endif
        go to 81
      endif
c mi and mj are both nonzero
       is=mi/iabs(mi)
       js=mj/iabs(mj)
c dx case
        if(is.eq.js)then
        is=mi*is
        js=mj*js
        if((is+1).eq.js)then
         mi=is
         mj=js
         sgn=1.
        else
         mi=-is
         mj=-js
         sgn=-1.
        endif
       if(ldiff)46,81,47
46     rad=rad*sgn*sqrt((lg+mi+1.)*(lg+mi)/4./(2.*lg+1.)/(2.*lg-1.))
       go to 48
47     rad=-rad*sgn*sqrt((lg-mj+1.)*(lg-mj)/4./(2.*lg+1.)/(2.*lg-1.))
48     tborn(ii,jj)=dx*rad
        if(itnum.eq.1.and.ialp.eq.1)then
        write(6,867)aki,akj,ii,jj,lch(isym,i),lch(isym,j),tmat(ii,jj)
     1  ,tborn(ii,jj)
      endif
       go to 81
       else
c dy cases
        write(6,*)'not yet'
       endif
      endif
81    continue
c get total born cross sections
c     if(itnum.eq.1.and.ialp.eq.1.and.ichan(isym).ne.jchan(isym))then
      if(itnum.eq.1.and.ialp.eq.1)then
      totcb(isym,iene)=0.
      ifirst=istart(isym)
      ilast=istop(isym)
      jfirst=jstart(isym)
      jlast=jstop(isym)
      do 63 i=ifirst,ilast
      ii=i+1-ifirst
      do 63 j=jfirst,jlast
      jj=j+1-jfirst
      totcb(isym,iene)=totcb(isym,iene)+abs(tborn(ii,jj))**2
63    continue
      endif
      endif
c
c
c accumulate scattering amplitude for a particular energy
c
      ifirst=istart(isym)
      ilast=istop(isym)
      jfirst=jstart(isym)
      jlast=jstop(isym)
      do 31 i=ifirst,ilast
      ii=i+1-ifirst
      li = lch(isym,i)
      mi = mch(isym,i)
      do 31 j=jfirst,jlast
      jj=j+1-jfirst
      lj = lch(isym,j)
      mj = mch(isym,j)
      do 32 iarg=1,narg
      f(iarg,iene) = ylm(iarg,li,mi)*tmat(ii,jj)*ylmp(iarg,lj,mj)
     1 *ai**(li-lj)
     #  + f(iarg,iene)
32    continue
c
c if polar case, accumulate born amplitude
c
      if(ifpolar.ne.0) then
      do 83 iarg=1,narg
      fborn(iarg,iene) = ylm(iarg,li,mi)*tborn(ii,jj)*ylmp(iarg,lj,mj)
     1 *ai**(li-lj)
     #  + fborn(iarg,iene)
83    continue
      endif
c
31    continue
24    continue
c      call close(8)
      close(8)
25    continue
c
c if polar case, compute analytic born amplitude
c
      if(ifpolar.ne.0)then
      do 54 iene=1,nenergy
         if (iene.le.nskip) go to 54
      ak=sqrt(twoe(iene))
      akp=sqrt(twoep(iene)-2.*erot)
      qmag=(ak*ak+akp*akp-2.*c*ak*akp)
      pre=sqrt(ak*akp/qmag)/2./pi
      qmag=sqrt(qmag)
      do 554 iarg=1,narg
      zdif=dz*(zpeu(iarg)*akp-zeu(iarg)*ak)
     1 +dx*(xpeu(iarg)*akp-xeu(iarg)*ak)
     1 +dy*(ypeu(iarg)*akp-yeu(iarg)*ak)
      zdif=zdif/qmag
      fb(iarg,iene)=-fudge*pre*ai*zdif
554    continue
54    continue
      endif
c
c accumulate cross section
c
      if(ifpolar.eq.0)then
      do 50 iene=1,nenergy
      do 50 iarg=1,narg
      xsec(iene)=xsec(iene)+ww(iarg)*cdabs(f(iarg,iene))**2
50    continue
      else
      do 51 iene=1,nenergy
      do 51 iarg=1,narg
      xsec(iene)=xsec(iene)+ww(iarg)*cdabs(f(iarg,iene)
     1 +fb(iarg,iene)-fborn(iarg,iene))**2
      xsecb(iene)=xsecb(iene)+ww(iarg)*cdabs(fborn(iarg,iene))**2
51    continue
      endif
c
c close loop on alpha quadrature blocks
1     continue
      write(6,157) theta
157   format(' scattering angle =',f15.8)
      do 80 iene=1,nenergy
c
c  the 1/(8*pi**2) is the factor for the average over euler angles
c  the 16*pi**2 is a factor of 4*pi in the T-matrix definition
c
      xsec(iene)=xsec(iene)/(8.*pi**2)/twoe(iene)
      xsec(iene) = xsec(iene)*16.*pi**2
      einc = twoe(iene)/2.
      write(6,101) einc,xsec(iene)
  101 format(' incident e =',e12.5,'  cross section =',e12.5)
      egy(iene)=einc
      xs(itnum,iene)=xsec(iene)
80    continue
c
c if polar case, finish computing xsec from born t-matrix
c
      if(ifpolar.ne.0) then
      do 91 iene=1,nenergy
      einc = twoe(iene)/2.
      xsecb(iene)=xsecb(iene)/(8.*pi**2)/twoe(iene)
      xsecb(iene) = xsecb(iene)*16.*pi**2
      write(6,324) einc,xsecb(iene)
      xsb(itnum,iene)=xsecb(iene)
  324 format(' incident e =',e12.5,'  born xsec from T =',e12.5)
91    continue
      endif
      go to 1000
  999 continue
      do 409 i=1,nenergy
      totan(i)=0.
      do 409 j=1,nsym
      totan(i)=totan(i)+4.*pi*totc(j,i)/twoe(i)
409   continue
      if(ifpolar.ne.0)then
      do 427 i=1,nenergy
      totanb(i)=0.
      do 427 j=1,nsym
      totanb(i)=totanb(i)+4.*pi*totcb(j,i)/twoe(i)
427   continue
      endif
c
      scalbhl=.529*.529
c
      ispl=0
      do 972 i=200,1,-1
         ispl=ispl+1
         aspl(ispl)=-1.+(i-1)*.010050251
         tspl(ispl)=acos(aspl(ispl))*57.29577951
 972  continue
      do 997 i=1,nenergy
      write(7,998)egy(i),(th(j),xs(j,i),j=1,itnum)
      call monder(itnum,th,xs(1,i),dspl,ascr,cspl)
      call pwcfev(1,itnum,th,xs(1,i),dspl,200,tspl,xss)
      tot=.5*(xss(1)+xss(200))
      totm=xss(200)
c      write(6,973)tspl(1),tot,totm
      do 765 ispl=2,199
         spll=xss(ispl)
      spill=spll*(1.-aspl(ispl))
      tot=tot+spll
      totm=totm+spill
c      write(6,973)tspl(ispl),spll,spill
 973  format(3e16.8)
765   continue
      tot=tot*.063146
      totm=totm*.063146
      anal=totan(i)
      if(ifpolar.eq.0)then
      write(6,764)egy(i),tot,anal,scalbhl*anal
      else
c      if(ichan(1).ne.jchan(1))then
      ak=sqrt(twoe(i))
      akp=sqrt(twoep(i))
         if(ichan(1).eq.jchan(1))akp=sqrt(twoe(i)-2.*erot)
         fulborn=8.*pi/3./twoe(i)*dipole**2*dlog((ak+akp)
     1 /(ak-akp))*fudge**2
      anal=anal-totanb(i)+fulborn
      write(6,764)egy(i),tot,anal,scalbhl*anal
      write(6,864)totan(i),fulborn,totanb(i)
 864  format(' total from T=',e12.4,' 3-D Born=',e12.4,
     $ ' partial wave Born=',e12.4)
c      else
c      write(6,766)egy(i),tot,scalbhl*tot
c      endif
      endif
      write(6,767)totm
 767  format("     momentum transfer cross section =",e12.4)
764   format(" energy =",f10.5,"  total cross section=",2e12.4,
     $ " total (angstroms**2) ",e12.4)
766   format(" energy =",f10.5,"  total cross section=",e12.4,
     $ " total (angstroms**2) ",e12.4)
997   continue
998   format(e16.8,/,(2f20.10))
      call exit
      end
