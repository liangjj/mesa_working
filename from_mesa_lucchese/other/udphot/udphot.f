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
      parameter (nmx=500,ltop=7,mutop=5,ith=50)
      parameter (nsymx=8,nlmmax=72)
      parameter(nenemax=20,nchanmx=10,nrmx=50,npmx=64*nrmx)
      dimension xeu(nmx),xpeu(nmx),yeu(nmx),ypeu(nmx),zeu(nmx),zpeu(nmx)
      dimension totc(nsymx,nenemax),totcb(nsymx,nenemax)
      real totanb(nenemax),totan(nenemax),kchan(nchanmx)
      dimension dummy(2),xx(64),wwt(64),scr(64),beeta(ith)
      dimension rr(nrmx),ng(nrmx),r(npmx),wt(npmx)
      dimension aj(0:ltop),ajp(0:ltop),aplus(ltop,nenemax)
      dimension aminus(ltop,nenemax)
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
      dimension cspl(ith+3),dspl(ith+3)
      complex fborn(nmx,nenemax),fb(nmx,nenemax)
      complex ai,cmp, tmat(nlmmax,nlmmax), f(nmx,nenemax)
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
      open(5,file='inpdcs')
      open(6,file='outdcs')
      open(7,file='pltdcs')
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
      read(5,*) lmax,mumax,iprt
      write(6,158) lmax,mumax
158   format( ' maximum l and m are ',2i5)
c
c Is this a dipole transition?
c
      read(5,*) ifpolar
      if(ifpolar.ne.0) then
      write(6,161)
161   format(' ** dipole flag is on **')
      read(5,*) dx,dy,dz,fudge
      mdip=1
      mxdip=1
      mydip=1
      dipole=sqrt(dx**2+dy**2+dz**2)
      if(abs(dz).gt.1.e-4)mdip=0
      if(abs(dx).gt.1.e-4)mxdip=0
      if(abs(dy).gt.1.e-4)mydip=0
      write(6,323) dipole
323   format(' transition moment = ',e12.5)
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
      read(5,*) nsym,nchan,ichan,jchan
      write(6,77)nchan,nsym,ichan,jchan
77    format(//" this is a",i3," channel calculation with",i3," symmetry
     1 components"/" cross sections for going from channel"i3,
     1 " to channel",i3)
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
      read(5,*) nenergy
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
       cmphi(i)=real(cmp)
       smphi(i)=aimag(cmp)
       cmp=(cmplx(cpkp(i),spkp(i)))**mu
       cmphip(i)=real(cmp)
       smphip(i)=aimag(cmp)
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
      open(8,file=name(isym),form='unformatted')
c
      rewind(8)
      do 24 iene=1,nenergy
66    continue
      read(8)icw,jcw,ni,nj,aki,akj
      twoe(iene)=aki*aki
      twoep(iene)=akj*akj
      read(8)((tmat(i,j),i=1,ni),j=1,nj)
      if(ichan.eq.icw.and.jchan.eq.jcw)go to 67
      go to 66
67    continue
222   format(4x,6e12.5)
      degi=sgni(isym)
      do 304 i=1,ni
      do 304 j=1,nj
304   tmat(i,j)=tmat(i,j)*degi
c get total cross sections
      if(itnum.eq.1.and.ialp.eq.1)then
c
c..photo
      if(iprt.gt.0) then
       write(6,*)' Tmatrix ',icw,jcw
      write(6,88001)((tmat(i,j),i=1,ni),j=1,nj)
88001 format(3(4x,f12.8,2x,f12.8))
      end if 
c..photo
c
      totc(isym,iene)=0.
      ifirst=istart(isym)
      ilast=istop(isym)
      jfirst=jstart(isym)
      jlast=jstop(isym)
      do 71 i=ifirst,ilast
      ii=i+1-ifirst
      do 71 j=jfirst,jlast
      jj=j+1-jfirst
      totc(isym,iene)=totc(isym,iene)+cabs(tmat(ii,jj))**2
71    continue
c
      if(iprt.gt.0) then
      write(6,*)' isym iene totc ',isym,iene,totc(isym,iene)
      end if
c
      end if
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
      if(itnum.eq.1.and.ialp.eq.1.and.ichan.ne.jchan)then
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
     #  + f(iarg,iene)
32    continue
c
c..photo     1 *ai**(li-lj)
c
c if polar case, accumulate born amplitude
c
      if(ifpolar.ne.0) then
      do 83 iarg=1,narg
      fborn(iarg,iene) = ylm(iarg,li,mi)*tborn(ii,jj)*ylmp(iarg,lj,mj)
     #  + fborn(iarg,iene)
83    continue
c
c..photo     1 *ai**(li-lj)
c
      endif
c
31    continue
24    continue
      call close(8)
25    continue
c
c if polar case, compute analytic born amplitude
c
      if(ifpolar.ne.0)then
      do 54 iene=1,nenergy
      ak=sqrt(twoe(iene))
      akp=sqrt(twoep(iene))
      qmag=(ak*ak+akp*akp-2.*c*ak*akp)
      pre=sqrt(ak*akp/qmag)/2./pi
      qmag=sqrt(qmag)
      do 54 iarg=1,narg
      zdif=dz*(zpeu(iarg)*akp-zeu(iarg)*ak)
     1 +dx*(xpeu(iarg)*akp-xeu(iarg)*ak)
     1 +dy*(ypeu(iarg)*akp-yeu(iarg)*ak)
      zdif=zdif/qmag
      fb(iarg,iene)=-fudge*pre*ai*zdif
54    continue
      endif
c
c accumulate cross section
c
      if(ifpolar.eq.0)then
      do 50 iene=1,nenergy
      do 50 iarg=1,narg
      xsec(iene)=xsec(iene)+ww(iarg)*cabs(f(iarg,iene))**2
50    continue
      else
      do 51 iene=1,nenergy
      do 51 iarg=1,narg
      xsec(iene)=xsec(iene)+ww(iarg)*cabs(f(iarg,iene)
     1 +fb(iarg,iene)-fborn(iarg,iene))**2
      xsecb(iene)=xsecb(iene)+ww(iarg)*cabs(fborn(iarg,iene))**2
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
c
         xsec(iene)=xsec(iene)/(2.*pi)
c
c..photo
c      xsec(iene)=xsec(iene)/(8.*pi**2)
c
c      xsec(iene) = xsec(iene)*16.*pi**2
c..c photo
c..c
c      xsec(iene) = xsec(iene)*16.*pi**2
c..c photo
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
cc      totan(i)=totan(i)+4.*pi*totc(j,i)/twoe(i)
c..photo
      totan(i)=totan(i)+totc(j,i)
c..photo
409   continue
cc
      if(ifpolar.ne.0)then
      do 427 i=1,nenergy
      totanb(i)=0.
      do 427 j=1,nsym
      totanb(i)=totanb(i)+4.*pi*totcb(j,i)/twoe(i)
427   continue
      endif
cc
c
      scalbhl=.529*.529
c
      do 997 i=1,nenergy
      if(totan(i).le.1.e-5)then
            do 9994 j=1,itnum
 9994       beeta(j)=0.
      else
      do 9992 j=1,itnum
      theta=th(j)*pi/180.
      xx2=cos(theta)
      pp2=.5*(3.*xx2*xx2-1.)
      beeta(j)=(xs(j,i)*4.*pi/totan(i)-1.)/pp2
 9992 continue
      endif
      write(7,998)egy(i),(th(j),beeta(j),xs(j,i),j=1,itnum)
      tot=spline1(th,xs(1,i),tspl,itnum,cspl,dspl,1)
      tot=0.
      tspl=180.
      tot=tot+.5*spline1(th,xs(1,i),tspl,itnum,cspl,dspl,2)
c..      totm=2.*tot
      tspl=0.
      tot=tot+.5*spline1(th,xs(1,i),tspl,itnum,cspl,dspl,2)
      do 765 ispl=2,199
      aspl=-1.+(ispl-1)*.01005
      tspl=acos(aspl)*57.29577951
      spll=spline1(th,xs(1,i),tspl,itnum,cspl,dspl,2)
      tot=tot+spll
c..      totm=totm+spll*(1.-aspl)
 765  continue
c..photo
      tot=tot*.063146
c..photo
c..      totm=totm*.063146
      anal=totan(i)
cc
      if(ifpolar.eq.0)then
      write(6,764)egy(i),tot,anal,scalbhl*anal
      else
c
      if(ichan.ne.jchan)then
      ak=sqrt(twoe(i))
      akp=sqrt(twoep(i))
      anal=anal-totanb(i)+8.*pi/3./twoe(i)*dipole**2*alog((ak+akp)
     1 /(ak-akp))*fudge**2
      write(6,764)egy(i),tot,anal,scalbhl*anal
      else
      write(6,766)egy(i),tot,scalbhl*tot
      endif
c
      endif
cc
c..      write(6,767)totm
 767  format("     momentum transfer cross section =",e12.4)
764   format(" energy =",f10.5,"  total cross section=",2e12.4,
     $ " total (angstroms**2) ",e12.4)
766   format(" energy =",f10.5,"  total cross section=",e12.4,
     $ " total (angstroms**2) ",e12.4)
997   continue
998   format(e16.8,/,(3f20.10))
      call exit
      end
      subroutine radial(lmax,np,r,wt,ak,akp,aplus,aminus)
      dimension aplus(1),aminus(1),r(1),wt(1)
      dimension aj(0:20),ajp(0:20),ay(0:20)
      do 1 i=1,np
      x1=ak*r(i)
      x2=akp*r(i)
      call sjymec(x1,aj,ay,ncal,lmax)
      call sjymec(x2,ajp,ay,ncal,lmax)
      do 1 l=1,lmax
      aplus(l)=aplus(l)+aj(l-1)*ajp(l)*wt(i)
      aminus(l)=aminus(l)+aj(l)*ajp(l-1)*wt(i)
1     continue
c     write(6,66)ak,akp
66    format(" ak,akp:",2e14.4)
c     do 12 i=1,lmax
c     write(6,100)i,aplus(i),aminus(i)
12    continue
100   format(i5,2e12.4)
      return
      end
      function spline1(x,y,z,nn,c,d,isw)
      dimension x(1),y(1),c(1),d(1),error(2)
    1 format(24h0** error in spline1... 2a8/)
      error(1)=8hunordere
      error(2)=8hd x-vals
    2 n=nn
      np1=n+1
      np2=n+2
      go to (4,19),isw
    3 return
c     section 4 calculates the spline coefficients c(1) is the
c     constant term, c(2) is the coeff of the linear term, and
c     c(3) thru c(n+2) are the spline coefficients.
    4 c(1)=y(1)
      d(1)=1.0
      c(np1)=0.0
      d(np1)=0.0
      c(np2)=0.0
      d(np2)=0.0
      do 5 i=2,n
      c(i)=y(i)-y(1)
    5 d(i)=x(i)-x(1)
      do 14 i=3,np2
      pivot=1.0/d(i-1)
      if(i.lt.np2) 6,9
    6 supd=x(i-1)-x(i-2)
      if(supd.lt.0.0) 7,8
    7 write(6,1) error(1),error(2)
      spline1=0.0
      return
    8 supd=supd**3
      go to 10
    9 supd=1.0
   10 dfact=supd*pivot
      cfact=c(i-1)*pivot
      if(i.gt.n) go to 30
      do 11 j=i,n
      v=x(j)-x(i-2)
      c(j)=c(j)-d(j)*cfact
   11 d(j)=(v**3)-d(j)*dfact
   30 if(i.lt.np2) 12,13
   12 c(np1)=c(np1)-d(np1)*cfact
      d(np1)=1.0-d(np1)*dfact
   13 c(np2)=c(np2)-d(np2)*cfact
   14 d(np2)=x(i-2)-d(np2)*dfact
      do 18 i=1,n
      j=np2-i
      if(j.eq.np1) 15,16
   15 v=1.0
      go to 17
   16 v=x(j)-x(j-1)
      v=v**3
   17 c(j+1)=c(j+1)/d(j+1)
   18 c(j)=c(j)-c(j+1)*v
      c(2)=c(2)/d(2)
      go to 22
c     section 19 evaluates the spline function
   19 spline1=c(1)+c(2)*(z-x(1))
      do 20 i=1,n
      v=z-x(i)
      if(v.gt.0.0) 20,21
   20 spline1=spline1+c(i+2)*(v**3)
   21 return
   22 spline1=0.0
      return
      end
      subroutine pts(nalpha,nbeta,ngamma,x,w)
      real nmula,nmulb,nmulg
      dimension xx(64),ww(64),x(64,3),w(64,3)
     1,dummy(2),scr(64)
      read(5,*)itype,jtype,ktype
      read(5,*)nalpha,nbeta,ngamma
      read(5,*)a1,a2
      read(5,*)b1,b2
      read(5,*)g1,g2
      pi=3.14159265358
      fac=pi/180.
      a1=fac*a1
      a2=fac*a2
      b1=fac*b1
      b2=fac*b2
      g1=fac*g1
      g2=fac*g2
      faca=a2-a1
      facb=b2-b1
      facg=g2-g1
      nmula=2.*pi/faca
      nmulb=pi/facb
      nmulg=2*pi/facg
      write(6,16)nmula,nmulb,nmulg
16    format(" multiplicative factors for euler angles  :",3f10.5)
      irun=0
      call gaussq(itype,nalpha,0.,0.,0,dummy,scr,xx,ww)
      ampa=(a2-a1)/2.
      abpa=(a2+a1)/2.
      do 1 i=1,nalpha
      x(i,1)=ampa*xx(i)+abpa
1     w(i,1)=nmula*ampa*ww(i)
      call gaussq(jtype,nbeta,0.,0.,0,dummy,scr,xx,ww)
      ampb=(cos(b1)-cos(b2))/2.
      abpb=(cos(b2)+cos(b1))/2.
      do 2 i=1,nbeta
      x(i,2)=ampb*xx(i)+abpb
2     w(i,2)=nmulb*ampb*ww(i)
      call gaussq(ktype,ngamma,0.,0.,0,dummy,scr,xx,ww)
      ampg=(g2-g1)/2.
      abpg=(g2+g1)/2.
      do 3 i=1,ngamma
      x(i,3)=ampg*xx(i)+abpg
3     w(i,3)=nmulg*ampg*ww(i)
      return
      end
      subroutine gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w)
c
c           this set of routines computes the nodes x(i) and weights
c        c(i) for gaussian-type quadrature rules with pre-assigned
c        nodes.  these are used when one wishes to approximate
c
c                 integral (from a to b)  f(x) w(x) dx
c
c                              n
c        by                   sum c  f(x )
c                             i=1  i    i
c
c        here w(x) is one of six possible non-negative weight
c        functions (listed below), and f(x) is the
c        function to be integrated.  gaussian quadrature is particularly
c        useful on infinite intervals (with appropriate weight
c        functions), since then other techniques often fail.
c
c           associated with each weight function w(x) is a set of
c        orthogonal polynomials.  the nodes x(i) are just the zeroes
c        of the proper n-th degree polynomial.
c
c     input parameters
c
c        kind     an integer between 0 and 6 giving the type of
c                 quadrature rule
c
c        kind = 0=  simpson's rule w(x) = 1 on (-1, 1) n must be odd.
c        kind = 1=  legendre quadrature, w(x) = 1 on (-1, 1)
c        kind = 2=  chebyshev quadrature of the first kind
c                   w(x) = 1/dsqrt(1 - x*x) on (-1, +1)
c        kind = 3=  chebyshev quadrature of the second kind
c                   w(x) = dsqrt(1 - x*x) on (-1, 1)
c        kind = 4=  hermite quadrature, w(x) = exp(-x*x) on
c                   (-infinity, +infinity)
c        kind = 5=  jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**
c                   beta on (-1, 1), alpha, beta .gt. -1.
c                   note= kind=2 and 3 are a special case of this.
c        kind = 6=  generalized laguerre quadrature, w(x) = exp(-x)*
c                   x**alpha on (0, +infinity), alpha .gt. -1
c
c        n        the number of points used for the quadrature rule
c        alpha    real parameter used only for gauss-jacobi and gauss-
c                 laguerre quadrature (otherwise use 0.).
c        beta     real parameter used only for gauss-jacobi quadrature--
c                 (otherwise use 0.).
c        kpts     (integer) normally 0, unless the left or right end-
c                 point (or both) of the interval is required to be a
c                 node (this is called gauss-radau or gauss-lobatto
c                 quadrature).  then kpts is the number of fixed
c                 endpoints (1 or 2).
c        endpts   real array of length 2.  contains the values of
c                 any fixed endpoints, if kpts = 1 or 2.
c        b        real scratch array of length n
c
c     output parameters (both arrays of length n)
c
c        t        will contain the desired nodes x(1),,,x(n)
c        w        will contain the desired weights c(1),,,c(n)
c
c     subroutines required
c
c        gbslve, class, and gbtql2 are provided. underflow may sometimes
c        occur, but it is harmless if the underflow interrupts are
c        turned off as they are on this machine.
c
c     accuracy
c
c        the routine was tested up to n = 512 for legendre quadrature,
c        up to n = 136 for hermite, up to n = 68 for laguerre, and up
c        to n = 10 or 20 in other cases.  in all but two instances,
c        comparison with tables in ref. 3 showed 12 or more significant
c        digits of accuracy.  the two exceptions were the weights for
c        hermite and laguerre quadrature, where underflow caused some
c        very small weights to be set to zero.  this is, of course,
c        completely harmless.
c
c     method
c
c           the coefficients of the three-term recurrence relation
c        for the corresponding set of orthogonal polynomials are
c        used to form a symmetric tridiagonal matrix, whose
c        eigenvalues (determined by the implicit ql-method with
c        shifts) are just the desired nodes.  the first components of
c        the orthonormalized eigenvectors, when properly scaled,
c        yield the weights.  this technique is much faster than using a
c        root-finder to locate the zeroes of the orthogonal polynomial.
c        for further details, see ref. 1.  ref. 2 contains details of
c        gauss-radau and gauss-lobatto quadrature only.
c
c     references
c
c        1.  golub, g. h., and welsch, j. h.,  calculation of gaussian
c            quadrature rules,  mathematics of computation 23 (april,
c            1969), pp. 221-230.
c        2.  golub, g. h.,  some modified matrix eigenvalue problems,
c            siam review 15 (april, 1973), pp. 318-334 (section 7).
c        3.  stroud and secrest, gaussian quadrature formulas, prentice-
c            hall, englewood cliffs, n.j., 1966.
c
c     ..................................................................
c
      implicit real*8 (a-h,o-z)
      real*8  muzero
      dimension  b(n),t(n),w(n),endpts(2)
      if(kind.eq.0) then
       if(2*(n/2).eq.n) then
        write(6,800) n
800     format(" n must be odd for simpson's rule ",i5)
        stop
      endif
        if(n.le.1) then
        t(1) = 0.
        w(1) = 2.0
        return
      endif
       h = 2.0/(n-1)
       t(1) = -1.0
       t(n) = 1.0
       w(1) = h/3.0
       w(n) = h/3.0
       nm1 = n-1
       do 801 i=2,nm1
       t(i) = t(i-1) + h
       w(i) = 4.0 - 2.0*(i-2*(i/2))
801    w(i) = w(i)*h/3.0
       return
      endif
c
      call class (kind, n, alpha, beta, b, t, muzero)
c
c           the matrix of coefficients is assumed to be symmetric.
c           the array t contains the diagonal elements, the array
c           b the off-diagonal elements.
c           make appropriate changes in the lower right 2 by 2
c           submatrix.
c
      if (kpts.eq.0)  go to 100
      if (kpts.eq.2)  go to  50
c
c           if kpts=1, only t(n) must be changed
c
      t(n) =gbslve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)
      go to 100
c
c           if kpts=2, t(n) and b(n-1) must be recomputed
c
   50 gam =gbslve(endpts(1), n, t, b)
      t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, t, b) - gam))
      b(n-1) =  dsqrt(t1)
      t(n) = endpts(1) + gam*t1
c
c           note that the indices of the elements of b run from 1 to n-1
c           and thus the value of b(n) is arbitrary.
c           now compute the eigenvalues of the symmetric tridiagonal
c           matrix, which has been modified as necessary.
c           the method used is a ql-type method with origin shifting
c
  100 w(1) = 1.0d0
      do 105 i = 2, n
  105    w(i) = 0.0d0
c
      call gbtql2 (n, t, b, w, ierr)
      do 110 i = 1, n
  110    w(i) = muzero * w(i) * w(i)
c
      return
      end
c
c
c
      function gbslve(shift, n, a, b)
c
c       this procedure performs elimination to solve for the
c       n-th component of the solution delta to the equation
c
c             (jn - shift*identity) * delta  = en,
c
c       where en is the vector of all zeroes except for 1 in
c       the n-th position.
c
c       the matrix jn is symmetric tridiagonal, with diagonal
c       elements a(i), off-diagonal elements b(i).  this equation
c       must be solved to obtain the appropriate changes in the lower
c       2 by 2 submatrix of coefficients for orthogonal polynomials.
c
c
      implicit real*8 (a-h,o-z)
      dimension  a(n),b(n)
c
      alpha = a(1) - shift
      nm1 = n - 1
      do 10 i = 2, nm1
   10    alpha = a(i) - shift - b(i-1)**2/alpha
      gbslve = 1.0d0  /alpha
      return
      end
c
c
c
      subroutine class(kind, n, alpha, beta, b, a, muzero)
c
c           this procedure supplies the coefficients a(j), b(j) of the
c        recurrence relation
c
c             b p (x) = (x - a ) p   (x) - b   p   (x)
c              j j            j   j-1       j-1 j-2
c
c        for the various classical (normalized) orthogonal polynomials,
c        and the zero-th moment
c
c             muzero = integral w(x) dx
c
c        of the given polynomial   weight function w(x).  since the
c        polynomials are orthonormalized, the tridiagonal matrix is
c        guaranteed to be symmetric.
c
c           the input parameter alpha is used only for laguerre and
c        jacobi polynomials, and the parameter beta is used only for
c        jacobi polynomials.  the laguerre and jacobi polynomials
c        require the gamma function.
c
c     ..................................................................
c
      implicit real*8 (a-h,o-z)
      dimension  a(n),b(n)
      real*8  muzero
      data pi / 3.141592653589793d0  /
c
      nm1 = n - 1
      go to (10, 20, 30, 40, 50, 60), kind
c
c              kind = 1=  legendre polynomials p(x)
c              on (-1, +1), w(x) = 1.
c
   10 muzero = 2.0d0
      do 11 i = 1, nm1
         a(i) = 0.0d0
         abi = i
   11    b(i) = abi/ dsqrt(4*abi*abi - 1.0d0  )
      a(n) = 0.0d0
      return
c
c              kind = 2=  chebyshev polynomials of the first kind t(x)
c              on (-1, +1), w(x) = 1 / dsqrt(1 - x*x)
c
   20 muzero = pi
      do 21 i = 1, nm1
         a(i) = 0.0d0
   21    b(i) = 0.5d0
      b(1) =  dsqrt(0.5d0  )
      a(n) = 0.0d0
      return
c
c              kind = 3=  chebyshev polynomials of the second kind u(x)
c              on (-1, +1), w(x) = dsqrt(1 - x*x)
c
   30 muzero = pi/2.0d0
      do 31 i = 1, nm1
         a(i) = 0.0d0
   31    b(i) = 0.5d0
      a(n) = 0.0d0
      return
c
c              kind = 4=  hermite polynomials h(x) on (-infinity,
c              +infinity), w(x) = exp(-x**2)
c
   40 muzero =  dsqrt(pi)
      do 41 i = 1, nm1
         a(i) = 0.0d0
   41    b(i) =  dsqrt(i/2.0d0  )
      a(n) = 0.0d0
      return
c
c              kind = 5=  jacobi polynomials p(alpha, beta)(x) on
c              (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and
c              beta greater than -1
c
   50 ab = alpha + beta
      abi = 2.0d0   + ab
      muzero = 2.0d0   ** (ab + 1.0d0  ) * gamfun(alpha + 1.0d0  ) * gam
     vfun(
     x beta + 1.0d0  ) / gamfun(abi)
      a(1) = (beta - alpha)/abi
      b(1) =  dsqrt(4.0d0  *(1.0d0  + alpha)*(1.0d0   + beta)/((abi + 1.
     v0d0  )*
     1  abi*abi))
      a2b2 = beta*beta - alpha*alpha
      do 51 i = 2, nm1
         abi = 2.0d0  *i + ab
         a(i) = a2b2/((abi - 2.0d0  )*abi)
   51    b(i) =  dsqrt (4.0d0  *i*(i + alpha)*(i + beta)*(i + ab)/
     1   ((abi*abi - 1)*abi*abi))
      abi = 2.0d0  *n + ab
      a(n) = a2b2/((abi - 2.0d0  )*abi)
      return
c
c              kind = 6=  laguerre polynomials l(alpha)(x) on
c              (0, +infinity), w(x) = exp(-x) * x**alpha, alpha greater
c              than -1.
c
   60 muzero = gamfun(alpha + 1.0d0  )
      do 61 i = 1, nm1
         a(i) = 2.0d0  *i - 1.0d0   + alpha
   61    b(i) =  dsqrt(i*(i + alpha))
      a(n) = 2.0d0  *n - 1 + alpha
      return
      end
c     ------------------------------------------------------------------
c
      subroutine gbtql2(n, d, e, z, ierr)
c
c     this subroutine is a translation of the algol procedure imtql2,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues and first components of the
c     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
c     method, and is adapted from the eispak routine imtql2
c
c     on input=
c
c        n is the order of the matrix;
c
c        d contains the diagonal elements of the input matrix;
c
c        e contains the subdiagonal elements of the input matrix
c          in its first n-1 positions.  e(n) is arbitrary;
c
c        z contains the first row of the identity matrix.
c
c      on output=
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1, 2, ..., ierr-1;
c
c        e has been destroyed;
c
c        z contains the first components of the orthonormal eigenvectors
c          of the symmetric tridiagonal matrix.  if an error exit is
c          made, z contains the eigenvectors associated with the stored
c          eigenvalues;
c
c        ierr is set to
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     ------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      integer i, j, k, l, m, n, ii, mml, ierr
      real*8  machep
      dimension  d(n),e(n),z(n)
c
c     ========== machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c                machep = 16.0d0**(-13) for long form arithmetic
c                on s360 ==========
       machep=1.0e-14
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
c     ========== look for small sub-diagonal element ==========
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            if ( dabs(e(m)) .le. machep * ( dabs(d(m)) +  dabs(d(m+1))))
     x         go to 120
  110    continue
c
  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
c     ========== form shift ==========
         g = (d(l+1) - p) / (2.0d0   * e(l))
         r =  dsqrt(g*g+1.0d0  )
         g = d(m) - p + e(l) / (g +  dsign(r, g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c     ========== for i=m-1 step -1 until l do -- ==========
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if ( dabs(f) .lt.  dabs(g)) go to 150
            c = g / f
            r =  dsqrt(c*c+1.0d0  )
            e(i+1) = f * r
            s = 1.0d0   / r
            c = c * s
            go to 160
  150       s = f / g
            r =  dsqrt(s*s+1.0d0  )
            e(i+1) = g * r
            c = 1.0d0   / r
            s = s * c
  160       g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0   * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
c     ========== form first component of vector ==========
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
            z(i) = c * z(i) - s * f
c
  200    continue
c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
  240 continue
c     ========== order eigenvalues and eigenvectors ==========
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         p = z(i)
         z(i) = z(k)
         z(k) = p
c
  300 continue
c
      go to 1001
c     ========== set error -- no convergence to an
c                eigenvalue after 30 iterations ==========
 1000 ierr = l
 1001 return
c     ========== last card of gbtql2 ==========
      end
      function  gamfun(z)
c  this is a procedure that evaluates gamma(z) for
c     0 lt z le 3 to 16 significant figures
c    it is based on a chebyshev-type polynomial
c   approximation given in h. werner and r. collinge, math. comput.
c    15 (1961), pp. 195-97.
c   approximations to the gamma function, accurate up to 18 significant
c   digits, may be found in the paper quoted above
c
c
c
      implicit real*8 (a-h,o-z)
      dimension  a(18)
c
       a(1)=1.0d0
       a(2)=.4227843350984678d0
       a(3)=.4118403304263672d0
      a(4)=.0815769192502609d0
      a(5)=.0742490106800904d0
      a(6)=-.0002669810333484d0
      a(7)=.0111540360240344d0
      a(8)=-.0028525821446197d0
      a(9)=.0021036287024598d0
      a(10)=-.0009184843690991d0
      a(11)=.0004874227944768d0
      a(12)=-.0002347204018919d0
      a(13)=.0001115339519666d0
      a(14)=-.0000478747983834d0
      a(15)=.0000175102727179d0
      a(16)=-.0000049203750904d0
      a(17)=.0000009199156407d0
      a(18)=-.0000000839940496d0
c
c
c
      if(z.le.1.0d0  ) go to 10
      if(z.le.2.0d0  ) go to 20
      t=z-2.0d0
      go to 30
10    t=z
      go to 30
20    t=z-1.0d0
30    p=a(18)
      do 40 k1=1,17
      k=18-k1
      p=t*p+a(k)
40    continue
c
      if(z.gt.2.0d0  ) go to 50
      if(z.gt.1.0d0  ) go to 60
      gamfun=p/(z*(z+1.0d0  ))
      return
60    gamfun=p/z
      return
50    gamfun=p
      return
      end
      subroutine plm(x,n,mu,lmax,p)
      parameter (nmx=500,ng7=7*nmx)
      common/stuff/scr(ng7)
      dimension p(nmx,0:1),x(nmx),p0(nmx),p1(nmx),p2(nmx)
      dimension arg(nmx),q0(nmx),q1(nmx),q2(nmx)
      equivalence (scr(1),arg(1)),
     1 (scr(1+1*nmx),p0(1)),(scr(1+2*nmx),p1(1)),
     1 (scr(1+3*nmx),q0(1)),(scr(1+4*nmx),q1(1)),
     1 (scr(1+5*nmx),p2(1)),(scr(1+6*nmx),q2(1))
c
c this routine calculates a vector of associated legengre polynomials P(l,mu) .
c for fixed mu and l=mu,,,lmax for an array of n points, x(i),i=1,,n.
c
c
c mu = 0 case
c
      if(mu.eq.0)then
       do 2 i=1,n
       p(i,0)=1.
2      p(i,1)=x(i)
       if(lmax.lt.2)return
       do 1 l=2,lmax
       do 3 i=1,n
3      p(i,l)=((2*l-1)*x(i)*p(i,l-1)-(l-1)*p(i,l-2))/l
1      continue
       return
      endif
c
c mu = 1 case
c
      if(mu.eq.1)then
       do 6 i=1,n
       arg(i)=sqrt(1.-x(i)*x(i))
6      p(i,1)=-arg(i)
       if(lmax.lt.2)return
       do 7 i=1,n
7      p(i,2)=-3.*x(i)*arg(i)
       if(lmax.lt.3)return
       do 8 l=3,lmax
       do 8 i=1,n
8      p(i,l)=((2*l-1)*x(i)*p(i,l-1)-l*p(i,l-2))/(l-1)
       return
      endif
c
c mu must be larger than 1
c
c
c recurr across to l=mu+1, with m=0
c
      mu1=mu+1
      do 4 i=1,n
      p0(i)=1.
4     p1(i)=x(i)
      do 5 l=2,mu1
      do 5 i=1,n
      p2(i)=((2*l-1)*x(i)*p1(i)-(l-1)*p0(i))/l
      p0(i)=p1(i)
      p1(i)=p2(i)
5     continue
c
c recurr across to l=mu+1, with m=1
c
      do 9 i=1,n
      arg(i)=sqrt(1.-x(i)*x(i))
      q0(i)=-arg(i)
9     q1(i)=-3.*arg(i)*x(i)
      do 10 l=3,mu1
      do 10 i=1,n
      q2(i)=((2*l-1)*x(i)*q1(i)-l*q0(i))/(l-1)
      q0(i)=q1(i)
10    q1(i)=q2(i)
c
c with l fixed at mu and mu+1, recurr down to m=mu
c
      l0=mu
      l1=mu1
      do 11 m=2,mu
      do 11 i=1,n
      p2(i)=-2.*(m-1)*x(i)*q0(i)/arg(i)-(l0-m+2)*(l0+m-1)*p0(i)
      q2(i)=-2.*(m-1)*x(i)*q1(i)/arg(i)-(l1-m+2)*(l1+m-1)*p1(i)
      p0(i)=q0(i)
      p1(i)=q1(i)
      q0(i)=p2(i)
11    q1(i)=q2(i)
c
c compute the desired vector and quit
c
      do 12 i=1,n
      p(i,mu)=q0(i)
12    continue
      if(lmax.eq.mu)return
      do 13 i=1,n
13    p(i,mu1)=q1(i)
      if(lmax.eq.mu1)return
      mu2=mu1+1
      do 14 l=mu2,lmax
      do 14 i=1,n
14    p(i,l)=((2*l-1)*x(i)*p(i,l-1)-(l+mu-1)*p(i,l-2))/(l-mu)
      return
      end
      subroutine sjymec(x,j,y,mn,nout)
c
c
c
c     jsubn(x),ysubn(x)      spherical bessel functions
c     do not use x=0.
c     jsub0(0)=1.   jsub n (0)=0.   for all n.gt.0
c     ysub n (0)= infinity  for all n.ge.0.
c     y sub n (x)=-j sub (-(n+1))  (x)        for all n.ge.0
c     mn  maximum number of subscript of bessel functions computed
c     nout   maximum subscript of desired bessel functions
c     recurrence techniques for bessel functions
c     fr.  mechel   math. of comp.  vol. 22   no. 101  jan. 1968
c
c
c
      real j,jx,jnmp,jnm,jnm1
      dimension j(0:1),y(0:1)
c
c     the dimension of j and y may be varied to suit the users needs
c     set dimension in calling routine  j((0,nout+5)),y((0,nout+5))
c
c
c
      data  bound /1.e13/, nadd /10/,  nxadd /10/
c
c     bound  nadd   nxadd   may be varied
c     nxadd controls the number of terms in the continued fraction
c
      mn=1
      j(0)=sin(x)/x
      y(0)=-cos(x)/x
      y(1)=y(0)/x-j(0)
      j(1)=j(0)/x+y(0)
      nmin=1
      nmax= nout
c
c     upward recursion valid for n<x for both bessel functions
c
      ix=x
      nx=min0(nout,ix)
      do 66 n=1,nx
      mn=n+1
      twon1=n+mn
      j(mn)=twon1*j(n)/x-j(n-1)
   66 y(mn)=twon1*y(n)/x-y(n-1)
      if(nout.le.mn) return
      jx=j(mn)
      mnout=mn
      nmin=mn
c
c     upward recursion for ysubn(x) to obtain good guess for jsubn(x)
c
      do 1 n=nmin,nout
      mn=n+1
      twon1=n+mn
    1 y(mn)=twon1*y(n)/x-y(n-1)
      ynm1=y(nout-1)
      ynm=y(nout)
      nmin=nout
      nmax=nout+nadd
   33 do 11 n=nmin,nmax
      mn=n+1
      twon1=n+mn
      ynmp=twon1*ynm/x-ynm1
      ynm1=ynm
      ynm=ynmp
   11 continue
      if(abs(ynmp).gt.bound) 2,22
   22 nmin=nmax+1
      nmax=nmax+nadd
      go to 33
    2 x2=x*x
      jnm=-1.0/(x2*ynmp)
      pmin1=1.
      qmin1=0.
      q0=1.0
c
c     develop continued fraction to obtain good guess for jsubn+1(x)
c
      b0=mn+mn+1
      p0=mn+mn+1
      do 4 n=1,nxadd
      b0=b0+2.
      p1=b0*p0-x2*pmin1
      q1=b0*q0-x2*qmin1
      pmin1=p0
      p0=p1
      qmin1=q0
      q0=q1
    4 continue
c
      gnmax1=p1/q1
      jnmp=x*jnm/gnmax1
      do 55 n=mn-1,nout+1,-1
      nm=n+1
      twon1=n+nm
      jnm1=twon1*jnm/x-jnmp
      jnmp=jnm
      jnm=jnm1
   55 continue
      j(nout+1)=jnmp
      j(nout)=jnm1
      do 5 n=nout,mnout+1,-1
      nm=n+1
      twon1=n+nm
    5 j(n-1)=twon1*j(n)/x-j(n+1)
      c=jx/j(mnout)
c
c     correct bessel functions with weighting factor
c
      if (x.lt.1.e-3) return
      do 6 n=mnout,nout
    6 j(n)=c*j(n)
c
      return
      end
