*deck @(#)init4x.f	5.1  11/6/94
      subroutine init4x(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,
     $                  ss)
      implicit real*8 (a-h,o-z)
      integer xor
      integer arr,tr1,tr2,asm,aos,os,wtw,wtx,wty,wab,ss,ssi,ssj,symorb
      integer bmax,orbfrm
      real*8 int(nmax),c(nwks),s(nwks),int1(2)
      dimension kadd(symorb),ladd(symorb),ijadd(numij),wtw(orbfrm,nsym)
      dimension wtx(orbfrm,nsym),wty(orbfrm),wab(orbfrm),ss(norbs)
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,nrow4o,nwksoc,nlevoc,norboc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
c  universal identity of the objects in these common
c     common /all/ val1,val2,val3,arr,itr1,itr2,ia,ja,itype,isegt
c    *,lvfrm1,nlwki,nlwkj,imax,imin
c     common /symm/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
c    #,             numsym(8)
      common /symm/ asm,js,is,mx(8),mn(8),os(8),numsym(8)
      common /tapes/itap20,itape5,itape6,itape8,itap12,itap03,itap04
     *,             itape3,itap05,itap06
      common /all/val1,val2,val3,arr,tr1,tr2,ia,ja,m,iseg,n,n1,n2
     *,           imax,imin
      common /minmax/ iming,imaxg,jming,jmaxg
      common/dry1/idry,jdry,ksegsv,indxs(100)
c
      real*8 z,val1,val2,val3,zl,zr
      sqrt2=sqrt(2.0d+00)
      sqrt3=sqrt(3.0d+00)
      sqt1p5=sqrt(1.5d+00)
      return
c
c
c
      entry shape4(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,
     $                  ss)
c
      aos=os(asm+1)
      go to (301,301,301,301,301,301,301,301,301,301,301,301,301,314
c     go to ( yz, xz, xy, yy, xx, zy, yx, yw, wz, wy, wx, xw, ww,---yy
c
     1,315,316,317,318,319),m
c       xx, ww, zz,---4 internals),m
c
      stop 'funny m'
c
  301 print 302
  302 format(' this entry is not possible')
      stop
c
c    yy    no  9    four externals
c***********************************************************************
  314 if (is.gt.js) return
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      m3=mx(is+1)
      m4=mn(is+1)
      if (m4.gt.n) return
      m5=imax
      m6=imin
      if (m5.gt.m1) m5=m1
      if (m6.lt.m2) m6=m2
      if (m5.lt.m6) return
      jj=ja+wty(m6)
      do 2000 i=m6,m5
      if (i.gt.jmaxg.or.i.lt.jming) go to 2000
      if (is.eq.js) m3=i-1
      if (m3.lt.m4) go to 2000
      kos=os(1)
      los=os(ss(i)+1)
      lad=ijadd(i*(i+1)/2)+kadd(i+kos)+3
      ii=ia+wty(m4)
      cjj=c(jj)
      sjj=s(jj)
      do 2001 j=m4,m3
      lad1=lad+ladd(j+los)
c     z=int(lad1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad1)=int(lad1)+(cjj*s(ii)+c(ii)*sjj)
c      d(jj,ii;i,j)
      iix=indxs(i)+j
      int1(iix)=int1(iix)+cjj*s(ii)
      iix=indxs(j)+i
      int1(iix)=int1(iix)+c(ii)*sjj
c
 2001 ii=ii+1
c     s(jj)=s(jj)+sjj
 2000 jj=jj+1
      return
c
c   xx   no 13 & 161 four externals
c***********************************************************************
  315 m5=imax
      m6=imin
      if (m6.lt.3) m6=3
      if (m5.lt.m6) go to 2003
      do 2004 i=m6,m5
      if (i.lt.jming.or.i.gt.jmaxg) go to 2004
      ssj=xor(ss(i),js)
      ssi=xor(ss(i),is)
      if (ssj.gt.ss(i)) go to 2004
      if (ssi.gt.ssj) go to 2004
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 2004
      if (ssj.eq.ss(i)) m1=i-1
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 2004
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 2004
      jj=ja+wtx(i,js+1)+wty(m2)
      iia=ia+wtx(i,is+1)
      kos=os(1)
      lad=ijadd(i*(i+1)/2)
      do 2005 j=m2,m1
      if (ssi.eq.ssj) m3=j-1
      if (m3.lt.m4) go to 2005
      los=os(ss(j)+1)
      lad1=lad+kadd(j+kos)
      ii=iia+wty(m4)
      cjj=c(jj)
      sjj=s(jj)
      do 2006 k=m4,m3
      lad2=lad1+ladd(k+los)
c     z=-int(lad2+1)+int(lad2+2)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=(cjj*s(ii)+c(ii)*sjj)
      int(lad2+1)=int(lad2+1)-z
      int(lad2+2)=int(lad2+2)+z
c
 2006 ii=ii+1
c     s(jj)=s(jj)+sjj
 2005 jj=jj+1
 2004 continue
 2003 continue
c
c     no 161
      if (n.lt.3) go to 8084
      do 7084 i=3,n
      ssj=xor(ss(i),js)
      ssi=xor(ss(i),is)
      if (ssj.gt.ss(i)) go to 7084
      if (ssi.gt.ssj) go to 7084
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 7084
      if (m2.lt.2) m2=2
      if (ssj.eq.ss(i)) m1=i-1
      if (m1.lt.m2) go to 7084
      m5=imax
      m6=imin
      if (m5.gt.m1) m5=m1
      if (m6.lt.m2) m6=m2
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 7084
      jj=ja+wtx(i,js+1)+wty(m6)
      iia=ia+wtx(i,is+1)
      kos=os(1)
      do 7085 j=m6,m5
      if (j.lt.jming.or.j.gt.jmaxg) go to 7085
      if (ssi.eq.ssj) m3=j-1
      if (m3.lt.m4) go to 7085
      los=os(ss(j)+1)
      lad3=ijadd(j*(j+1)/2)+kadd(j+kos)+3
      ii=iia+wty(m4)
      cjj=c(jj)
      sjj=s(jj)
      do 7086 k=m4,m3
c     z=int(lad3+ladd(k+los))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad3+ladd(k+los))=int(lad3+ladd(k+los))+(cjj*s(ii)+c(ii)*sjj)
c
c      d(jj,ii;j,k)
      iix=indxs(j)+k
      int1(iix)=int1(iix)+cjj*s(ii)
      iix=indxs(k)+j
      int1(iix)=int1(iix)+c(ii)*sjj
 7086 ii=ii+1
c     s(jj)=s(jj)+sjj
 7085 jj=jj+1
 7084 continue
 8084 continue
c   no  14a
c
      m5=imax
      m6=imin
      if (m6.lt.3) m6=3
      if (m5.lt.m6) go to 8007
      do 2007 i=m6,m5
      ssj=xor(ss(i),js)
      ssi=xor(is,ssj)
      if (ssj.gt.ssi) go to 2007
      if (ssi.gt.ss(i)) go to 2007
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 2007
      if (m2.lt.2) m2=2
      if (ssi.eq.ss(i)) m1=i-1
      if (m1.lt.m2) go to 2007
      m3=mx(ssj+1)
      m4=mn(ssj+1)
      if (m4.gt.n) go to 2007
      jja=ja+wtx(i,js+1)
      kos=os(1)
      los1=os(ss(i)+1)
c     lad3=ijadd(i*(i+1)/2)+kadd(i+kos)+3
      do 2008 j=m2,m1
      if (j.lt.jming.or.j.gt.jmaxg) go to 2008
      if (ssi.eq.ssj) m3=j-1
      if (m3.lt.m4) go to 2008
c     z1=int(lad3+ladd(j+los1))
c     z1=0.0d+00
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      lad1=ijadd(i*(i-1)/2+j)
      ii=ia+wtx(j,is+1)+wty(m4)
      jj=jja+wty(m4)
      do 2009 k=m4,m3
      lsm=xor(ksm,ss(k))
      los=os(lsm+1)
      lad2=lad1+kadd(k+kos)+ladd(k+los)
c     z=-int(lad2+1)+int(lad2+2)+z1
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=(c(ii)*s(jj)+c(jj)*s(ii))
c     z1=z1+z
      int(lad2+1)=int(lad2+1)-z
      int(lad2+2)=int(lad2+2)+z
c
      ii=ii+1
 2009 jj=jj+1
c     int(lad3+ladd(j+los1))=int(lad3+ladd(j+los1))+z1
c
 2008 continue
 2007 continue
 8007 continue
c   no  14b
c
      m5=imax
      m6=imin
      if (m6.lt.3) m6=3
      if (m5.lt.m6) go to 4007
      do 3007 i=m6,m5
      if (i.lt.jming.or.i.gt.jmaxg) go to 3007
      ssj=xor(ss(i),js)
      ssi=xor(is,ssj)
      if (ssj.gt.ssi) go to 3007
      if (ssi.gt.ss(i)) go to 3007
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 3007
      if (m2.lt.2) m2=2
      if (ssi.eq.ss(i)) m1=i-1
      if (m1.lt.m2) go to 3007
      m3=mx(ssj+1)
      m4=mn(ssj+1)
      if (m4.gt.n) go to 3007
      jja=ja+wtx(i,js+1)
      kos=os(1)
      los1=os(ss(i)+1)
      lad3=ijadd(i*(i+1)/2)+kadd(i+kos)+3
      do 3008 j=m2,m1
      if (ssi.eq.ssj) m3=j-1
      if (m3.lt.m4) go to 3008
c     z1=int(lad3+ladd(j+los1))
      z1=0.0d+00
      zr=z1
      zl=z1
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      lad1=ijadd(i*(i-1)/2+j)
      ii=ia+wtx(j,is+1)+wty(m4)
      jj=jja+wty(m4)
      do 3009 k=m4,m3
      lsm=xor(ksm,ss(k))
      los=os(lsm+1)
c     lad2=lad1+kadd(k+kos)+ladd(k+los)
c     z=-int(lad2+1)+int(lad2+2)+z1
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z1=z1+(c(ii)*s(jj)+c(jj)*s(ii))
      zl=zl+c(ii)*s(jj)
      zr=zr+c(jj)*s(ii)
c     int(lad2+1)=int(lad2+1)-z
c     int(lad2+2)=int(lad2+2)+z
c
      ii=ii+1
 3009 jj=jj+1
      int(lad3+ladd(j+los1))=int(lad3+ladd(j+los1))+z1
c
c      d(jj,ii;i,j)
      iix=indxs(i)+j
      int1(iix)=int1(iix)+zr
      iix=indxs(j)+i
      int1(iix)=int1(iix)+zl
 3008 continue
 3007 continue
 4007 continue
c
c    no 16   both
c
      if (m5.lt.m6) go to 8010
      do 2010 i=m6,m5
      if (i.lt.jming.or.i.gt.jmaxg) go to 2010
      ssj=xor(ss(i),js)
      ssi=xor(is,ssj)
      if (ssi.gt.ssj) go to 2010
      if (ssj.gt.ss(i)) go to 2010
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 2010
      if (ssj.eq.ss(i)) m1=i-1
      if (m2.lt.2) m2=2
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 2010
      jj=ja+wtx(i,js+1)+wty(m2)
      kos=os(1)
      los=os(ss(i)+1)
      lad2=ijadd(i*(i+1)/2)+kadd(i+kos)
      do 2011 j=m2,m1
      if (ssi.eq.ssj) m3=j-1
      if (m3.lt.m4) go to 2011
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      lad3=ijadd(i*(i-1)/2+j)+kadd(j+kos)
      ii=ia+wtx(j,is+1)+wty(m4)
      cjj=c(jj)
      sjj=s(jj)
      do 2012 k=m4,m3
      lad4=lad2+ladd(k+los)
      lad5=lad3+ladd(k+los)
c     z=-int(lad4+3)-int(lad5+2)+int(lad5+1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=(cjj*s(ii)+c(ii)*sjj)
      int(lad4+3)=int(lad4+3)-z
      int(lad5+2)=int(lad5+2)-z
      int(lad5+1)=int(lad5+1)+z
c
c      d(jj,ii;i,k)
      iix=indxs(i)+k
      int1(iix)=int1(iix)-cjj*s(ii)
      iix=indxs(k)+i
      int1(iix)=int1(iix)-c(ii)*sjj
 2012 ii=ii+1
c     s(jj)=s(jj)+sjj
 2011 jj=jj+1
 2010 continue
 8010 continue
c
c    no  20
c
      if (m5.lt.4) go to 8013
      m6=imin
      if (m6.lt.4) m6=4
      do 2013 i=m6,m5
      ssj=xor(ss(i),js)
      if (ssj.gt.ss(i)) go to 2013
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 2013
      if (ssj.eq.ss(i)) m1=i-1
      if (m2.lt.3) m2=3
      if (m1.lt.m2) go to 2013
      jj=ja+wtx(i,js+1)+wty(m2)
      do 2014 j=m2,m1
      if (j.lt.jming.or.j.gt.jmaxg) go to 2014
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      lad1=ijadd(i*(i-1)/2+j)
      cjj=c(jj)
      sjj=s(jj)
      do 2015 k=2,j-1
      ssi=xor(ss(k),is)
      if (ssi.gt.ss(k)) go to 2015
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 2015
      if (ss(k).eq.ssi) m3=k-1
      if (m3.lt.m4) go to 2015
      ii=ia+wtx(k,is+1)+wty(m4)
      lsm=xor(ksm,ss(k))
      los=os(lsm+1)
      lad2=lad1+kadd(k+kos)
      do 2016 l=m4,m3
      lad3=lad2+ladd(l+los)
c     z=int(lad3+1)-int(lad3+3)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=(cjj*s(ii)+c(ii)*sjj)
      int(lad3+1)=int(lad3+1)+z
      int(lad3+3)=int(lad3+3)-z
c
 2016 ii=ii+1
 2015 continue
c     s(jj)=s(jj)+sjj
 2014 jj=jj+1
 2013 continue
 8013 continue
c
c      no  21
c
      if(m5.lt.m6) go to 8017
      do 2017 i=m6,m5
      ssj=xor(ss(i),js)
      if (ssj.gt.ss(i)) go to 2017
      jja=ja+wtx(i,js+1)
      do 2018 j=3,i-1
      if (j.lt.jming.or.j.gt.jmaxg) go to 2018
      ssi=xor(ss(j),is)
      if (ssi.gt.ss(j)) go to 2018
      if (ssj.gt.ssi) go to 2018
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 2018
      if (ssi.eq.ss(j)) m1=j-1
      if (m2.lt.2) m2=2
      if (m1.lt.m2) goto 2018
      ii=ia+wtx(j,is+1)+wty(m2)
      m3=mx(ssj+1)
      m4=mn(ssj+1)
      if (m4.gt.n) go to 2018
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      lad1=ijadd(i*(i-1)/2+j)
      do 2019 k=m2,m1
      if (ssj.eq.ssi) m3=k-1
      if (m3.lt.m4) go to 2019
      jj=jja+wty(m4)
      lad2=lad1+kadd(k+kos)
      lsm=xor(ksm,ss(k))
      los=os(lsm+1)
      cii=c(ii)
      sii=s(ii)
      do 2020 l=m4,m3
      lad3=lad2+ladd(l+los)
c     z=-int(lad3+1)+int(lad3+2)
c     sii=sii+z*c(jj)
c     s(jj)=s(jj)+z*cii
      z=(cii*s(jj)+c(jj)*sii)
      int(lad3+1)=int(lad3+1)-z
      int(lad3+2)=int(lad3+2)+z
c
 2020 jj=jj+1
c     s(ii)=s(ii)+sii
 2019 ii=ii+1
 2018 continue
 2017 continue
 8017 continue
c
c   no  22
c
      if(m5.lt.m6) return
      do 2021 i=m6,m5
      ssj=xor(ss(i),js)
      if (ssj.gt.ss(i)) go to 2021
      jja=ja+wtx(i,js+1)
      do 2022 j=3,i-1
      if (j.lt.jming.or.j.gt.jmaxg) go to 2022
      ssi=xor(ss(j),is)
      if (ssj.gt.ss(j))go to 2022
      if (ssi.gt.ssj) go to 2022
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 2022
      if (ssj.eq.ss(j)) m1=j-1
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 2022
      jj=jja+wty(m2)
      iia=ia+wtx(j,is+1)
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 2022
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      lad1=ijadd(i*(i-1)/2+j)
      do 2023 k=m2,m1
      if (ssi.eq.ssj) m3=k-1
      if (m3.lt.m4) go to 2023
      ii=iia+wty(m4)
      lsm=xor(ksm,ss(k))
      los=os(lsm+1)
      lad2=lad1+kadd(k+kos)
      cjj=c(jj)
      sjj=s(jj)
      do 2024 l=m4,m3
      lad3=lad2+ladd(l+los)
c     z=-int(lad3+3)+int(lad3+2)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=(cjj*s(ii)+c(ii)*sjj)
      int(lad3+3)=int(lad3+3)-z
      int(lad3+2)=int(lad3+2)+z
c
 2024 ii=ii+1
c     s(jj)=s(jj)+sjj
 2023 jj=jj+1
 2022 continue
 2021 continue
      return
c
c   ww   no 69a four externals
c***********************************************************************
  316 m5=imax
      m6=imin
      if (m6.lt.2) m6=2
      if (m5.lt.m6) go to 6586
      if (is.ne.0.or.js.ne.0) go to 6586
      do 2031 i=m6,m5
      if (i.lt.jming.or.i.gt.jmaxg) go to 2031
      jj=ja+wab(i)
      lad=ijadd(i*(i+1)/2)
      cjj=c(jj)
      kos=os(1)
      sjj=s(jj)
      do 2032 j=1,i-1
      los=os(ss(j)+1)
      lad1=kadd(j+kos)+ladd(j+los)+lad
      ii=ia+wab(j)
c     z=int(lad1+1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad1+1)=int(lad1+1)+(cjj*s(ii)+c(ii)*sjj)
c
 2032 continue
c     s(jj)=s(jj)+sjj
 2031 continue
c
c    no   72a
c
 6586 if (js.ne.0) go to  6635
      if (m5.lt.m6) go to 6635
      do 2033 i=m6,m5
      if (i.lt.jming.or.i.gt.jmaxg) go to 2033
      ssi=xor(ss(i),is)
      if (ssi.gt.ss(i)) go to 2033
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 2033
      if (ssi.eq.ss(i)) m1=i-1
      if (m1.lt.m2) goto 2033
      jj=ja+wab(i)
      ii=ia+wtw(i,is+1)+wty(m2)
      kos=os(1)
      los=os(ss(i)+1)
      lad=ijadd(i*(i+1)/2)+kadd(i+kos)
      cjj=c(jj)
      sjj=s(jj)
      do 2034 j=m2,m1
      lad1=lad+ladd(j+los)
c     z=sqrt2*(int(lad1+3)+int(lad1+1))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=(cjj*s(ii)+c(ii)*sjj)*sqrt2
      int(lad1+3)=int(lad1+3)+z
      int(lad1+1)=int(lad1+1)+z
c      d(jj,ii;i,j)
      iix=indxs(i)+j
      int1(iix)=int1(iix)+sqrt2*cjj*s(ii)
      iix=indxs(j)+i
      int1(iix)=int1(iix)+sqrt2*c(ii)*sjj
c
 2034 ii=ii+1
c     s(jj)=s(jj)+sjj
 2033 continue
c
c    no  73
c
 6635 m6=imin
      if (m6.lt.3) m6=3
      if (m5.lt.m6) go to 6636
      if (is.ne.0) go to 6636
      do 2035 i=m6,m5
      ssj=xor(ss(i),js)
      if (ssj.gt.ss(i)) go to 2035
      jja=ja+wtw(i,js+1)
      m3=mx(ssj+1)
      m4=mn(ssj+1)
      if (m4.gt.n) go to 2035
      do 2036 j=2,i-1
      if (j.lt.jming.or.j.gt.jmaxg) go to 2036
      if (ssj.gt.ss(j)) go to 2036
      if (ssj.eq.ss(j)) m3=j-1
      if (m3.lt.m4) go to 2036
      jj=jja+wty(m4)
      ii=ia+wab(j)
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      los=os(ss(i)+1)
      lad1=ijadd(i*(i-1)/2+j)+kadd(j+kos)
      cii=c(ii)
      sii=s(ii)
      do 2037 k=m4,m3
      lad2=lad1+ladd(k+los)
c     z=sqrt2*int(lad2+1)
c     sii=sii+z*c(jj)
c     s(jj)=s(jj)+z*cii
      int(lad2+1)=int(lad2+1)+sqrt2*(cii*s(jj)+c(jj)*sii)
c
 2037 jj=jj+1
c     s(ii)=s(ii)+sii
 2036 continue
 2035 continue
c
c   no  74a
c
 6636 if (js.ne.0) go to 6637
      if (m5.lt.m6) go to 6637
      do 2038 i=m6,m5
      if (i.lt.jming.or.i.gt.jmaxg) go to 2038
      jj=ja+wab(i)
      lad=ijadd(i*(i+1)/2)
      cjj=c(jj)
      kos=os(1)
      sjj=s(jj)
      do 2039 j=2,i-1
      ssi=xor(ss(j),is)
      if (ssi.gt.ss(j)) go to 2039
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 2039
      if (ssi.eq.ss(j)) m1=j-1
      if (m1.lt.m2) go to 2039
      ii=ia+wtw(j,is+1)+wty(m2)
      los=os(ss(j)+1)
      lad1=lad+kadd(j+kos)
      do 2040 k=m2,m1
      lad2=lad1+ladd(k+los)
c     z=sqrt2*int(lad2+1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad2+1)=int(lad2+1)+sqrt2*(cjj*s(ii)+c(ii)*sjj)
c
 2040 ii=ii+1
 2039 continue
c     s(jj)=s(jj)+sjj
 2038 continue
c
c   no  76
c
 6637 if (is.ne.0) go to 6638
      if (m5.lt.m6) go to 8041
      do 2041 i=m6,m5
      ssj=xor(ss(i),js)
      if (ssj.gt.ss(i)) go to 2041
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 2041
      if (ssj.eq.ss(i)) m1=i-1
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 2041
      jj=ja+wtw(i,js+1)+wty(m2)
      do 2042 j=m2,m1
      if (j.lt.jming.or.j.gt.jmaxg) go to 2042
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      lad1=ijadd(i*(i-1)/2+j)
      cjj=c(jj)
      sjj=s(jj)
      do 2043 k=1,j-1
      ii=ia+wab(k)
      lsm=xor(ksm,ss(k))
      los=os(lsm+1)
      lad2=lad1+kadd(k+kos)+ladd(k+los)
c     z=sqrt2*int(lad2+1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad2+1)=int(lad2+1)+sqrt2*(cjj*s(ii)+c(ii)*sjj)
c
 2043 continue
c     s(jj)=s(jj)+sjj
 2042 jj=jj+1
c2042 continue
 2041 continue
 8041 continue
c
c    no  78
c
      m6=imin
      if (m6.lt.2) m6=2
      if (m5.lt.m6) go to 6638
      do 2044 i=m6,m5
      if (i.lt.jming.or.i.gt.jmaxg) go to 2044
      ssj=xor(ss(i),js)
      if (ssj.gt.ss(i)) go to 2044
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 2044
      if (ssj.eq.ss(i)) m1=i-1
      if (m1.lt.m2) go to 2044
      jj=ja+wtw(i,js+1)+wty(m2)
      kos=os(1)
      los=os(ss(i)+1)
      lad=ijadd(i*(i+1)/2)+kadd(i+kos)
      do 2045 j=m2,m1
      lad1=lad+ladd(j+los)
      ii=ia+wab(j)
c     z=sqrt2*(int(lad1+3)+int(lad1+2))
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=sqrt2*(c(ii)*s(jj)+c(jj)*s(ii))
      int(lad1+3)=int(lad1+3)+z
      int(lad1+2)=int(lad1+2)+z
c      d(jj,ii;i,j)
      iix=indxs(i)+j
      int1(iix)=int1(iix)+sqrt2*c(jj)*s(ii)
      iix=indxs(j)+i
      int1(iix)=int1(iix)+sqrt2*c(ii)*s(jj)
c
 2045 jj=jj+1
 2044 continue
c
c    no   84a'
c
 6638 m6=imin
      if (m6.lt.3) m6=3
      if (m5.lt.m6) goto 8054
      do 3048 i=m6,m5
      ssj=xor(ss(i),js)
      ssi=xor(is,ssj)
      if (ssj.gt.ssi) go to 3048
      if (ssi.gt.ss(i)) go to 3048
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 3048
      m3=mx(ssj+1)
      m4=mn(ssj+1)
      if (m4.gt.n) go to 3048
      if (ssi.eq.ss(i)) m1=i-1
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 3048
      jja=ja+wtw(i,js+1)
      kos1=os(1)
      los1=os(ss(i)+1)
c     lad3=ijadd(i*(i+1)/2)+kadd(i+kos1)+3
      do 3049 j=m2,m1
      if (j.lt.jming.or.j.gt.jmaxg) go to 3049
      if (ssi.eq.ssj) m3=j-1
      if (m3.lt.m4) go to 3049
c     z1=int(lad3+ladd(j+los1))
c     z1=0.0d+00
      ii=ia+wtw(j,is+1)+wty(m4)
      jj=jja+wty(m4)
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      lad1=ijadd(i*(i-1)/2+j)
      do 3050 k=m4,m3
      lsm=xor(ksm,ss(k))
      los=os(lsm+1)
      lad2=lad1+kadd(k+kos)+ladd(k+los)
c     z=int(lad2+1)+int(lad2+2)+z1
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=(c(ii)*s(jj)+c(jj)*s(ii))
c     z1=z1+z
      int(lad2+1)=int(lad2+1)+z
      int(lad2+2)=int(lad2+2)+z
c
      ii=ii+1
 3050 jj=jj+1
c
c     int(lad3+ladd(j+los1))=int(lad3+ladd(j+los1))+z1
c
 3049 continue
 3048 continue
c
c    no   84a''
c
      m6=imin
      if (m6.lt.3) m6=3
      if (m5.lt.m6) goto 8054
      do 2048 i=m6,m5
      if (i.lt.jming.or.i.gt.jmaxg) go to 2048
      ssj=xor(ss(i),js)
      ssi=xor(is,ssj)
      if (ssj.gt.ssi) go to 2048
      if (ssi.gt.ss(i)) go to 2048
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 2048
      m3=mx(ssj+1)
      m4=mn(ssj+1)
      if (m4.gt.n) go to 2048
      if (ssi.eq.ss(i)) m1=i-1
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 2048
      jja=ja+wtw(i,js+1)
      kos1=os(1)
      los1=os(ss(i)+1)
      lad3=ijadd(i*(i+1)/2)+kadd(i+kos1)+3
      do 2049 j=m2,m1
      if (ssi.eq.ssj) m3=j-1
      if (m3.lt.m4) go to 2049
c     z1=int(lad3+ladd(j+los1))
      z1=0.0d+00
      zl=z1
      zr=z1
      ii=ia+wtw(j,is+1)+wty(m4)
      jj=jja+wty(m4)
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
c     lad1=ijadd(i*(i-1)/2+j)
      do 2050 k=m4,m3
      lsm=xor(ksm,ss(k))
      los=os(lsm+1)
c     lad2=lad1+kadd(k+kos)+ladd(k+los)
c     z=int(lad2+1)+int(lad2+2)+z1
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z1=z1+(c(ii)*s(jj)+c(jj)*s(ii))
      zl=zl+c(jj)*s(ii)
      zr=zr+c(ii)*s(jj)
c     int(lad2+1)=int(lad2+1)+z
c     int(lad2+2)=int(lad2+2)+z
c
      ii=ii+1
 2050 jj=jj+1
c
      int(lad3+ladd(j+los1))=int(lad3+ladd(j+los1))+z1
c      d(jj,ii;i,j)
      iix=indxs(i)+j
      int1(iix)=int1(iix)+zl
      iix=indxs(j)+i
      int1(iix)=int1(iix)+zr
c
 2049 continue
 2048 continue
c
c   no  88a'
c
      do 3051 i=m6,m5
      ssj=xor(ss(i),js)
      ssi=xor(is,ssj)
      if (ssi.gt.ssj) go to 3051
      if (ssj.gt.ss(i)) go to 3051
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 3051
      if (ssj.eq.ss(i)) m1=i-1
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 3051
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 3051
      jj=ja+wtw(i,js+1)+wty(m2)
      kos1=os(1)
      los1=os(ss(i)+1)
c     lad1=ijadd(i*(i+1)/2)+kadd(i+kos1)
      do 3052 j=m2,m1
      if (j.lt.jming.or.j.gt.jmaxg) go to 3052
      if (ssi.eq.ssj) m3=j-1
      if (m3.lt.m4) go to 3052
      ii=ia+wtw(j,is+1)+wty(m4)
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      lad2=ijadd(i*(i-1)/2+j)+kadd(j+kos)
      cjj=c(jj)
      sjj=s(jj)
      do 3053 k=m4,m3
c     z=int(lad1+ladd(k+los1)+3)+int(lad2+ladd(k+los1)+2)
c     z=z+int(lad2+ladd(k+los1)+1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=(cjj*s(ii)+c(ii)*sjj)
c     int(lad1+ladd(k+los1)+3)=int(lad1+ladd(k+los1)+3)+z
      lad2p=lad2+ladd(k+los1)
      int(lad2p+1)=int(lad2p+1)+z
      int(lad2p+2)=int(lad2p+2)+z
c
 3053 ii=ii+1
c     s(jj)=s(jj)+sjj
 3052 jj=jj+1
 3051 continue
c
c   no  88a''
c
      do 2051 i=m6,m5
      if (i.lt.jming.or.i.gt.jmaxg) go to 2051
      ssj=xor(ss(i),js)
      ssi=xor(is,ssj)
      if (ssi.gt.ssj) go to 2051
      if (ssj.gt.ss(i)) go to 2051
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 2051
      if (ssj.eq.ss(i)) m1=i-1
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 2051
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 2051
      jj=ja+wtw(i,js+1)+wty(m2)
      kos1=os(1)
      los1=os(ss(i)+1)
      lad1=ijadd(i*(i+1)/2)+kadd(i+kos1)
      do 2052 j=m2,m1
      if (ssi.eq.ssj) m3=j-1
      if (m3.lt.m4) go to 2052
      ii=ia+wtw(j,is+1)+wty(m4)
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
c     lad2=ijadd(i*(i-1)/2+j)+kadd(j+kos)
      cjj=c(jj)
      sjj=s(jj)
      do 2053 k=m4,m3
c     z=int(lad1+ladd(k+los1)+3)+int(lad2+ladd(k+los1)+2)
c     z=z+int(lad2+ladd(k+los1)+1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=(cjj*s(ii)+c(ii)*sjj)
      int(lad1+ladd(k+los1)+3)=int(lad1+ladd(k+los1)+3)+z
c     lad2p=lad2+ladd(k+los1)
c     int(lad2p+1)=int(lad2p+1)+z
c     int(lad2p+2)=int(lad2p+2)+z
c      d(jj,ii;i,k)
      iix=indxs(i)+k
      int1(iix)=int1(iix)+cjj*s(ii)
      iix=indxs(k)+i
      int1(iix)=int1(iix)+c(ii)*sjj
c
 2053 ii=ii+1
c     s(jj)=s(jj)+sjj
 2052 jj=jj+1
 2051 continue
c
c   no   89 & 160
c
      do 2054 i=m6,m5
      if (i.lt.jming.or.i.gt.jmaxg) go to 2054
      ssj=xor(ss(i),js)
      ssi=xor(ss(i),is)
      if (ssj.gt.ss(i)) go to 2054
      if (ssi.gt.ssj) go to 2054
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 2054
      if (ssj.eq.ss(i)) m1=i-1
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 2054
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 2054
      jj=ja+wtw(i,js+1)+wty(m2)
      iia=ia+wtw(i,is+1)
      kos=os(1)
      lad=ijadd(i*(i+1)/2)
      do 2055 j=m2,m1
      if (ssi.eq.ssj) m3=j-1
      if (m3.lt.m4) go to 2055
      ii=iia+wty(m4)
      lad1=lad+kadd(j+kos)
      los=os(ss(j)+1)
      cjj=c(jj)
      sjj=s(jj)
      do 2056 k=m4,m3
      lad2=lad1+ladd(k+los)
c     z=int(lad2+1)+int(lad2+2)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=(cjj*s(ii)+c(ii)*sjj)
      int(lad2+1)=int(lad2+1)+z
      int(lad2+2)=int(lad2+2)+z
c
 2056 ii=ii+1
c     s(jj)=s(jj)+sjj
 2055 jj=jj+1
 2054 continue
 8054 continue
c
c  no 160
      if (n.lt.3) go to 7080
      do 7087 i=3,n
      ssj=xor(ss(i),js)
      ssi=xor(ss(i),is)
      if (ssj.gt.ss(i)) go to 7087
      if (ssi.gt.ssj) go to 7087
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 7087
      if (ssj.eq.ss(i)) m1=i-1
      if (m2.lt.2) m2=2
      m5=imax
      m6=imin
      if (m5.gt.m1) m5=m1
      if (m6.lt.m2) m6=m2
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 7087
      jj=ja+wtw(i,js+1)+wty(m6)
      iia=ia+wtw(i,is+1)
      kos=os(1)
      do 7088 j=m6,m5
      if (j.lt.jming.or.j.gt.jmaxg) go to 7088
      los=os(ss(j)+1)
      lad3=ijadd(j*(j+1)/2)+kadd(j+kos)+3
      if (ssi.eq.ssj) m3=j-1
      if (m3.lt.m4) go to 7088
      ii=iia+wty(m4)
      cjj=c(jj)
      sjj=s(jj)
      do 7089 k=m4,m3
c     z=int(lad3+ladd(k+los))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad3+ladd(k+los))=int(lad3+ladd(k+los))+(cjj*s(ii)+c(ii)*sjj)
c      d(jj,ii;j,k)
      iix=indxs(j)+k
      int1(iix)=int1(iix)+cjj*s(ii)
      iix=indxs(k)+j
      int1(iix)=int1(iix)+c(ii)*sjj
c
 7089 ii=ii+1
c     s(jj)=s(jj)+sjj
 7088 jj=jj+1
 7087 continue
 7080 continue
c     no    99a
c
      m5=imax
      m6=imin
      if (m6.lt.4) m6=4
      if (m5.lt.m6) return
      do 2057 i=m6,m5
      ssj=xor(ss(i),js)
      if (ssj.gt.ss(i)) go to 2057
      jja=ja+wtw(i,js+1)
      do 2058 j=3,i-1
      if (j.lt.jming.or.j.gt.jmaxg) go to 2058
      ssi=xor(ss(j),is)
      if(ssj.gt.ss(j)) go to 2058
      if (ssi.gt.ssj) go to 2058
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 2058
      if (ssj.eq.ss(j)) m1=j-1
      if (m2.lt.2) m2=2
      jj=jja+wty(m2)
      iia=ia+wtw(j,is+1)
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 2058
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      lad1=ijadd(i*(i-1)/2+j)
      do 2059 k=m2,m1
      if (ssi.eq.ssj) m3=k-1
      if (m3.lt.m4) go to 2059
      ii=iia+wty(m4)
      lsm=xor(ksm,ss(k))
      los=os(lsm+1)
      lad2=lad1+kadd(k+kos)
      cjj=c(jj)
      sjj=s(jj)
      do 2060 l=m4,m3
      lad3=lad2+ladd(l+los)
c     z=int(lad3+3)+int(lad3+2)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=(cjj*s(ii)+c(ii)*sjj)
      int(lad3+3)=int(lad3+3)+z
      int(lad3+2)=int(lad3+2)+z
c
 2060 ii=ii+1
c     s(jj)=s(jj)+sjj
 2059 jj=jj+1
 2058 continue
 2057 continue
c
c     no   100a
c
      do 2061 i=m6,m5
      ssj=xor(ss(i),js)
      if (ssj.gt.ss(i)) go to 2061
      jja=ja+wtw(i,js+1)
      do 2062 j=3,i-1
      if (j.lt.jming.or.j.gt.jmaxg) go to 2062
      ssi=xor(ss(j),is)
      if (ssi.gt.ss(j)) go to 2062
      if (ssj.gt.ssi) go to 2062
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 2062
      if (ssi.eq.ss(j)) m1=j-1
      if (m2.lt.2) m2=2
      m3=mx(ssj+1)
      m4=mn(ssj+1)
      if (m4.gt.n) go to 2062
      ii=ia+wtw(j,is+1)+wty(m2)
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      lad1=ijadd(i*(i-1)/2+j)
      do 2063 k=m2,m1
      if (ssj.eq.ssi) m3=k-1
      if (m3.lt.m4) goto 2063
      jj=jja+wty(m4)
      lsm=xor(ksm,ss(k))
      los=os(lsm+1)
      lad2=lad1+kadd(k+kos)
      cii=c(ii)
      sii=s(ii)
      do 2064 l=m4,m3
      lad3=lad2+ladd(l+los)
c     z=int(lad3+1)+int(lad3+2)
c     sii=sii+z*c(jj)
c     s(jj)=s(jj)+z*cii
      z=(cii*s(jj)+c(jj)*sii)
      int(lad3+1)=int(lad3+1)+z
      int(lad3+2)=int(lad3+2)+z
c
 2064 jj=jj+1
c     s(ii)=s(ii)+sii
 2063 ii=ii+1
 2062 continue
 2061 continue
c
c   no   101a
c
      do 2065 i=m6,m5
      ssj=xor(ss(i),js)
      if (ssj.gt.ss(i)) go to 2065
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) go to 2065
      if (ssj.eq.ss(i)) m1=i-1
      if (m2.lt.3) m2=3
      if (m1.lt.m2) goto 2065
      jj=ja+wtw(i,js+1)+wty(m2)
      do 2066 j=m2,m1
      if (j.lt.jming.or.j.gt.jmaxg) go to 2066
      ksm=xor(ss(i),ss(j))
      kos=os(ksm+1)
      lad1=ijadd(i*(i-1)/2+j)
      cjj=c(jj)
      sjj=s(jj)
      do 2067 k=2,j-1
      ssi=xor(ss(k),is)
      if (ssi.gt.ss(k)) go to 2067
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 2067
      if (ss(k).eq.ssi) m3=k-1
      if (m3.lt.m4) goto 2067
      ii=ia+wtw(k,is+1)+wty(m4)
      lsm=xor(ksm,ss(k))
      los=os(lsm+1)
      lad2=lad1+kadd(k+kos)
      do 2068 l=m4,m3
      lad3=lad2+ladd(l+los)
c     z=int(lad3+1)+int(lad3+3)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=(cjj*s(ii)+c(ii)*sjj)
      int(lad3+1)=int(lad3+1)+z
      int(lad3+3)=int(lad3+3)+z
c
 2068 ii=ii+1
 2067 continue
c     s(jj)=s(jj)+sjj
 2066 jj=jj+1
c2066 continue
 2065 continue
      return
c
c     closed internal loop
c***********************************************************************
  318 continue
      stop ' illegal entry 318 in init4x '
c         original statement wasnt mult by two--see below
cbl   z=(dot(c(ia),s(ja),n1)+dot(c(ja),s(ia),n1))*val1
cbl   int(arr+tr1)=int(arr+tr1)+z
cbl   int(arr+tr2)=int(arr+tr2)+z*val2
cbl   return
c
c     ----- internal case for tracks of (3,2,1) -----
c***********************************************************************
  319 continue
c      stop ' illegal entry 319 in init4x '
c     z=dot(c(ia),s(ja),n1)*val1
      z=(dot(c(ia),s(ja),n1)+dot(c(ja),s(ia),n1))*val1
      int(arr+1)=int(arr+1)+z
      int(arr+2)=int(arr+2)+z
      int(arr+3)=int(arr+3)+z
      return
c
c***********************************************************************
c
  317 z=val1*(int(arr+tr1)+val2*int(arr+tr2))
c   zz   arrival,for completeness
      stop ' zz arrival in external'
cbl   return
      end
