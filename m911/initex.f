*deck @(#)initex.f	5.1  11/6/94
      subroutine initex(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,
     $                  ss)
c
c
c
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
      common/count3/miseg(40,20),mxm,mxseg,msout
c
      real*8 z,val1,val2,val3,zr,zl
      sqrt2=sqrt(2.0d+00)
      sqrt3=sqrt(3.0d+00)
      sqt1p5=sqrt(1.5d+00)
      return
c
c
      entry shapes(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,
     $             ss)
c
c
      aos=os(asm+1)
      mxm=max(mxm,m)
      mxseg=max(mxseg,iseg)
      if((iseg.le.40).and.(m.le.20))then
      miseg(iseg,m)=miseg(iseg,m)+1
      else
      msout=msout+1
      endif
      go to (999,999,304,999,999,310,307,308,303,306,309,302,301,314
c     go to (301,302,303,304,305,306,307,308,309,310,311,312,313,314
c     go to ( yz, xz, xy, yy, xx, zy, yx, yw, wz, wy, wx, xw, ww,---yy
c
     1,314,314,317,318,319),m
c       xx, ww, zz,---4 internals),m
c
      stop 'funny m'
c
c     ----- yz no 1 iseg=22 -----
c***********************************************************************
  301 if (iseg.eq.16) go to 3011
      if (iseg.eq.18) go to 3012
      if (iseg.ne.22) go to 999
      if (js.ne.0) return
      jj=ja
      m1=mx(is+1)
      m2=mn(is+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      ii=ia+wty(m2)
      cjj=c(jj)
      sjj=s(jj)
      do 152 i=m2,m1
      lad=arr+ladd(i+aos)
c     z=val1*(int(lad+tr1)+val2*int(lad+tr2))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad+tr1)=int(lad+tr1)+z
      int(lad+tr2)=int(lad+tr2)+z*val2
c
  152 ii=ii+1
c     s(jj)=s(jj)+sjj
      return
c
c    yz   n0 1  iseg=16
c
 3011 jj=ja
      if (js.ne.0) return
      m1=mx(is+1)
      m2=mn(is+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      ii=ia+wty(m2)
      cjj=c(jj)
      sjj=s(jj)
      do 151 i=m2,m1
      lad=arr+ladd(i+aos)
c     z=val1*(int(lad+3)+int(lad+1))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad+3)=int(lad+3)+z
      int(lad+1)=int(lad+1)+z
c
c      d(jj,ii;idry,i)
      iix=indxs(idry)+i
      int1(iix)=int1(iix)+val1*cjj*s(ii)
      iix=indxs(i)+idry
      int1(iix)=int1(iix)+val1*c(ii)*sjj
  151 ii=ii+1
c     s(jj)=s(jj)+sjj
      return
c
c  yz  no 1  iseg=18
c
 3012 if (js.ne.0) return
      m1=mx(is+1)
      m2=mn(is+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      ii=ia+wty(m2)
      jj=ja
      cjj=c(jj)
      sjj=s(jj)
      do 150 i=m2,m1
      lad=arr+ladd(i+aos)
c     z=val1*int(lad+3)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad+3)=int(lad+3)+val1*(cjj*s(ii)+c(ii)*sjj)
c
c      d(jj,ii;idry,i)
      iix=indxs(idry)+i
      int1(iix)=int1(iix)+val1*cjj*s(ii)
      iix=indxs(i)+idry
      int1(iix)=int1(iix)+val1*c(ii)*sjj
  150 ii=ii+1
c     s(jj)=s(jj)+sjj
      return
c
c   xz  no 2 iseg=13
c***********************************************************************
  302 if (iseg.ne.13) go to 999
      if (js.ne.0) return
      jj=ja
      cjj=c(jj)
      sjj=s(jj)
      if (n.lt.2) return
      do 250 i=2,n
      ssi=xor(ss(i),is)
      if (ssi.gt.ss(i)) go to 250
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 250
      if (ssi.eq.ss(i)) m1=i-1
      if (m1.lt.m2) go to 250
      ii=ia+wtx(i,is+1)+wty(m2)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)
      do 251 j=m2,m1
      lad1=lad+ladd(j+los)
c     z=val1*(int(lad1+1)-int(lad1+3))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad1+1)=int(lad1+1)+z
      int(lad1+3)=int(lad1+3)-z
c
  251 ii=ii+1
  250 continue
c     s(jj)=s(jj)+sjj
      return
c
c   xy   no 3 iseg=18
c***********************************************************************
  303 continue
      if (n.lt.2) return
      if (iseg.eq.16) go to 3031
      if (iseg.eq.3) go to 3033
      if (iseg.eq.22) go to 3034
      if (iseg.ne.18) go to 999
      ssi=xor(js,is)
      if (ssi.gt.js) go to 6601
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6601
      if (m1.gt.n) m1=n
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 6601
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6601
      jj=ja+wty(m2)
      do 365 i=m2,m1
      if (ssi.eq.js) m3=i-1
      if (m3.lt.m4) go to 365
      ii=ia+wtx(i,is+1)+wty(m4)
      cjj=c(jj)
      sjj=s(jj)
      do 366 j=m4,m3
      lad=arr+ladd(j+aos)
c     z=-val1*int(lad+3)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad+3)=int(lad+3)-val1*(cjj*s(ii)+c(ii)*sjj)
c      d(jj,ii;idry,j)
      iix=indxs(idry)+j
      int1(iix)=int1(iix)-val1*cjj*s(ii)
      iix=indxs(j)+idry
      int1(iix)=int1(iix)-val1*c(ii)*sjj
c
  366 ii=ii+1
c     s(jj)=s(jj)+sjj
  365 jj=jj+1
c
c  no 150
c
 6601 ssi=xor(js,is)
      if (js.gt.ssi) return
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) return
      do 5352 i=m2,m1
      if (js.eq.ssi) m3=i-1
      if (m3.lt.m4) go to 5352
      ii=ia+wtx(i,is+1)+wty(m4)
      jj=ja+wty(m4)
c     z=val1*int(arr+ladd(i+aos)+3)
      z=0.0d+00
      zl=0.0d+00
      zr=0.0d+00
      do 5353 j=m4,m3
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=z+(c(jj)*s(ii)+c(ii)*s(jj))
      zr=zr+c(ii)*s(jj)
      zl=zl+c(jj)*s(ii)
c
      ii=ii+1
 5353 jj=jj+1
      int(arr+ladd(i+aos)+3)=int(arr+ladd(i+aos)+3)+z*val1
c      d(jj,ii;idry,i)
      iix=indxs(idry)+i
      int1(iix)=int1(iix)+val1*zl
      iix=indxs(i)+idry
      int1(iix)=int1(iix)+val1*zr
c
 5352 continue
      return
c
c  xy  no 3 iseg=16
c
 3031 ssi=xor(js,is)
      if (ssi.gt.js) go to 6602
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6602
      if (m1.gt.n) m1=n
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 6602
      if (m2.lt.2) m2=2
      jj=ja+wty(m2)
      do 363 i=m2,m1
      cjj=c(jj)
      sjj=s(jj)
      if (ssi.eq.js) m3=i-1
      if (m3.lt.m4) go to 363
      ii=ia+wtx(i,is+1)+wty(m4)
      do 364 j=m4,m3
      lad=arr+ladd(j+aos)
c     z=-val1*(int(lad+3)+int(lad+1))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=-val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad+3)=int(lad+3)+z
      int(lad+1)=int(lad+1)+z
c      d(jj,ii;idry,j)
      iix=indxs(idry)+j
      int1(iix)=int1(iix)-val1*cjj*s(ii)
      iix=indxs(j)+idry
      int1(iix)=int1(iix)-val1*c(ii)*sjj
c
  364 ii=ii+1
c     s(jj)=s(jj)+sjj
  363 jj=jj+1
c
c  no 150
c
 6602 ssi=xor(js,is)
      if (js.gt.ssi) return
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) return
      do 6352 i=m2,m1
      if (js.eq.ssi) m3=i-1
      if (m3.lt.m4) go to 6352
      ii=ia+wtx(i,is+1)+wty(m4)
      jj=ja+wty(m4)
      lad=arr+ladd(i+aos)
c     z=val1*(int(lad+3)+int(lad+1))
      z=0.0d+00
      zr=0.d0
      zl=0.d0
      do 6353 j=m4,m3
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=z+c(jj)*s(ii)+c(ii)*s(jj)
      zl=zl+c(jj)*s(ii)
      zr=zr+c(ii)*s(jj)
      ii=ii+1
 6353 jj=jj+1
      z=z*val1
      zl=zl*val1
      zr=zr*val1
      int(lad+3)=int(lad+3)+z
      int(lad+1)=int(lad+1)+z
c      d(jj,ii;idry,i)
      iix=indxs(idry)+i
      int1(iix)=int1(iix)+zl
      iix=indxs(i)+idry
      int1(iix)=int1(iix)+zr
c
 6352 continue
      return
c
c  xy  no 3 iseg=22
c
 3034 ssi=xor(js,is)
      if (ssi.gt.js) go to 6603
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6603
      if (m1.gt.n) m1=n
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 6603
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6603
      jj=ja+wty(m2)
      do 380 i=m2,m1
      if (ssi.eq.js) m3=i-1
      if (m3.lt.m4) go to 380
      ii=ia+wtx(i,is+1)+wty(m4)
      cjj=c(jj)
      sjj=s(jj)
      do 381 j=m4,m3
      lad=arr+ladd(j+aos)
c     z=-val1*(int(lad+tr1)+val2*int(lad+tr2))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=-val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad+tr1)=int(lad+tr1)+z
      int(lad+tr2)=int(lad+tr2)+z*val2
c
  381 ii=ii+1
c     s(jj)=s(jj)+sjj
  380 jj=jj+1
c
c no 150
c
 6603 ssi=xor(js,is)
      if (js.gt.ssi) return
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) return
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) return
      do 7352 i=m2,m1
      if (js.eq.ssi) m3=i-1
      if (m3.lt.m4) go to 7352
      ii=ia+wtx(i,is+1)+wty(m4)
      jj=ja+wty(m4)
      lad=arr+ladd(i+aos)
c     z=val1*(int(lad+tr1)+val2*int(lad+tr2))
      z=0.0d+00
      do 7353 j=m4,m3
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=z+c(jj)*s(ii)+c(ii)*s(jj)
c
      ii=ii+1
 7353 jj=jj+1
      z=z*val1
      int(lad+tr1)=int(lad+tr1)+z
      int(lad+tr2)=int(lad+tr2)+z*val2
c
 7352 continue
      return
c
c  xy  no 3  iseg=3
c
 3033 ssi=xor(js,is)
      if (ssi.gt.js) go to 6604
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6604
      if (m1.gt.n) m1=n
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 6604
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6604
      jj=ja+wty(m2)
      do 350 i=m2,m1
      if (i.lt.jming.or.i.gt.jmaxg) go to 350
      if (ssi.eq.js) m3=i-1
      if (m3.lt.m4) go to 350
      ii=ia+wtx(i,is+1)+wty(m4)
      ksm=xor(asm,ss(i))
      kos=os(ksm+1)
      lad=ijadd(arr+i)+kadd(i+kos)
      cjj=c(jj)
      sjj=s(jj)
      do 351 j=m4,m3
      lad1=lad+ladd(j+aos)
c     z=-val1*(int(lad1+2)-int(lad1+1))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=-val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad1+2)=int(lad1+2)+z
      int(lad1+1)=int(lad1+1)-z
c
  351 ii=ii+1
c     s(jj)=s(jj)+sjj
  350 jj=jj+1
c
c  no 4
c
 6604 ssi=xor(js,is)
      if (js.gt.ssi) go to 6605
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 6605
      if (m1.gt.n) m1=n
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) go to 6605
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6605
      do 352 i=m2,m1
      if (i.lt.jming.or.i.gt.jmaxg) go to 352
      if (js.eq.ssi) m3=i-1
      if (m3.lt.m4) go to 352
      jj=ja+wty(m4)
      ii=ia+wtx(i,is+1)+wty(m4)
      ksm=xor(asm,ss(i))
      kos=os(ksm+1)
      lad=ijadd(arr+i)
      do 353 j=m4,m3
      lsm=xor(ksm,ss(j))
      los=os(lsm+1)
      lad1=lad+kadd(j+kos)+ladd(j+los)
c     z=-val1*(int(lad1+1)-int(lad1+2))
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=-val1*(c(jj)*s(ii)+c(ii)*s(jj))
      int(lad1+1)=int(lad1+1)+z
      int(lad1+2)=int(lad1+2)-z
c
      ii=ii+1
  353 jj=jj+1
  352 continue
c
c    no   5
c
 6605 do 354 i=3,n
      if (i.lt.jming.or.i.gt.jmaxg) go to 354
      ssi=xor(ss(i),is)
      if (js.gt.ssi) go to 354
      if (ssi.gt.ss(i)) go to 354
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 354
      if (ssi.eq.ss(i).or.m1.gt.n) m1=i-1
      if (m2.lt.2) m2=2
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) go to 354
      ii=ia+wtx(i,is+1)+wty(m2)
      ksm=xor(asm,ss(i))
      kos=os(ksm+1)
      lad=ijadd(arr+i)
      do 355 j=m2,m1
      if (js.eq.ssi) m3=j-1
      if (m3.lt.m4) go to 355
      jj=ja+wty(m4)
      lsm=xor(ksm,ss(j))
      los=os(lsm+1)
      lad1=lad+kadd(j+kos)
      cii=c(ii)
      sii=s(ii)
      do 356 k=m4,m3
      lad2=lad1+ladd(k+los)
c     z=-val1*(int(lad2+1)-int(lad2+2))
c     s(jj)=s(jj)+z*cii
c     sii=sii+z*c(jj)
      z=-val1*(cii*s(jj)+c(jj)*sii)
      int(lad2+1)=int(lad2+1)+z
      int(lad2+2)=int(lad2+2)-z
c
  356 jj=jj+1
c     s(ii)=s(ii)+sii
  355 ii=ii+1
  354 continue
c
c    no   6
c
      do 357 i=3,n
      if (i.lt.jming.or.i.gt.jmaxg) go to 357
      ssi=xor(ss(i),is)
      if (js.gt.ss(i)) go to 357
      if (ssi.gt.js) go to 357
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 357
      if (js.eq.ss(i).or.m1.gt.n) m1=i-1
      if (m2.lt.2) m2=2
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 357
      jj=ja+wty(m2)
      iia=ia+wtx(i,is+1)
      ksm=xor(asm,ss(i))
      kos=os(ksm+1)
      lad=ijadd(arr+i)
      do 358 j=m2,m1
      if (ssi.eq.js) m3=j-1
      if (m3.lt.m4) go to 358
      ii=iia+wty(m4)
      lsm=xor(ksm,ss(j))
      los=os(lsm+1)
      lad1=lad+kadd(j+kos)
      cjj=c(jj)
      sjj=s(jj)
      do 359 k=m4,m3
      lad2=lad1+ladd(k+los)
c     z=-val1*(int(lad2+3)-int(lad2+2))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=-val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad2+3)=int(lad2+3)+z
      int(lad2+2)=int(lad2+2)-z
c
  359 ii=ii+1
c     s(jj)=s(jj)+sjj
  358 jj=jj+1
  357 continue
c
c   no 7
c
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.3) m2=3
      if (m1.lt.m2) return
      jj=ja+wty(m2)
      do 360 i=m2,m1
      if (i.lt.jming.or.i.gt.jmaxg) go to 360
      ksm=xor(asm,ss(i))
      kos=os(ksm+1)
      lad=ijadd(arr+i)
      cjj=c(jj)
      sjj=s(jj)
      do 361 j=2,i-1
      ssi=xor(is,ss(j))
      if (ssi.gt.ss(j))go to 361
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 361
      if (ssi.eq.ss(j)) m3=j-1
      if (m3.lt.m4) go to 361
      ii=ia+wtx(j,is+1)+wty(m4)
      lsm=xor(ksm,ss(j))
      los=os(lsm+1)
      lad1=lad+kadd(j+kos)
      do 362 k=m4,m3
      lad2=lad1+ladd(k+los)
c     z=val1*(int(lad2+1)-int(lad2+3))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad2+1)=int(lad2+1)+z
      int(lad2+3)=int(lad2+3)-z
c
  362 ii=ii+1
  361 continue
c     s(jj)=s(jj)+sjj
  360 jj=jj+1
c 360 continue
      return
c
c   yy   no 9  iseg=4
c***********************************************************************
  304 if (iseg.eq.5) go to 3043
      if (iseg.eq.6) go to 3041
      if (iseg.eq.8) go to 3042
      if (iseg.ne.4) go to 999
      val4=-val1/sqrt2
      if (is.gt.js) go to 6606
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6606
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      m3=mx(is+1)
      m4=mn(is+1)
      if (m4.gt.n) go to 6606
      jj=ja+wty(m2)
      do 450 i=m2,m1
      if (is.eq.js) m3=i-1
      if (m3.lt.m4) go to 450
      ii=ia+wty(m4)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)
      cjj=c(jj)
      sjj=s(jj)
      do 451 j=m4,m3
      lad1=lad+ladd(j+los)+1
c     z=val4*(int(lad1)-2.*int(lad1+1))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=val4*(cjj*s(ii)+c(ii)*sjj)
      int(lad1)=int(lad1)+z
      int(lad1+1)=int(lad1+1)-2*z
c
  451 ii=ii+1
c     s(jj)=s(jj)+sjj
  450 jj=jj+1
c
c   no 10
c
 6606 if(ia.eq.ja)return
      if (is.ne.js) return
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      ii=ia+wty(m2)
      jj=ja+wty(m2)
      do 460 i=m2,m1
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+ladd(i+los)+kadd(i+aos)
c     z=val4*(int(lad+1)-2.*int(lad+2))
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=val4*(c(jj)*s(ii)+c(ii)*s(jj))
      int(lad+1)=int(lad+1)+z
      int(lad+2)=int(lad+2)-2*z
c
      ii=ii+1
  460 jj=jj+1
      return
c
c yy  no 9 iseg=8
c
 3042 val4=val1*sqt1p5
      if (is.gt.js) go to 6607
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6607
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6607
      m3=mx(is+1)
      m4=mn(is+1)
      if (m4.gt.n) go to 6607
      jj=ja+wty(m2)
      do 420 i=m2,m1
      if (is.eq.js) m3=i-1
      if (m3.lt.m4) go to 420
      ii=ia+wty(m4)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)
      cjj=c(jj)
      sjj=s(jj)
      do 421 j=m4,m3
      lad1=lad+ladd(j+los)+tr1
c     z=val4*int(lad1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad1)=int(lad1)+val4*(cjj*s(ii)+c(ii)*sjj)
c
  421 ii=ii+1
c     s(jj)=s(jj)+sjj
  420 jj=jj+1
c
c  no 10
c
 6607 if(ia.eq.ja) return
      if (is.ne.js) go to 6598
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6598
      if (m1.gt.n) m1=n
      ii=ia+wty(m2)
      jj=ja+wty(m2)
      do 430 i=m2,m1
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+ladd(i+los)+kadd(i+aos)
c     z=val4*int(lad+1)
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      int(lad+1)=int(lad+1)+val4*(c(jj)*s(ii)+c(ii)*s(jj))
c
      ii=ii+1
  430 jj=jj+1
c
c  no 8
c
 6598 if (js.gt.is) return
      m1=mx(is+1)
      m2=mn(is+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) return
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) return
      ii=ia+wty(m2)
      do 440 i=m2,m1
      if (is.eq.js) m3=i-1
      if (m3.lt.m4) go to 440
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)
      jj=ja+wty(m4)
      do 441 j=m4,m3
      lad1=lad+ladd(j+los)
c     z=val4*int(lad1+1)
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      int(lad1+1)=int(lad1+1)+val4*(c(jj)*s(ii)+c(ii)*s(jj))
c
  441 jj=jj+1
  440 ii=ii+1
      return
c
c  yy  no 10  iseg=5
c
 3043 continue
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      ii=ia+wty(m2)
      jj=ja+wty(m2)
      xx=-val1/sqrt2+val3*sqt1p5
      if (abs(xx).lt.1.0d-06) go to 490
      val5=sqrt2*val1
      if(ia.eq.ja) go to 491
      if (is.ne.js) go to 491
      do 493 i=m2,m1
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+ladd(i+los)+kadd(i+aos)
c     z=xx*int(lad+1)+val5*int(lad+2)
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=c(ii)*s(jj)+c(jj)*s(ii)
      int(lad+1)=int(lad+1)+xx*z
      int(lad+2)=int(lad+2)+val5*z
c
      ii=ii+1
  493 jj=jj+1
      go to 491
  490 val5=sqrt2*val1
      if(ia.eq.ja) go to 492
      if (is.ne.js) go to 492
      do 494 i=m2,m1
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+ladd(i+los)+kadd(i+aos)
c     z=val5*int(lad+2)
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      int(lad+2)=int(lad+2)+val5*(c(jj)*s(ii)+c(ii)*s(jj))
c
      ii=ii+1
  494 jj=jj+1
      go to 492
c
c   no   9
c
  491 if (is.gt.js) return
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) return
      m3=mx(is+1)
      m4=mn(is+1)
      if (m4.gt.n) return
      jj=ja+wty(m2)
      do 495 i=m2,m1
      if (is.eq.js) m3=i-1
      if (m3.lt.m4) go to 495
      ii=ia+wty(m4)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)
      cjj=c(jj)
      sjj=s(jj)
      do 496 j=m4,m3
      lad1=lad+ladd(j+los)
c     z=xx*int(lad1+1)+val5*int(lad1+2)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=cjj*s(ii)+c(ii)*sjj
      int(lad1+1)=int(lad1+1)+z*xx
      int(lad1+2)=int(lad1+2)+z*val5
c
  496 ii=ii+1
c     s(jj)=s(jj)+sjj
  495 jj=jj+1
      return
c
  492 if (is.gt.js) return
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) return
      m3=mx(is+1)
      m4=mn(is+1)
      if (m4.gt.n) return
      jj=ja+wty(m2)
      do 497 i=m2,m1
      if (is.eq.js) m3=i-1
      if (m3.lt.m4) go to 497
      ii=ia+wty(m4)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)
      cjj=c(jj)
      sjj=s(jj)
      do 498 j=m4,m3
      lad1=lad+ladd(j+los)+2
c     z=val5*int(lad1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad1)=int(lad1)+val5*(cjj*s(ii)+c(ii)*sjj)
c
  498 ii=ii+1
c     s(jj)=s(jj)+sjj
  497 jj=jj+1
      return
c
c   yy   no 8  iseg=6
c
 3041 xx=-val1/sqrt2+val3*sqt1p5
      if (abs(xx).lt.1.0d-06) go to 480
      val5=val1*sqrt2
      if (js.gt.is) go to 481
      m1=mx(is+1)
      m2=mn(is+1)
      if (m2.gt.n) go to 481
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 481
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) go to 481
      ii=ia+wty(m2)
      do 472 i=m2,m1
      if (is.eq.js) m3=i-1
      if (m3.lt.m4) go to 472
      jj=ja+wty(m4)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)
      cii=c(ii)
      sii=s(ii)
      do 473 j=m4,m3
      lad1=lad+ladd(j+los)
c     z=xx*int(lad1+1)+val5*int(lad1+2)
c     sii=sii+z*c(jj)
c     s(jj)=s(jj)+z*cii
      z=cii*s(jj)+c(jj)*sii
      int(lad1+1)=int(lad1+1)+z*xx
      int(lad1+2)=int(lad1+2)+z*val5
c
  473 jj=jj+1
c     s(ii)=s(ii)+sii
  472 ii=ii+1
      go to 481
c
  480 val5=sqrt2*val1
      if (js.gt.is) go to 483
      m1=mx(is+1)
      m2=mn(is+1)
      if (m2.gt.n) go to 483
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 483
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) go to 483
      ii=ia+wty(m2)
      do 474 i=m2,m1
      if (is.eq.js) m3=i-1
      if (m3.lt.m4) go to 474
      jj=ja+wty(m4)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)
      cii=c(ii)
      sii=s(ii)
      do 475 j=m4,m3
      lad1=lad+ladd(j+los)
c     z=val5*int(lad1+2)
c     sii=sii+z*c(jj)
c     s(jj)=s(jj)+z*cii
      int(lad1+2)=int(lad1+2)+val5*(cii*s(jj)+c(jj)*sii)
c
  475 jj=jj+1
c     s(ii)=s(ii)+sii
  474 ii=ii+1
      go to 483
c
c    no  10
c
  481 if (is.ne.js) go to 482
      if(ia.eq.ja) go to 482
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 482
      if (m1.gt.n) m1=n
      ii=ia+wty(m2)
      jj=ja+wty(m2)
      do 470 i=m2,m1
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+ladd(i+los)+kadd(i+aos)
c     z=xx*int(lad+1)+val5*int(lad+2)
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=s(ii)*c(jj)+s(jj)*c(ii)
      int(lad+1)=int(lad+1)+z*xx
      int(lad+2)=int(lad+2)+z*val5
c
      ii=ii+1
  470 jj=jj+1
      go to 482
c
  483 if (is.ne.js) go to 484
      if (ia.eq.ja) go to 484
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 484
      if (m1.gt.n) m1=n
      ii=ia+wty(m2)
      jj=ja+wty(m2)
      do 471 i=m2,m1
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+ladd(i+los)+kadd(i+aos)
c     z=val5*int(lad+2)
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      int(lad+2)=int(lad+2)+val5*(c(ii)*s(jj)+c(jj)*s(ii))
c
      ii=ii+1
  471 jj=jj+1
      go to 484
c
c  no  9
c
  482 if (is.gt.js) return
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) return
      m3=mx(is+1)
      m4=mn(is+1)
      if (m4.gt.n) return
      jj=ja+wty(m2)
      do 476 i=m2,m1
      if (is.eq.js) m3=i-1
      if (m3.lt.m4) go to 476
      ii=ia+wty(m4)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)
      cjj=c(jj)
      sjj=s(jj)
      do 477 j=m4,m3
      lad1=lad+ladd(j+los)
c     z=xx*int(lad1+3)+val5*int(lad1+2)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=cjj*s(ii)+c(ii)*sjj
      int(lad1+3)=int(lad1+3)+z*xx
      int(lad1+2)=int(lad1+2)+z*val5
c
  477 ii=ii+1
c     s(jj)=s(jj)+sjj
  476 jj=jj+1
      return
c
  484 if (is.gt.js) return
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) return
      m3=mx(is+1)
      m4=mn(is+1)
      if (m4.gt.n) return
      jj=ja+wty(m2)
      do 478 i=m2,m1
      if (is.eq.js) m3=i-1
      if (m3.lt.m4) go to 478
      ii=ia+wty(m4)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)
      cjj=c(jj)
      sjj=s(jj)
      do 479 j=m4,m3
      lad1=lad+ladd(j+los)
c     z=val5*int(lad1+2)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
c         following line appears to be wrong
c     int(lad1+2)=int(lad1+2)+2*val5*cjj*c(ii)
      int(lad1+2)=int(lad1+2)+val5*(cjj*s(ii)+c(ii)*sjj)
c
  479 ii=ii+1
c     s(jj)=s(jj)+sjj
  478 jj=jj+1
      return
c
c   wy  no 27  iseg=2
c***********************************************************************
  310 if (iseg.eq.15) go to 3101
      if (iseg.eq.17) go to 3102
      if (iseg.eq.21) go to 3103
      if (iseg.ne.2) go to 999
c
c   no 32
c
c
c
      val4=val1*sqrt2
      if (is.ne.0) go to 6587
      do 1001 i=2,n
      if (i.lt.jming.or.i.gt.jmaxg) go to 1001
      if (js.gt.ss(i)) go to 1001
      ii=ia+wab(i)
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 1001
      if (js.eq.ss(i).or.m1.gt.n) m1=i-1
      jj=ja+wty(m2)
      ksm=xor(asm,ss(i))
      kos=os(ksm+1)
      lad=ijadd(arr+i)+kadd(i+kos)
      cii=c(ii)
      sii=s(ii)
      do 1002 j=m2,m1
      lad1=lad+ladd(j+aos)
c     z=int(lad1+1)*val4
c     sii=sii+z*c(jj)
c     s(jj)=s(jj)+z*cii
      int(lad1+1)=int(lad1+1)+val4*(cii*s(jj)+c(jj)*sii)
c
 1002 jj=jj+1
c     s(ii)=s(ii)+sii
 1001 continue
c
c   no  28
c
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6599
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6599
      jj=ja+wty(m2)
      do 1003 i=m2,m1
      if (i.lt.jming.or.i.gt.jmaxg) go to 1003
      ksm=xor(asm,ss(i))
      kos=os(ksm+1)
      lad=ijadd(arr+i)
      cjj=c(jj)
      sjj=s(jj)
      do 1004 j=1,i-1
      ii=ia+wab(j)
      lsm=xor(ksm,ss(j))
      los=os(lsm+1)
      lad1=lad+kadd(j+kos)+ladd(j+los)
c     z=val4*int(lad1+1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad1+1)=int(lad1+1)+val4*(cjj*s(ii)+c(ii)*sjj)
c
 1004 continue
c     s(jj)=s(jj)+sjj
 1003 jj=jj+1
c
c   no 30
c
 6587 m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6599
      if (m1.gt.n) m1=n
      if (m2.lt.3) m2=3
      if (m1.lt.m2) go to 6599
      jj=ja+wty(m2)
      do 1005 i=m2,m1
      if (i.lt.jming.or.i.gt.jmaxg) go to 1005
      ksm=xor(asm,ss(i))
      kos=os(ksm+1)
      lad=ijadd(arr+i)
      cjj=c(jj)
      sjj=s(jj)
      do 1006 j=2,i-1
      ssi=xor(is,ss(j))
      if (ssi.gt.ss(j))go to 1006
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 1006
      if (ssi.eq.ss(j)) m3=j-1
      if (m3.lt.m4) go to 1006
      ii=ia+wtw(j,is+1)+wty(m4)
      lsm=xor(ksm,ss(j))
      los=os(lsm+1)
      lad1=lad+kadd(j+kos)
      do 1007 k=m4,m3
      lad2=lad1+ladd(k+los)
c     z=val1*(int(lad2+1)+int(lad2+3))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad2+1)=int(lad2+1)+z
      int(lad2+3)=int(lad2+3)+z
c
 1007 ii=ii+1
 1006 continue
c     s(jj)=s(jj)+sjj
 1005 jj=jj+1
c1005 continue
c
c   no 32
c
 6599 if(n.lt.3) go to 1018
      do 1008 i=3,n
      if (i.lt.jming.or.i.gt.jmaxg) go to 1008
      ssi=xor(ss(i),is)
      if (js.gt.ssi) go to 1008
      if (ssi.gt.ss(i)) go to 1008
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 1008
      if (ssi.eq.ss(i).or.m1.gt.n) m1=i-1
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 1008
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) go to 1008
      ii=ia+wtw(i,is+1)+wty(m2)
      ksm=xor(asm,ss(i))
      kos=os(ksm+1)
      lad=ijadd(arr+i)
      do 1009 j=m2,m1
      if (js.eq.ssi) m3=j-1
      if (m3.lt.m4) go to 1009
      jj=ja+wty(m4)
      lsm=xor(ksm,ss(j))
      los=os(lsm+1)
      lad1=lad+kadd(j+kos)
      cii=c(ii)
      sii=s(ii)
      do 1010 k=m4,m3
      lad2=lad1+ladd(k+los)
c     z=val1*(int(lad2+1)+int(lad2+2))
c     sii=sii+z*c(jj)
c     s(jj)=s(jj)+z*cii
      z=val1*(cii*s(jj)+c(jj)*sii)
      int(lad2+1)=int(lad2+1)+z
      int(lad2+2)=int(lad2+2)+z
c
 1010 jj=jj+1
c     s(ii)=s(ii)+sii
 1009 ii=ii+1
 1008 continue
 1018 continue
c
c   no 33
c
      ssi=xor(js,is)
      if (ssi.gt.js)go to  6608
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6608
      if (m1.gt.n) m1=n
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 6608
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6608
      jj=ja+wty(m2)
      do 1011 i=m2,m1
      if (i.lt.jming.or.i.gt.jmaxg) go to 1011
      if (ssi.eq.js) m3=i-1
      if (m3.lt.m4) go to 1011
      ii=ia+wtw(i,is+1)+wty(m4)
      ksm=xor(asm,ss(i))
      kos=os(ksm+1)
      lad=ijadd(arr+i)+kadd(i+kos)
      cjj=c(jj)
      sjj=s(jj)
      do 1012 j=m4,m3
      lad1=lad+ladd(j+aos)
c     z=val1*(int(lad1+2)+int(lad1+1))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad1+2)=int(lad1+2)+z
      int(lad1+1)=int(lad1+1)+z
c
 1012 ii=ii+1
c     s(jj)=s(jj)+sjj
 1011 jj=jj+1
c
c   no 36
c
 6608 ssi=xor(js,is)
      if (js.gt.ssi) go to 6609
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 6609
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) go to 6609
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6609
      do 1013 i=m2,m1
      if (i.lt.jming.or.i.gt.jmaxg) go to 1013
      if (js.eq.ssi) m3=i-1
      if (m3.lt.m4) go to 1013
      jj=ja+wty(m4)
      ii=ia+wtw(i,is+1)+wty(m4)
      ksm=xor(asm,ss(i))
      kos=os(ksm+1)
      lad=ijadd(arr+i)
      do 1014 j=m4,m3
      lsm=xor(ksm,ss(j))
      los=os(lsm+1)
      lad1=lad+kadd(j+kos)+ladd(j+los)
c     z=val1*(int(lad1+1)+int(lad1+2))
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=val1*(c(ii)*s(jj)+c(jj)*s(ii))
      int(lad1+1)=int(lad1+1)+z
      int(lad1+2)=int(lad1+2)+z
c
      ii=ii+1
 1014 jj=jj+1
 1013 continue
c
c    no 38
c
 6609 do 1015 i=3,n
      if (i.lt.jming.or.i.gt.jmaxg) go to 1015
      ssi=xor(ss(i),is)
      if (js.gt.ss(i)) go to 1015
      if (ssi.gt.js) go to 1015
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 1015
      if (js.eq.ss(i).or.m1.gt.n) m1=i-1
      if (m2.lt.2) m2=2
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 1015
      jj=ja+wty(m2)
      iia=ia+wtw(i,is+1)
      ksm=xor(asm,ss(i))
      kos=os(ksm+1)
      lad=ijadd(arr+i)
      do 1016 j=m2,m1
      if (ssi.eq.js) m3=j-1
      if (m3.lt.m4) go to 1016
      ii=iia+wty(m4)
      lsm=xor(ksm,ss(j))
      los=os(lsm+1)
      lad1=lad+kadd(j+kos)
      cjj=c(jj)
      sjj=s(jj)
      do 1017 k=m4,m3
      lad2=lad1+ladd(k+los)
c     z=val1*(int(lad2+3)+int(lad2+2))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad2+3)=int(lad2+3)+z
      int(lad2+2)=int(lad2+2)+z
c
 1017 ii=ii+1
c     s(jj)=s(jj)+sjj
 1016 jj=jj+1
 1015 continue
      return
c
c   wy   no 26   iseg=21
c
 3103 val4=val1*sqrt2
      if (is.ne.0) go to 6610
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6610
      if (m1.gt.n) m1=n
      jj=ja+wty(m2)
      do 1020 i=m2,m1
      ii=ia+wab(i)
      lad=arr+ladd(i+aos)
c     z=val4*(int(lad+tr1)+val2*int(lad+tr2))
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=val4*(c(ii)*s(jj)+c(jj)*s(ii))
      int(lad+tr1)=int(lad+tr1)+z
      int(lad+tr2)=int(lad+tr2)+z*val2
c
 1020 jj=jj+1
c
c   no 33
c
 6610 ssi=xor(js,is)
      if (ssi.gt.js) go to 6611
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6611
      if (m1.gt.n) m1=n
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 6611
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6611
      jj=ja+wty(m2)
      do 1021 i=m2,m1
      if (ssi.eq.js) m3=i-1
      if (m3.lt.m4) go to 1021
      ii=ia+wtw(i,is+1)+wty(m4)
      cjj=c(jj)
      sjj=s(jj)
      do 1022 j=m4,m3
      lad=arr+ladd(j+aos)
c     z=val1*(int(lad+tr1)+val2*int(lad+tr2))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad+tr1)=int(lad+tr1)+z
      int(lad+tr2)=int(lad+tr2)+z*val2
c
 1022 ii=ii+1
c     s(jj)=s(jj)+sjj
 1021 jj=jj+1
c
c no 149
c
 6611 ssi=xor(js,is)
      if (js.gt.ssi) return
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) return
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) return
      do 6013 i=m2,m1
      if (js.eq.ssi) m3=i-1
      if (m3.lt.m4) go to 6013
      ii=ia+wtw(i,is+1)+wty(m4)
      jj=ja+wty(m4)
      lad=arr+ladd(i+aos)
c     z=val1*(int(lad+tr1)+val2*int(lad+tr2))
      z=0.0d+00
      do 6014 j=m4,m3
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=z+c(ii)*s(jj)+c(jj)*s(ii)
      ii=ii+1
 6014 jj=jj+1
      z=z*val1
      int(lad+tr1)=int(lad+tr1)+z
      int(lad+tr2)=int(lad+tr2)+z*val2
c
 6013 continue
      return
c
c    wy  no 26   iseg=15
c
 3101 val4=val1*sqrt2
      if (is.ne.0) go to 6612
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6612
      if (m1.gt.n) m1=n
      jj=ja+wty(m2)
      do 1030 i=m2,m1
      ii=ia+wab(i)
      lad=arr+ladd(i+aos)
c     z=val4*(int(lad+1)+int(lad+2)+int(lad+3))
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=val4*(c(ii)*s(jj)+c(jj)*s(ii))
      int(lad+1)=int(lad+1)+z
      int(lad+2)=int(lad+2)+z
      int(lad+3)=int(lad+3)+z
c
c         d(jj,ii;idry,i)
      iix=indxs(idry)+i
      int1(iix)=int1(iix)+val4*c(jj)*s(ii)
      iix=indxs(i)+idry
      int1(iix)=int1(iix)+val4*c(ii)*s(jj)
 1030 jj=jj+1
c
c    no 33
c
 6612 ssi=xor(js,is)
      if (ssi.gt.js) go to 6613
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6613
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 6613
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6613
      jj=ja+wty(m2)
      do 1031 i=m2,m1
      if (ssi.eq.js) m3=i-1
      if (m3.lt.m4) go to 1031
      ii=ia+wtw(i,is+1)+wty(m4)
      cjj=c(jj)
      sjj=s(jj)
      do 1032 j=m4,m3
      lad=arr+ladd(j+aos)
c     z=val1*(int(lad+3)+int(lad+1))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad+3)=int(lad+3)+z
      int(lad+1)=int(lad+1)+z
c
c         d(jj,ii;idry,j)
      iix=indxs(idry)+j
      int1(iix)=int1(iix)+val1*cjj*s(ii)
      iix=indxs(j)+idry
      int1(iix)=int1(iix)+val1*c(ii)*sjj
 1032 ii=ii+1
c     s(jj)=s(jj)+sjj
 1031 jj=jj+1
c
c no 149
c
 6613 ssi=xor(js,is)
      if (js.gt.ssi) return
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) return
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) return
      do 7013 i=m2,m1
      if (js.eq.ssi) m3=i-1
      if (m3.lt.m4) go to 7013
      ii=ia+wtw(i,is+1)+wty(m4)
      jj=ja+wty(m4)
      lad=arr+ladd(i+aos)
c     z=val1*(int(lad+3)+int(lad+1))
      z=0.0d+00
      zr=z
      zl=z
      do 7014 j=m4,m3
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=z+c(ii)*s(jj)+c(jj)*s(ii)
      zl=zl+c(ii)*s(jj)
      zr=zr+c(jj)*s(ii)
      ii=ii+1
 7014 jj=jj+1
      z=z*val1
      zr=zr*val1
      zl=zl*val1
      int(lad+3)=int(lad+3)+z
      int(lad+1)=int(lad+1)+z
c
c         d(jj,ii;idry,i)
      iix=indxs(idry)+i
      int1(iix)=int1(iix)+zr
      iix=indxs(i)+idry
      int1(iix)=int1(iix)+zl
 7013 continue
      return
c
c   wy  no 26 iseg=17
c
 3102 val4=val1*sqrt2
      if (is.ne.0) go to 6614
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6614
      if (m1.gt.n) m1=n
      jj=ja+wty(m2)
      do 1023 i=m2,m1
      ii=ia+wab(i)
      lad=arr+ladd(i+aos)
c     z=val4*(int(lad+3)+int(lad+2))
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=val4*(c(ii)*s(jj)+c(jj)*s(ii))
      int(lad+3)=int(lad+3)+z
      int(lad+2)=int(lad+2)+z
c      d(jj,ii;idry,i)
      iix=indxs(idry)+i
      int1(iix)=int1(iix)+val4*c(jj)*s(ii)
      iix=indxs(i)+idry
      int1(iix)=int1(iix)+val4*c(ii)*s(jj)
c
 1023 jj=jj+1
c
c   no  33
c
 6614 ssi=xor(js,is)
      if (ssi.gt.js) go to 6615
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6615
      m3=mx(ssi+1)
      m4=mn(ssi+1)
      if (m4.gt.n) go to 6615
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6615
      jj=ja+wty(m2)
      do 1024 i=m2,m1
      if (ssi.eq.js) m3=i-1
      if (m3.lt.m4) goto 1024
      ii=ia+wtw(i,is+1)+wty(m4)
      cjj=c(jj)
      sjj=s(jj)
      do 1025 j=m4,m3
      lad=arr+ladd(j+aos)
c     z=val1*int(lad+3)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad+3)=int(lad+3)+val1*(cjj*s(ii)+c(ii)*sjj)
c      d(jj,ii;idry,j)
      iix=indxs(idry)+j
      int1(iix)=int1(iix)+val1*cjj*s(ii)
      iix=indxs(j)+idry
      int1(iix)=int1(iix)+val1*c(ii)*sjj
c
 1025 ii=ii+1
c     s(jj)=s(jj)+sjj
 1024 jj=jj+1
c
c no 149
c
 6615 ssi=xor(js,is)
      if (js.gt.ssi) return
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) return
      m3=mx(js+1)
      m4=mn(js+1)
      if (m4.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) return
      do 5013 i=m2,m1
      if (js.eq.ssi) m3=i-1
      if (m3.lt.m4) go to 5013
      ii=ia+wtw(i,is+1)+wty(m4)
      jj=ja+wty(m4)
      lad=arr+ladd(i+aos)
c     z=val1*int(lad+3)
      z=0.0d+00
      zl=0.d0
      zr=0.d0
      do 5014 j=m4,m3
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=z+c(ii)*s(jj)+c(jj)*s(ii)
      zl=zl+c(jj)*s(ii)
      zr=zr+c(ii)*s(jj)
c
      ii=ii+1
 5014 jj=jj+1
      int(lad+3)=int(lad+3)+z*val1
c      d(jj,ii;idry,i)
      iix=indxs(idry)+i
      int1(iix)=int1(iix)+val1*zl
      iix=indxs(i)+idry
      int1(iix)=int1(iix)+val1*zr
c
 5013 continue
      return
c
c   zy   no 17   iseg=20
c***********************************************************************
  306 if (iseg.ne.20) go to 999
      if (is.ne.0) return
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      jj=ja+wty(m2)
      ii=ia
      cii=c(ii)
      sii=s(ii)
      do 601 i=m2,m1
      lad=arr+ladd(i+aos)
c     z=val1*(int(lad+tr1)+val2*int(lad+tr2))
c     sii=sii+z*c(jj)
c     s(jj)=s(jj)+z*cii
      z=val1*(cii*s(jj)+c(jj)*sii)
      int(lad+tr1)=int(lad+tr1)+z
      int(lad+tr2)=int(lad+tr2)+z*val2
c
  601 jj=jj+1
c     s(ii)=s(ii)+sii
      return
c
c     yx   no  18   iseg=20
c***********************************************************************
  307 if (iseg.ne.20) go to 999
      ssj=xor(is,js)
      if (ssj.gt.is) go to 6630
      m1=mx(is+1)
      m2=mn(is+1)
      if (m2.gt.n) go to 6630
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6630
      ii=ia+wty(m2)
      m3=mx(ssj+1)
      m4=mn(ssj+1)
      if (m4.gt.n) go to 6630
      do 701 i=m2,m1
      if (ssj.eq.is) m3=i-1
      if (m3.lt.m4) go to 701
      jj=ja+wtx(i,js+1)+wty(m4)
      cii=c(ii)
      sii=s(ii)
      do 702 j=m4,m3
      lad=arr+ladd(j+aos)
c     z=-val1*(int(lad+tr1)+val2*int(lad+tr2))
c     sii=sii+z*c(jj)
c     s(jj)=s(jj)+z*cii
      z=-val1*(c(ii)*s(jj)+c(jj)*s(ii))
      int(lad+tr1)=int(lad+tr1)+z
      int(lad+tr2)=int(lad+tr2)+z*val2
c
  702 jj=jj+1
c     s(ii)=s(ii)+sii
  701 ii=ii+1
c
c no 136
c
 6630 ssj=xor(is,js)
      if (is.gt.ssj) return
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) return
      m3=mx(is+1)
      m4=mn(is+1)
      if (m4.gt.n) return
      do 9080 i=m2,m1
      if (is.eq.ssj) m3=i-1
      if (m3.lt.m4) go to 9080
      ii=ia+wty(m4)
      jj=ja+wtx(i,js+1)+wty(m4)
      lad=arr+ladd(i+aos)
c     z=val1*(int(lad+tr1)+val2*int(lad+tr2))
      z=0.0d+00
      do 9081 j=m4,m3
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=z+c(ii)*s(jj)+c(jj)*s(ii)
      ii=ii+1
 9081 jj=jj+1
      z=z*val1
      int(lad+tr1)=int(lad+tr1)+z
      int(lad+tr2)=int(lad+tr2)+z*val2
c
 9080 continue
      return
c
c    wz   no 23  iseg=11
c***********************************************************************
  309 if (iseg.eq.10) go to 3091
      if (iseg.eq.14) return
      if (iseg.ne.11) go to 999
      val4=val1*sqrt2
      if (is.ne.0.or.js.ne.0) go to 6631
      jj=ja
      cjj=c(jj)
      sjj=s(jj)
      do 901 i=1,n
      ii=ia+wab(i)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)+ladd(i+los)
c     z=val4*int(lad+1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad+1)=int(lad+1)+val4*(cjj*s(ii)+c(ii)*sjj)
c
  901 continue
c     s(jj)=s(jj)+sjj
c
c     no 25
c
 6631 if (js.ne.0) return
      jj=ja
      cjj=c(jj)
      sjj=s(jj)
      do 902 i=2,n
      ssi=xor(ss(i),is)
      if (ssi.gt.ss(i)) go to 902
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 902
      if (ssi.eq.ss(i)) m1=i-1
      ii=ia+wtw(i,is+1)+wty(m2)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)
      do 903 j=m2,m1
      lad1=lad+ladd(j+los)
c     z=val1*(int(lad1+1)+int(lad1+3))
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      z=val1*(cjj*s(ii)+c(ii)*sjj)
      int(lad1+1)=int(lad1+1)+z
      int(lad1+3)=int(lad1+3)+z
c
  903 ii=ii+1
  902 continue
c     s(jj)=s(jj)+sjj
      return
c
c   yw   no 19   iseg=19
c***********************************************************************
  308 if (iseg.ne.19) go to 999
      val4=val1*sqrt2
      if (js.ne.0) go to 6632
      m1=mx(is+1)
      m2=mn(is+1)
      if (m2.gt.n) go to 6632
      if (m1.gt.n) m1=n
      ii=ia+wty(m2)
      do 801 i=m2,m1
      jj=ja+wab(i)
      lad=arr+ladd(i+aos)
c     z=val4*(int(lad+tr1)+val2*int(lad+tr2))
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=val4*(c(ii)*s(jj)+c(jj)*s(ii))
      int(lad+tr1)=int(lad+tr1)+z
      int(lad+tr2)=int(lad+tr2)+z*val2
c
      ii=ii+1
  801 continue
c
c      no 19a
c
 6632 ssj=xor(is,js)
      if (ssj.gt.is) go to 6633
      m1=mx(is+1)
      m2=mn(is+1)
      if (m2.gt.n) go to 6633
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) go to 6633
      ii=ia+wty(m2)
      m3=mx(ssj+1)
      m4=mn(ssj+1)
      if (m4.gt.n) go to 6633
      do 802 i=m2,m1
      if (ssj.eq.is) m3=i-1
      if (m3.lt.m4) go to 802
      jj=ja+wtw(i,js+1)+wty(m4)
      cii=c(ii)
      sii=s(ii)
      do 803 j=m4,m3
      lad=arr+ladd(j+aos)
c     z=val1*(int(lad+tr1)+val2*int(lad+tr2))
c     sii=sii+z*c(jj)
c     s(jj)=s(jj)+z*cii
      z=val1*(cii*s(jj)+c(jj)*sii)
      int(lad+tr1)=int(lad+tr1)+z
      int(lad+tr2)=int(lad+tr2)+z*val2
c
  803 jj=jj+1
c     s(ii)=s(ii)+sii
  802 ii=ii+1
c
c no 137
c
 6633 ssj=xor(is,js)
      if (is.gt.ssj) return
      m1=mx(ssj+1)
      m2=mn(ssj+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      if (m2.lt.2) m2=2
      if (m1.lt.m2) return
      m3=mx(is+1)
      m4=mn(is+1)
      if (m4.gt.n) return
      do 9082 i=m2,m1
      if (is.eq.ssj) m3=i-1
      if (m3.lt.m4) go to 9082
      ii=ia+wty(m4)
      jj=ja+wtw(i,js+1)+wty(m4)
      lad=arr+ladd(i+aos)
c     z=val1*(int(lad+tr1)+val2*int(lad+tr2))
      z=0.0d+00
      do 9083 j=m4,m3
c     s(jj)=s(jj)+z*c(ii)
c     s(ii)=s(ii)+z*c(jj)
      z=z+c(ii)*s(jj)+s(ii)*c(jj)
c
      ii=ii+1
 9083 jj=jj+1
      z=z*val1
      int(lad+tr1)=int(lad+tr1)+z
      int(lad+tr2)=int(lad+tr2)+z*val2
c
 9082 continue
      return
c
c   wz   no 23 iseg=10
c
 3091 if (is.ne.0.or.js.ne.0) go to 6634
      jj=ja
      cjj=c(jj)
      sjj=s(jj)
      do 904 i=1,n
      ii=ia+wab(i)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)+ladd(i+los)
c     z=val1*int(lad+1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad+1)=int(lad+1)+val1*(cjj*s(ii)+c(ii)*sjj)
c
  904 continue
c     s(jj)=s(jj)+sjj
c
c      no  25
c
 6634 val4=sqrt2*val1
      if (js.ne.0) return
      jj=ja
      cjj=c(jj)
      sjj=s(jj)
      do 905 i=2,n
      ssi=xor(ss(i),is)
      if (ssi.gt.ss(i)) go to 905
      m1=mx(ssi+1)
      m2=mn(ssi+1)
      if (m2.gt.n) go to 905
      if (ssi.eq.ss(i)) m1=i-1
      ii=ia+wtw(i,is+1)+wty(m2)
      lsm=xor(asm,ss(i))
      los=os(lsm+1)
      lad=arr+kadd(i+aos)
      do 906 j=m2,m1
      lad1=lad+ladd(j+los)
c     z=val4*int(lad1+1)
c     s(ii)=s(ii)+z*cjj
c     sjj=sjj+z*c(ii)
      int(lad1+1)=int(lad1+1)+val4*(cjj*s(ii)+c(ii)*sjj)
c
  906 ii=ii+1
  905 continue
c     s(jj)=s(jj)+sjj
      return
c
c    yy    no  9    four externals
c***********************************************************************
c
  999 write (itape6,998) m,iseg
  998 format (2i4)
      return
c
c     closed internal loop
c***********************************************************************
  318 continue
      z=(dot(c(ia),s(ja),n1)+dot(c(ja),s(ia),n1) )*val1
      int(arr+tr1)=int(arr+tr1)+z
      int(arr+tr2)=int(arr+tr2)+z*val2
c         d(ia,ja;idry,jdry)
c         ia .ne. ja  only should occur here
      if( (ksegsv.lt.173).or.(ksegsv.gt.200))return
      if((tr1.ne.3).or.(idry.eq.jdry))then
      write(itape6,40401) idry,jdry,tr1,ksegsv
40401 format( ' idry jdry tr1 ksegsv ',2i6)
      stop ' tr1 .ne 3 or idry eq jdry  in initex '
      else
      iix=indxs(idry)+jdry
      int1(iix)=int1(iix)+val1*dot(c(ja),s(ia),n1)
      iix=indxs(jdry)+idry
      int1(iix)=int1(iix)+val1*dot(c(ia),s(ja),n1)
      endif
      return
c
c***********************************************************************
  314 print 315
  315 format(' incorrect entry to this subroutine. this is for 4x.')
      stop
c
c     ----- internal case for tracks of (3,2,1) -----
c***********************************************************************
  319 continue
      z=(dot(c(ia),s(ja),n1)+dot(c(ja),s(ia),n1) )*val1
      int(arr+1)=int(arr+1)+z
      int(arr+2)=int(arr+2)+z
      int(arr+3)=int(arr+3)+z
c         d(ia,ja;idry,jdry)
c         ia .ne. ja  only should occur here
      if( (ksegsv.lt.173).or.(ksegsv.gt.200))return
      iix=indxs(idry)+jdry
      int1(iix)=int1(iix)+val1*dot(c(ja),s(ia),n1)
      iix=indxs(jdry)+idry
      int1(iix)=int1(iix)+val1*dot(c(ia),s(ja),n1)
      return
c
c***********************************************************************
c
  317 z=val1*(int(arr+tr1)+val2*int(arr+tr2))
c   zz   arrival,for completeness
      stop ' zz arrival in external'
cbl   return
      end
