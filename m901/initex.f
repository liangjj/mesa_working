*deck @(#)initex.f	5.1  11/6/94
      subroutine initex(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      implicit real*8 (a-h,o-z)
c
      integer xor
      integer arr,tr1,tr2,asm,aos,os,wtw,wtx,wty,wab,ss,ssi,ssj,symorb
      integer bmax,orbfrm
      dimension kadd(symorb),ladd(symorb),ijadd(numij),wtw(orbfrm,nsym)
      dimension wtx(orbfrm,nsym),wty(orbfrm),wab(orbfrm),ss(norbs)
      real*8 int(nmax),c(nwks,mxvc),s(nwks,mxvc)
      real*8 z,val1,val2,val3
c
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
c     universal identity of the objects in these common
c     common /all/ val1,val2,val3,arr,itr1,itr2,ia,ja,itype,isegt
c    *,lvfrm1,nlwki,nlwkj,imax,imin
c     common /symq/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
c    #,             numsym(8)
      common /symq/ asm,js,is,mx(8),mn(8),os(8),numsym(8)
      common /io/     itape5,itape6
      common /all/val1,val2,val3,arr,tr1,tr2,ia,ja,m,iseg,n,n1,n2
     *,           imax,imin
      common /minmax/ iming,imaxg,jming,jmaxg
      common /nvectr/ nvc,mxvc
c
      save sqrt2,sqt1p5
c
      sqrt2=sqrt(2.0d+00)
      sqt1p5=sqrt(1.5d+00)
      return
c
c
      entry shapes(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
c
c
      aos=os(asm+1)
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
      m1=mx(is+1)
      m2=mn(is+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      do 8001 ivc=1,nvc
         ii=ia+wty(m2)
         jj=ja
         cjj=c(jj,ivc)*val1
         sjj=0.0
         lad=arr+ladd(m2+aos)
         do 152 i=m2,m1
            z=int(lad+tr1)+val2*int(lad+tr2)
            s(ii,ivc)=s(ii,ivc)+z*cjj
            sjj=sjj+z*c(ii,ivc)
            ii=ii+1
            lad=lad+3
  152    continue
         s(jj,ivc)=s(jj,ivc)+sjj*val1
 8001 continue
      return
c
c    yz   n0 1  iseg=16
c
 3011 continue
      if (js.ne.0) return
      m1=mx(is+1)
      m2=mn(is+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      do 8002 ivc=1,nvc
         ii=ia+wty(m2)
         jj=ja
         cjj=c(jj,ivc)*val1
         sjj=0.0
         lad=arr+ladd(m2+aos)
         do 151 i=m2,m1
            z=int(lad+3)+int(lad+1)
            s(ii,ivc)=s(ii,ivc)+z*cjj
            sjj=sjj+z*c(ii,ivc)
            ii=ii+1
            lad=lad+3
  151    continue
         s(jj,ivc)=s(jj,ivc)+sjj*val1
 8002 continue
      return
c
c  yz  no 1  iseg=18
c
 3012 if (js.ne.0) return
      m1=mx(is+1)
      m2=mn(is+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      do 8003 ivc=1,nvc
         ii=ia+wty(m2)
         jj=ja
         cjj=c(jj,ivc)*val1
         sjj=0.0
         lad=arr+ladd(m2+aos)
         do 150 i=m2,m1
            z=int(lad+3)
            s(ii,ivc)=s(ii,ivc)+z*cjj
            sjj=sjj+z*c(ii,ivc)
            ii=ii+1
            lad=lad+3
  150    continue
         s(jj,ivc)=s(jj,ivc)+sjj*val1
 8003 continue
      return
c
c   xz  no 2 iseg=13
c***********************************************************************
  302 if (iseg.ne.13) go to 999
      if (js.ne.0) return
      do 8004 ivc=1,nvc
         jj=ja
         cjj=c(jj,ivc)*val1
         sjj=0.0
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
            lad=arr+kadd(i+aos)+ladd(m2+los)
            do 251 j=m2,m1
               z=int(lad+1)-int(lad+3)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  251       continue
  250    continue
         s(jj,ivc)=s(jj,ivc)+sjj*val1
 8004 continue
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
      do 8005 ivc=1,nvc
         jj=ja+wty(m2)
         do 365 i=m2,m1
            if (ssi.eq.js) m3=i-1
            if (m3.lt.m4) go to 365
            ii=ia+wtx(i,is+1)+wty(m4)
            cjj=-c(jj,ivc)*val1
            sjj=0.0
            lad=arr+ladd(m4+aos)
            do 366 j=m4,m3
               z=int(lad+3)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  366       continue
            s(jj,ivc)=s(jj,ivc)-sjj*val1
  365    jj=jj+1
 8005 continue
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
         do 8006 ivc=1,nvc
            ii=ia+wtx(i,is+1)+wty(m4)
            jj=ja+wty(m4)
            z=val1*int(arr+ladd(i+aos)+3)
cdir$ ivdep
            do 5353 j=m4,m3
               s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
               s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
               ii=ii+1
               jj=jj+1
 5353       continue
 8006    continue
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
      do 8007 ivc=1,nvc
         jj=ja+wty(m2)
         do 363 i=m2,m1
            cjj=-c(jj,ivc)*val1
            sjj=0.0
            if (ssi.eq.js) m3=i-1
            if (m3.lt.m4) go to 363
            ii=ia+wtx(i,is+1)+wty(m4)
            lad=arr+ladd(m4+aos)
            do 364 j=m4,m3
               z=int(lad+3)+int(lad+1)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  364       continue
            s(jj,ivc)=s(jj,ivc)-sjj*val1
  363    jj=jj+1
 8007 continue
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
         do 8008 ivc=1,nvc
            ii=ia+wtx(i,is+1)+wty(m4)
            jj=ja+wty(m4)
            lad=arr+ladd(i+aos)
            z=val1*(int(lad+3)+int(lad+1))
cdir$ ivdep
            do 6353 j=m4,m3
               s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
               s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
               ii=ii+1
               jj=jj+1
 6353       continue
 8008    continue
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
      do 8009 ivc=1,nvc
         jj=ja+wty(m2)
         do 380 i=m2,m1
            if (ssi.eq.js) m3=i-1
            if (m3.lt.m4) go to 380
            ii=ia+wtx(i,is+1)+wty(m4)
            cjj=-c(jj,ivc)*val1
            sjj=0.0
            lad=arr+ladd(m4+aos)
            do 381 j=m4,m3
               z=int(lad+tr1)+val2*int(lad+tr2)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  381       continue
            s(jj,ivc)=s(jj,ivc)-sjj*val1
  380    jj=jj+1
 8009 continue
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
         do 8010 ivc=1,nvc
            ii=ia+wtx(i,is+1)+wty(m4)
            jj=ja+wty(m4)
            lad=arr+ladd(i+aos)
            z=val1*(int(lad+tr1)+val2*int(lad+tr2))
cdir$ ivdep
            do 7353 j=m4,m3
               s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
               s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
               ii=ii+1
               jj=jj+1
 7353       continue
 8010    continue
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
      do 8011 ivc=1,nvc
         jj=ja+wty(m2)
         do 350 i=m2,m1
            if (ssi.eq.js) m3=i-1
            if (m3.lt.m4) go to 350
            ii=ia+wtx(i,is+1)+wty(m4)
            ksm=xor(asm,ss(i))
            kos=os(ksm+1)
            lad=ijadd(arr+i)+kadd(i+kos)+ladd(m4+los)
            cjj=-c(jj,ivc)*val1
            sjj=0.0
            do 351 j=m4,m3
               z=int(lad+2)-int(lad+1)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  351       continue
            s(jj,ivc)=s(jj,ivc)-sjj*val1
  350    jj=jj+1
 8011 continue
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
         if (js.eq.ssi) m3=i-1
         if (m3.lt.m4) go to 352
         do 8012 ivc=1,nvc
            jj=ja+wty(m4)
            ii=ia+wtx(i,is+1)+wty(m4)
            ksm=xor(asm,ss(i))
            kos=os(ksm+1)
            lad=ijadd(arr+i)
cdir$ ivdep
            do 353 j=m4,m3
               lsm=xor(ksm,ss(j))
               los=os(lsm+1)
               lad1=lad+kadd(j+kos)+ladd(j+los)
               z=-val1*(int(lad1+1)-int(lad1+2))
               s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
               s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
               ii=ii+1
               jj=jj+1
  353       continue
 8012    continue
  352 continue
c
c    no   5
c
 6605 do 354 i=3,n
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
         do 8013 ivc=1,nvc
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
               lad1=lad+kadd(j+kos)+ladd(m4+los)
               cii=-c(ii,ivc)*val1
               sii=0.0
               do 356 k=m4,m3
                  z=int(lad1+1)-int(lad1+2)
                  s(jj,ivc)=s(jj,ivc)+z*cii
                  sii=sii+z*c(jj,ivc)
                  jj=jj+1
                  lad1=lad1+3
  356          continue
               s(ii,ivc)=s(ii,ivc)-sii*val1
  355       ii=ii+1
 8013    continue
  354 continue
c
c    no   6
c
      do 357 i=3,n
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
         do 8014 ivc=1,nvc
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
               lad1=lad+kadd(j+kos)+ladd(m4+los)
               cjj=-c(jj,ivc)*val1
               sjj=0.0
               do 359 k=m4,m3
                  z=int(lad1+3)-int(lad1+2)
                  s(ii,ivc)=s(ii,ivc)+z*cjj
                  sjj=sjj+z*c(ii,ivc)
                  ii=ii+1
                  lad1=lad1+3
  359          continue
               s(jj,ivc)=s(jj,ivc)-sjj*val1
  358       jj=jj+1
 8014    continue
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
      do 8015 ivc=1,nvc
         jj=ja+wty(m2)
         do 360 i=m2,m1
            ksm=xor(asm,ss(i))
            kos=os(ksm+1)
            lad=ijadd(arr+i)
            cjj=c(jj,ivc)*val1
            sjj=0.0
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
               lad1=lad+kadd(j+kos)+ladd(m4+los)
               do 362 k=m4,m3
                  z=int(lad1+1)-int(lad1+3)
                  s(ii,ivc)=s(ii,ivc)+z*cjj
                  sjj=sjj+z*c(ii,ivc)
                  ii=ii+1
                  lad1=lad1+3
  362          continue
  361       continue
            s(jj,ivc)=s(jj,ivc)+sjj*val1
            jj=jj+1
  360    continue
 8015 continue
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
      do 8016 ivc=1,nvc
         jj=ja+wty(m2)
         do 450 i=m2,m1
            if (is.eq.js) m3=i-1
            if (m3.lt.m4) go to 450
            ii=ia+wty(m4)
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+kadd(i+aos)+ladd(m4+los)+1
            cjj=c(jj,ivc)*val4
            sjj=0.0
            do 451 j=m4,m3
               z=int(lad)-2*int(lad+1)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  451       continue
            s(jj,ivc)=s(jj,ivc)+sjj*val4
  450    jj=jj+1
 8016 continue
c
c   no 10
c
 6606 if(ia.eq.ja)return
      if (is.ne.js) return
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) return
      if (m1.gt.n) m1=n
      do 8017 ivc=1,nvc
         ii=ia+wty(m2)
         jj=ja+wty(m2)
cdir$ ivdep
         do 460 i=m2,m1
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+ladd(i+los)+kadd(i+aos)
            z=val4*(int(lad+1)-2.*int(lad+2))
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
            ii=ii+1
            jj=jj+1
  460    continue
 8017 continue
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
      do 8018 ivc=1,nvc
         jj=ja+wty(m2)
         do 420 i=m2,m1
            if (is.eq.js) m3=i-1
            if (m3.lt.m4) go to 420
            ii=ia+wty(m4)
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+kadd(i+aos)+ladd(m4+los)+tr1
            cjj=c(jj,ivc)*val4
            sjj=0.0
            do 421 j=m4,m3
               z=int(lad)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  421       continue
            s(jj,ivc)=s(jj,ivc)+sjj*val4
  420    jj=jj+1
 8018 continue
c
c  no 10
c
 6607 if(ia.eq.ja) return
      if (is.ne.js) go to 6598
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.gt.n) go to 6598
      if (m1.gt.n) m1=n
      do 8019 ivc=1,nvc
         ii=ia+wty(m2)
         jj=ja+wty(m2)
cdir$ ivdep
         do 430 i=m2,m1
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+ladd(i+los)+kadd(i+aos)
            z=val4*int(lad+1)
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
            ii=ii+1
            jj=jj+1
  430    continue
 8019 continue
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
      do 8020 ivc=1,nvc
         ii=ia+wty(m2)
         do 440 i=m2,m1
            if (is.eq.js) m3=i-1
            if (m3.lt.m4) go to 440
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
c
c bug fixed 4 february 1987   pws at lanl
c   in following, the 'm4' was erroneously 'j'
c
            lad=arr+kadd(i+aos)+ladd(m4+los)
            jj=ja+wty(m4)
            cii=c(ii,ivc)*val4
            sii=0.0
            do 441 j=m4,m3
               z=int(lad+1)
               s(jj,ivc)=s(jj,ivc)+z*cii
               sii=sii+z*c(jj,ivc)
               jj=jj+1
               lad=lad+3
  441       continue
            s(ii,ivc)=s(ii,ivc)+sii*val4
  440    ii=ii+1
 8020 continue
c
      return
c
c  yy  no 10  iseg=5
c
 3043 continue
      do 8021 ivc=1,nvc
         m1=mx(js+1)
         m2=mn(js+1)
         if (m2.gt.n) go to 8021
         if (m1.gt.n) m1=n
         ii=ia+wty(m2)
         jj=ja+wty(m2)
         xx=-val1/sqrt2+val3*sqt1p5
         if (abs(xx).lt.1.0d-06) go to 490
         val5=sqrt2*val1
         if(ia.eq.ja) go to 491
         if (is.ne.js) go to 491
cdir$ ivdep
         do 493 i=m2,m1
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+ladd(i+los)+kadd(i+aos)
            z=xx*int(lad+1)+val5*int(lad+2)
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
            ii=ii+1
            jj=jj+1
  493    continue
         go to 491
  490    val5=sqrt2*val1
         if(ia.eq.ja) go to 492
         if (is.ne.js) go to 492
cdir$ ivdep
         do 494 i=m2,m1
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+ladd(i+los)+kadd(i+aos)
            z=val5*int(lad+2)
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
            ii=ii+1
            jj=jj+1
  494    continue
         go to 492
c
c   no   9
c
  491    if (is.gt.js) go to 8021
         m1=mx(js+1)
         m2=mn(js+1)
         if (m2.gt.n) go to 8021
         if (m1.gt.n) m1=n
         if (m2.lt.2) m2=2
         if (m1.lt.m2) go to 8021
         m3=mx(is+1)
         m4=mn(is+1)
         if (m4.gt.n) go to 8021
         jj=ja+wty(m2)
         do 495 i=m2,m1
            if (is.eq.js) m3=i-1
            if (m3.lt.m4) go to 495
            ii=ia+wty(m4)
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+kadd(i+aos)+ladd(m4+los)
            cjj=c(jj,ivc)
            sjj=0.0
            do 496 j=m4,m3
               z=xx*int(lad+1)+val5*int(lad+2)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  496       continue
            s(jj,ivc)=s(jj,ivc)+sjj
  495    jj=jj+1
         go to 8021
c
  492    if (is.gt.js) go to 8021
         m1=mx(js+1)
         m2=mn(js+1)
         if (m2.gt.n) go to 8021
         if (m1.gt.n) m1=n
         if (m2.lt.2) m2=2
         if (m1.lt.m2) go to 8021
         m3=mx(is+1)
         m4=mn(is+1)
         if (m4.gt.n) go to 8021
         jj=ja+wty(m2)
         do 497 i=m2,m1
            if (is.eq.js) m3=i-1
            if (m3.lt.m4) go to 497
            ii=ia+wty(m4)
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+kadd(i+aos)+ladd(m4+los)+2
            cjj=c(jj,ivc)*val5
            sjj=0.0
            do 498 j=m4,m3
               z=int(lad)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  498       continue
            s(jj,ivc)=s(jj,ivc)+sjj*val5
  497    jj=jj+1
 8021 continue
      return
c
c   yy   no 8  iseg=6
c
 3041 continue
      do 8022 ivc=1,nvc
         xx=-val1/sqrt2+val3*sqt1p5
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
            lad=arr+kadd(i+aos)+ladd(m4+los)
            cii=c(ii,ivc)
            sii=0.0
            do 473 j=m4,m3
               z=xx*int(lad+1)+val5*int(lad+2)
               sii=sii+z*c(jj,ivc)
               s(jj,ivc)=s(jj,ivc)+z*cii
               jj=jj+1
               lad=lad+3
  473       continue
            s(ii,ivc)=s(ii,ivc)+sii
  472    ii=ii+1
         go to 481
c
  480    val5=sqrt2*val1
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
            lad=arr+kadd(i+aos)+ladd(m4+los)
            cii=c(ii,ivc)*val5
            sii=0.0
            do 475 j=m4,m3
               z=int(lad+2)
               sii=sii+z*c(jj,ivc)
               s(jj,ivc)=s(jj,ivc)+z*cii
               jj=jj+1
               lad=lad+3
  475       continue
            s(ii,ivc)=s(ii,ivc)+sii*val5
  474    ii=ii+1
         go to 483
c
c    no  10
c
  481    if (is.ne.js) go to 482
         if(ia.eq.ja) go to 482
         m1=mx(js+1)
         m2=mn(js+1)
         if (m2.gt.n) go to 482
         if (m1.gt.n) m1=n
         ii=ia+wty(m2)
         jj=ja+wty(m2)
cdir$ ivdep
         do 470 i=m2,m1
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+ladd(i+los)+kadd(i+aos)
            z=xx*int(lad+1)+val5*int(lad+2)
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
            ii=ii+1
            jj=jj+1
  470    continue
         go to 482
c
  483    if (is.ne.js) go to 484
         if (ia.eq.ja) go to 484
         m1=mx(js+1)
         m2=mn(js+1)
         if (m2.gt.n) go to 484
         if (m1.gt.n) m1=n
         ii=ia+wty(m2)
         jj=ja+wty(m2)
cdir$ ivdep
         do 471 i=m2,m1
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+ladd(i+los)+kadd(i+aos)
            z=val5*int(lad+2)
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
            ii=ii+1
            jj=jj+1
  471    continue
         go to 484
c
c  no  9
c
  482    if (is.gt.js) go to 8022
         m1=mx(js+1)
         m2=mn(js+1)
         if (m2.gt.n) go to 8022
         if (m1.gt.n) m1=n
         if (m2.lt.2) m2=2
         if (m1.lt.m2) go to 8022
         m3=mx(is+1)
         m4=mn(is+1)
         if (m4.gt.n) go to 8022
         jj=ja+wty(m2)
         do 476 i=m2,m1
            if (is.eq.js) m3=i-1
            if (m3.lt.m4) go to 476
            ii=ia+wty(m4)
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+kadd(i+aos)+ladd(m4+los)
            cjj=c(jj,ivc)
            sjj=0.0
            do 477 j=m4,m3
               z=xx*int(lad+3)+val5*int(lad+2)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  477       continue
            s(jj,ivc)=s(jj,ivc)+sjj
  476    jj=jj+1
         go to 8022
c
  484    if (is.gt.js) go to 8022
         m1=mx(js+1)
         m2=mn(js+1)
         if (m2.gt.n) go to 8022
         if (m1.gt.n) m1=n
         if (m2.lt.2) m2=2
         if (m1.lt.m2) go to 8022
         m3=mx(is+1)
         m4=mn(is+1)
         if (m4.gt.n) go to 8022
         jj=ja+wty(m2)
         do 478 i=m2,m1
            if (is.eq.js) m3=i-1
            if (m3.lt.m4) go to 478
            ii=ia+wty(m4)
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+kadd(i+aos)+ladd(m4+los)
            cjj=c(jj,ivc)*val5
            sjj=0.0
            do 479 j=m4,m3
               z=int(lad+2)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  479       continue
            s(jj,ivc)=s(jj,ivc)+sjj*val5
  478    jj=jj+1
 8022 continue
c
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
      do 8023 ivc=1,nvc
         val4=val1*sqrt2
         if (is.ne.0) go to 6587
         do 1001 i=2,n
            if (js.gt.ss(i)) go to 1001
            ii=ia+wab(i)
            m1=mx(js+1)
            m2=mn(js+1)
            if (m2.gt.n) go to 1001
            if (js.eq.ss(i).or.m1.gt.n) m1=i-1
            jj=ja+wty(m2)
            ksm=xor(asm,ss(i))
            kos=os(ksm+1)
            lad=ijadd(arr+i)+kadd(i+kos)+ladd(m2+aos)
            cii=c(ii,ivc)*val4
            sii=0.0
            do 1002 j=m2,m1
               z=int(lad+1)
               sii=sii+z*c(jj,ivc)
               s(jj,ivc)=s(jj,ivc)+z*cii
               jj=jj+1
               lad=lad+3
 1002       continue
            s(ii,ivc)=s(ii,ivc)+sii*val4
 1001    continue
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
            ksm=xor(asm,ss(i))
            kos=os(ksm+1)
            lad=ijadd(arr+i)
            cjj=c(jj,ivc)*val4
            sjj=0.0
cdir$ ivdep
            do 1004 j=1,i-1
               ii=ia+wab(j)
               lsm=xor(ksm,ss(j))
               los=os(lsm+1)
               lad1=lad+kadd(j+kos)+ladd(j+los)
               z=int(lad1+1)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
 1004       continue
            s(jj,ivc)=s(jj,ivc)+sjj*val4
            jj=jj+1
 1003    continue
c
c   no 30
c
 6587    m1=mx(js+1)
         m2=mn(js+1)
         if (m2.gt.n) go to 6599
         if (m1.gt.n) m1=n
         if (m2.lt.3) m2=3
         if (m1.lt.m2) go to 6599
         jj=ja+wty(m2)
         do 1005 i=m2,m1
            ksm=xor(asm,ss(i))
            kos=os(ksm+1)
            lad=ijadd(arr+i)
            cjj=c(jj,ivc)*val1
            sjj=0.0
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
               lad1=lad+kadd(j+kos)+ladd(m4+los)
               do 1007 k=m4,m3
                  z=int(lad1+1)+int(lad1+3)
                  s(ii,ivc)=s(ii,ivc)+z*cjj
                  sjj=sjj+z*c(ii,ivc)
                  ii=ii+1
                  lad1=lad1+3
 1007          continue
 1006       continue
            s(jj,ivc)=s(jj,ivc)+sjj*val1
            jj=jj+1
 1005    continue
c
c   no 32
c
 6599    if(n.lt.3) go to 1018
         do 1008 i=3,n
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
               lad1=lad+kadd(j+kos)+ladd(m4+los)
               cii=c(ii,ivc)*val1
               sii=0.0
               do 1010 k=m4,m3
                  z=int(lad1+1)+int(lad1+2)
                  sii=sii+z*c(jj,ivc)
                  s(jj,ivc)=s(jj,ivc)+z*cii
                  jj=jj+1
                  lad1=lad1+3
 1010          continue
               s(ii,ivc)=s(ii,ivc)+sii*val1
 1009          ii=ii+1
 1008    continue
 1018    continue
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
            if (ssi.eq.js) m3=i-1
            if (m3.lt.m4) go to 1011
            ii=ia+wtw(i,is+1)+wty(m4)
            ksm=xor(asm,ss(i))
            kos=os(ksm+1)
            lad=ijadd(arr+i)+kadd(i+kos)+ladd(m4+los)
            cjj=c(jj,ivc)*val1
            sjj=0.0
               do 1012 j=m4,m3
               z=int(lad+2)+int(lad+1)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
 1012       continue
            s(jj,ivc)=s(jj,ivc)+sjj*val1
 1011    jj=jj+1
c
c   no 36
c
 6608    ssi=xor(js,is)
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
            if (js.eq.ssi) m3=i-1
            if (m3.lt.m4) go to 1013
            jj=ja+wty(m4)
            ii=ia+wtw(i,is+1)+wty(m4)
            ksm=xor(asm,ss(i))
            kos=os(ksm+1)
            lad=ijadd(arr+i)
cdir$ ivdep
            do 1014 j=m4,m3
               lsm=xor(ksm,ss(j))
               los=os(lsm+1)
               lad1=lad+kadd(j+kos)+ladd(j+los)
               z=val1*(int(lad1+1)+int(lad1+2))
               s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
               s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
               ii=ii+1
               jj=jj+1
 1014       continue
 1013    continue
c
c    no 38
c
 6609    do 1015 i=3,n
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
               lad1=lad+kadd(j+kos)+ladd(m4+los)
               cjj=c(jj,ivc)*val1
               sjj=0.0
               do 1017 k=m4,m3
                  z=int(lad1+3)+int(lad1+2)
                  s(ii,ivc)=s(ii,ivc)+z*cjj
                  sjj=sjj+z*c(ii,ivc)
                  ii=ii+1
                  lad1=lad1+3
 1017          continue
               s(jj,ivc)=s(jj,ivc)+sjj*val1
 1016       jj=jj+1
 1015    continue
 8023 continue
c
      return
c
c   wy   no 26   iseg=21
c
 3103 continue
      do 8024 ivc=1,nvc
         val4=val1*sqrt2
         if (is.ne.0) go to 6610
         m1=mx(js+1)
         m2=mn(js+1)
         if (m2.gt.n) go to 6610
         if (m1.gt.n) m1=n
         jj=ja+wty(m2)
cdir$ ivdep
         do 1020 i=m2,m1
            ii=ia+wab(i)
            lad=arr+ladd(i+aos)
            z=val4*(int(lad+tr1)+val2*int(lad+tr2))
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
            jj=jj+1
 1020    continue
c
c   no 33
c
 6610    ssi=xor(js,is)
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
            cjj=c(jj,ivc)*val1
            sjj=0.0
            lad=arr+ladd(m4+aos)
            do 1022 j=m4,m3
               z=int(lad+tr1)+val2*int(lad+tr2)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
 1022       continue
            s(jj,ivc)=s(jj,ivc)+sjj*val1
 1021    jj=jj+1
c
c no 149
c
 6611    ssi=xor(js,is)
         if (js.gt.ssi) go to 8024
         m1=mx(ssi+1)
         m2=mn(ssi+1)
         if (m2.gt.n) go to 8024
         if (m1.gt.n) m1=n
         if (m2.lt.2) m2=2
         if (m1.lt.m2) go to 8024
         m3=mx(js+1)
         m4=mn(js+1)
         if (m4.gt.n) go to 8024
         do 6013 i=m2,m1
            if (js.eq.ssi) m3=i-1
            if (m3.lt.m4) go to 6013
            ii=ia+wtw(i,is+1)+wty(m4)
            jj=ja+wty(m4)
            lad=arr+ladd(i+aos)
            z=val1*(int(lad+tr1)+val2*int(lad+tr2))
cdir$ ivdep
            do 6014 j=m4,m3
               s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
               s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
               ii=ii+1
               jj=jj+1
 6014       continue
 6013    continue
 8024 continue
c
      return
c
c    wy  no 26   iseg=15
c
 3101 continue
      do 8025 ivc=1,nvc
         val4=val1*sqrt2
         if (is.ne.0) go to 6612
         m1=mx(js+1)
         m2=mn(js+1)
         if (m2.gt.n) go to 6612
         if (m1.gt.n) m1=n
         jj=ja+wty(m2)
cdir$ ivdep
         do 1030 i=m2,m1
            ii=ia+wab(i)
            lad=arr+ladd(i+aos)
            z=val4*(int(lad+1)+int(lad+2)+int(lad+3))
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
            jj=jj+1
 1030    continue
c
c    no 33
c
 6612    ssi=xor(js,is)
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
            cjj=c(jj,ivc)*val1
            sjj=0.0
            lad=arr+ladd(m4+aos)
            do 1032 j=m4,m3
               z=int(lad+3)+int(lad+1)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
 1032       continue
            s(jj,ivc)=s(jj,ivc)+sjj*val1
 1031    jj=jj+1
c
c no 149
c
 6613    ssi=xor(js,is)
         if (js.gt.ssi) go to 8025
         m1=mx(ssi+1)
         m2=mn(ssi+1)
         if (m2.gt.n)go to 8025
         m3=mx(js+1)
         m4=mn(js+1)
         if (m4.gt.n) go to 8025
         if (m1.gt.n) m1=n
         if (m2.lt.2) m2=2
         if (m1.lt.m2) go to 8025
         do 7013 i=m2,m1
            if (js.eq.ssi) m3=i-1
            if (m3.lt.m4) go to 7013
            ii=ia+wtw(i,is+1)+wty(m4)
            jj=ja+wty(m4)
            lad=arr+ladd(i+aos)
            z=val1*(int(lad+3)+int(lad+1))
cdir$ ivdep
            do 7014 j=m4,m3
               s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
               s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
               ii=ii+1
               jj=jj+1
 7014       continue
 7013    continue
 8025 continue
      return
c
c   wy  no 26 iseg=17
c
 3102 continue
      do 8026 ivc=1,nvc
         val4=val1*sqrt2
         if (is.ne.0) go to 6614
         m1=mx(js+1)
         m2=mn(js+1)
         if (m2.gt.n) go to 6614
         if (m1.gt.n) m1=n
         jj=ja+wty(m2)
         lad=arr+ladd(m2+aos)
cdir$ ivdep
         do 1023 i=m2,m1
            ii=ia+wab(i)
            z=val4*(int(lad+3)+int(lad+2))
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
            jj=jj+1
            lad=lad+3
 1023    continue
c
c   no  33
c
 6614    ssi=xor(js,is)
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
            cjj=c(jj,ivc)*val1
            sjj=0.0
            lad=arr+ladd(m4+aos)
            do 1025 j=m4,m3
               z=int(lad+3)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
 1025       continue
            s(jj,ivc)=s(jj,ivc)+sjj*val1
 1024    jj=jj+1
c
c no 149
c
 6615    ssi=xor(js,is)
         if (js.gt.ssi) go to 8026
         m1=mx(ssi+1)
         m2=mn(ssi+1)
         if (m2.gt.n) go to 8026
         m3=mx(js+1)
         m4=mn(js+1)
         if (m4.gt.n) go to 8026
         if (m1.gt.n) m1=n
         if (m2.lt.2) m2=2
         if (m1.lt.m2) go to 8026
         do 5013 i=m2,m1
            if (js.eq.ssi) m3=i-1
            if (m3.lt.m4) go to 5013
            ii=ia+wtw(i,is+1)+wty(m4)
            jj=ja+wty(m4)
            lad=arr+ladd(i+aos)
            z=val1*int(lad+3)
cdir$ ivdep
            do 5014 j=m4,m3
               s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
               s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
               ii=ii+1
               jj=jj+1
 5014       continue
 5013    continue
 8026 continue
c
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
      do 8027 ivc=1,nvc
         jj=ja+wty(m2)
         ii=ia
         cii=c(ii,ivc)*val1
         sii=0.0
         lad=arr+ladd(m2+aos)
            do 601 i=m2,m1
            z=int(lad+tr1)+val2*int(lad+tr2)
            sii=sii+z*c(jj,ivc)
            s(jj,ivc)=s(jj,ivc)+z*cii
            jj=jj+1
            lad=lad+3
  601    continue
         s(ii,ivc)=s(ii,ivc)+sii*val1
 8027 continue
c
      return
c
c     yx   no  18   iseg=20
c***********************************************************************
c
  307 continue
      do 8028 ivc=1,nvc
         if (iseg.ne.20) go to 999
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
            cii=-c(ii,ivc)*val1
            sii=0.0
            lad=arr+ladd(m4+aos)
            do 702 j=m4,m3
               z=int(lad+tr1)+val2*int(lad+tr2)
               sii=sii+z*c(jj,ivc)
               s(jj,ivc)=s(jj,ivc)+z*cii
               jj=jj+1
               lad=lad+3
  702       continue
            s(ii,ivc)=s(ii,ivc)-sii*val1
  701    ii=ii+1
c
c no 136
c
 6630    ssj=xor(is,js)
         if (is.gt.ssj) go to 8028
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.gt.n) go to 8028
         if (m1.gt.n) m1=n
         if (m2.lt.2) m2=2
         if (m1.lt.m2) go to 8028
         m3=mx(is+1)
         m4=mn(is+1)
         if (m4.gt.n) go to 8028
         do 9080 i=m2,m1
            if (is.eq.ssj) m3=i-1
            if (m3.lt.m4) go to 9080
            ii=ia+wty(m4)
            jj=ja+wtx(i,js+1)+wty(m4)
            lad=arr+ladd(i+aos)
            z=val1*(int(lad+tr1)+val2*int(lad+tr2))
cdir$ ivdep
            do 9081 j=m4,m3
               s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
               s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
               ii=ii+1
               jj=jj+1
 9081       continue
 9080    continue
 8028 continue
c
      return
c
c    wz   no 23  iseg=11
c***********************************************************************
  309 if (iseg.eq.10) go to 3091
      if (iseg.eq.14) return
      if (iseg.ne.11) go to 999
c
      do 8029 ivc=1,nvc
         val4=val1*sqrt2
         if (is.ne.0.or.js.ne.0) go to 6631
         jj=ja
         cjj=c(jj,ivc)*val4
         sjj=0.0
cdir$ ivdep
         do 901 i=1,n
            ii=ia+wab(i)
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+kadd(i+aos)+ladd(i+los)
            z=int(lad+1)
            s(ii,ivc)=s(ii,ivc)+z*cjj
            sjj=sjj+z*c(ii,ivc)
  901    continue
         s(jj,ivc)=s(jj,ivc)+sjj*val4
c
c     no 25
c
 6631    if (js.ne.0) go to 8029
         jj=ja
         cjj=c(jj,ivc)*val1
         sjj=0.0
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
            lad=arr+kadd(i+aos)+ladd(m2+los)
            do 903 j=m2,m1
               z=int(lad+1)+int(lad+3)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  903       continue
  902    continue
         s(jj,ivc)=s(jj,ivc)+sjj*val1
 8029 continue
c
      return
c
c   yw   no 19   iseg=19
c***********************************************************************
  308 if (iseg.ne.19) go to 999
      do 8030 ivc=1,nvc
         val4=val1*sqrt2
         if (js.ne.0) go to 6632
         m1=mx(is+1)
         m2=mn(is+1)
         if (m2.gt.n) go to 6632
         if (m1.gt.n) m1=n
         ii=ia+wty(m2)
cdir$ ivdep
         do 801 i=m2,m1
            jj=ja+wab(i)
            lad=arr+ladd(i+aos)
            z=val4*(int(lad+tr1)+val2*int(lad+tr2))
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
            ii=ii+1
  801    continue
c
c      no 19a
c
 6632    ssj=xor(is,js)
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
            cii=c(ii,ivc)*val1
            sii=0.0
            lad=arr+ladd(m4+aos)
            do 803 j=m4,m3
               z=int(lad+tr1)+val2*int(lad+tr2)
               sii=sii+z*c(jj,ivc)
               s(jj,ivc)=s(jj,ivc)+z*cii
               jj=jj+1
               lad=lad+3
  803       continue
            s(ii,ivc)=s(ii,ivc)+sii*val1
  802    ii=ii+1
c
c no 137
c
 6633    ssj=xor(is,js)
         if (is.gt.ssj) go to 8030
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.gt.n) go to 8030
         if (m1.gt.n) m1=n
         if (m2.lt.2) m2=2
         if (m1.lt.m2) go to 8030
         m3=mx(is+1)
         m4=mn(is+1)
         if (m4.gt.n) go to 8030
         do 9082 i=m2,m1
            if (is.eq.ssj) m3=i-1
            if (m3.lt.m4) go to 9082
            ii=ia+wty(m4)
            jj=ja+wtw(i,js+1)+wty(m4)
            lad=arr+ladd(i+aos)
            z=val1*(int(lad+tr1)+val2*int(lad+tr2))
cdir$ ivdep
            do 9083 j=m4,m3
               s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
               s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
               ii=ii+1
               jj=jj+1
 9083       continue
 9082    continue
 8030 continue
c
      return
c
c   wz   no 23 iseg=10
c
 3091 if (is.ne.0.or.js.ne.0) go to 6634
      do 8031 ivc=1,nvc
         jj=ja
         cjj=c(jj,ivc)
         sjj=0.0
cdir$ ivdep
         do 904 i=1,n
            ii=ia+wab(i)
            lsm=xor(asm,ss(i))
            los=os(lsm+1)
            lad=arr+kadd(i+aos)+ladd(i+los)
            z=val1*int(lad+1)
            s(ii,ivc)=s(ii,ivc)+z*cjj
            sjj=sjj+z*c(ii,ivc)
  904    continue
         s(jj,ivc)=s(jj,ivc)+sjj
 8031 continue
c
c      no  25
c
 6634 val4=sqrt2*val1
      if (js.ne.0) return
      do 8032 ivc=1,nvc
         jj=ja
         cjj=c(jj,ivc)*val4
         sjj=0.0
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
            lad=arr+kadd(i+aos)+ladd(m2+los)
            do 906 j=m2,m1
               z=int(lad+1)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
  906       continue
  905    continue
         s(jj,ivc)=s(jj,ivc)+sjj*val4
 8032 continue
c
      return
c
c    yy    no  9    four externals
c***********************************************************************
c
  999 write (6,998) m,iseg
  998 format (2i4)
      return
c
c     closed internal loop
c***********************************************************************
  318 z=val1*(int(arr+tr1)+val2*int(arr+tr2))
      ii=ia
      jj=ja
      do 4000 i=1,n1
cdir$ ivdep
         do 8033 ivc=1,nvc
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
 8033    continue
         ii=ii+1
 4000 jj=jj+1
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
      z=val1*(int(arr+1)+int(arr+2)+int(arr+3))
      ii=ia
      jj=ja
      do 3999 i=1,n1
cdir$ ivdep
         do 8034 ivc=1,nvc
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
 8034    continue
         ii=ii+1
 3999 jj=jj+1
      return
c
c***********************************************************************
c
  317 z=val1*(int(arr+tr1)+val2*int(arr+tr2))
c   zz   arrival,for completeness
      stop ' zz arrival in external'
      end
