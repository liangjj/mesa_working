*deck @(#)initst.f	5.1  11/6/94
      subroutine initst(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c
c
c
      implicit real*8 (a-h,o-z)
c
      integer xor
      integer arr,tr1,tr2,asm,os,wtw,wtx,wty,wab,ss,symorb
      integer bmax,orbfrm
      real*8 int(nmax),c(nwks),s(nwks)
c
      dimension ci(1),cj(1),si(1),sj(1),h1(1),h2(1)
      dimension kadd(symorb),ladd(symorb),ijadd(numij),wtw(orbfrm,nsym)
      dimension wtx(orbfrm,nsym),wty(orbfrm),wab(orbfrm),ss(norbs)
      dimension ijxx(numij ),nklxx(nsym),ijww(numij )
      dimension nklww(nsym)
c
      common /coefs/  a,b,intal,intau,intb,intad,intbd
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
      common /count/  icount,ixx4,ixx5,ixx6,ixx8,iww4,iww5,iww6,iww8
     #,               iwx7,ixw9,ixy3,ixy16,ixy18,ixy22,iwy2,iwy15,iwy17
     #,               iwy21,iyx20,iyw19
      common /xx4xbl/   nijvir
      common /minmax/ iming,imaxg,jming,jmaxg
      real*8 val1,val2,val3
c
      save sqrt2,sqt1p5,asqrt2,hsqrt3
c
      sqrt2=sqrt(2.0d+00)
      sqt1p5=sqrt(1.5d+00)
      asqrt2=1.0d+00/sqrt2
      hsqrt3=0.5d+00*sqrt(3.0d+00)
      return
c
c***************************** xx 4x ***********************************
c
      entry xx4x(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $           ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c..rlm     entry xx4x
c
      do 402 i=imin,imax
         ip=i*(i-1)/2
         jsy=xor(is,ss(i))
         jmin=mn(jsy+1)
         if (jmin.ge.i) go to 401
         nkl=nklxx(is+1)
         jmax=mx(jsy+1)
         if (jmax.ge.i) jmax=i-1
         iia=ia+wtx(i,is+1)
         do 400 j=jmin,jmax
            if (j.lt.jming.or.j.gt.jmaxg) go to 400
c           call dot(t,int(ijxx(ip+j)+1),c(ia),nkl)
            s(iia)=s(iia)+dot(int(ijxx(ip+j)+1),c(ia),nkl)
  400    iia=iia+1
  401    continue
  402 continue
      return
c
c***************************** ww 4x ***********************************
c
      entry ww4x(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $           ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c..rlm     entry ww4x
c
      do 414 i=imin,imax
      ip=i*(i-1)/2
      jsy=xor(is,ss(i))
      jmin=mn(jsy+1)
      if (jmin.gt.i) go to 413
      nkl=nklww(is+1)
      jmax=mx(jsy+1)
      if (jmax.gt.i) jmax=i
      if (i.ne.jmin) go to 410
      iia=ia+wab(i)
      go to 411
  410 continue
      iia=ia+wtw(i,is+1)
  411 continue
      do 412 j=jmin,jmax
      if (j.lt.jming.or.j.gt.jmaxg) go to 412
c     call dot(t,int(ijww(ip+j)+1),c(ia),nkl)
      s(iia)=s(iia)+dot(int(ijww(ip+j)+1),c(ia),nkl)
  412 iia=iia+1
  413 continue
  414 continue
      return
c
c******************************** yx ***********************************
c
c
      entry yx(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $         ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c..rlm     entry yx
c
      iyx20=iyx20+1
      a=val1
      b=val1*val2
c
c     ----- dm specific mod -----
c
c     a=a+a
c     b=b+b
c
      intal=tr1
      intb=tr2
      if (is.ne.0) go to 710
      if (js.ne.0) go to 702
c
c     ----- is.eq.0, js.eq.0 -----
c
      numi=numsym(1)
      if (numi.le.1) return
c     call hi    (h1,0,numi,int,ladd)
c     call square(c(ja),ci,numi)
c     call apbc  (s(ia),h1,ci,1,numi,numi)
c     call ebc   (si,h1,c(ia),numi,1,numi)
c     call fold  (s(ja),si,numi)
      call square(s(ja),sj,numi)
      call ebc   (h1,sj,c(ia),numi,numi,1)
c         for tden-dry
      call square(c(ja),cj,numi)
      call apbc   (h1,cj,s(ia),numi,numi,1)
c
      call hi    (h1,0,numi,int,ladd)
c
      return
c
c     ----- is.eq.0, js.ne.0 -----
c
  702 continue
      numi=numsym(js+1)
      if (numi.lt.1) return
      numj=numsym(1)
      if (numj.lt.1) return
      jas=ja+wtx(mn(js+1),js+1)
c     call hi    (h1,js,numi,int,ladd)
c     call apbct (s(ia),h1,c(jas),1,numi,numj)
c     call apbtct(s(jas),c(ia),h1,numj,1,numi)
      call ebtc  (h1,s(jas),c(ia),numi,numj,1)
c         for tden-dry
      call apbtc  (h1,c(jas),s(ia),numi,numj,1)
c
      call hi    (h1,js,numi,int,ladd)
c
      return
c
c     ----- is.ne.0 -----
c
  710 continue
      if (js.ne.0) go to 711
c
c     ----- is.ne.0, js.eq.0 -----
c
      numi=numsym(is+1)
      if (numi.le.1) return
      jas=ja+wtx(mn(is+1)+1,js+1)
c     call hi    (h1,is,numi,int,ladd)
c     call square(c(jas),ci,numi)
c     call apbc  (s(ia),h1,ci,1,numi,numi)
c     call ebc   (si,h1,c(ia),numi,1,numi)
c     call fold  (s(jas),si,numi)
      call square(s(jas),sj,numi)
      call ebc   (h1,sj,c(ia),numi,numi,1)
c         for tden -- dry
      call square(c(jas),cj,numi)
      call apbc   (h1,cj,s(ia),numi,numi,1)
c
      call hi    (h1,is,numi,int,ladd)
c
      return
c
c     ----- is.ne.0, js.ne.0 -----
c
  711 continue
      isym=xor(is,js)
      numj=numsym(is+1)
      if (numj.lt.1) return
      numi=numsym(isym+1)
      if (numi.lt.1) return
      if (isym.lt.is) go to 712
      jas=ja+wtx(mn(isym+1),js+1)
c     call hi    (h1,isym,numi,int,ladd)
c     call apbct (s(ia),h1,c(jas),1,numi,numj)
c     call apbtct(s(jas),c(ia),h1,numj,1,numi)
      call ebtc  (h1,s(jas),c(ia),numi,numj,1)
c         for tden-dry
      call apbtc  (h1,c(jas),s(ia),numi,numj,1)
c
      call hi    (h1,isym,numi,int,ladd)
c
      return
  712 continue
      jas=ja+wtx(mn(is+1),js+1)
c     call hi    (h1,isym,numi,int,ladd)
c     call ambc  (s(ia),h1,c(jas),1,numi,numj)
c     call ambc  (s(jas),h1,c(ia),numi,1,numj)
      call embc  (h1,s(jas),c(ia),numi,numj,1)
c         for tden-dry
      call ambc  (h1,c(jas),s(ia),numi,numj,1)
c
      call hi    (h1,isym,numi,int,ladd)
c
      return
c
c****************************** xy *************************************
c
      entry xy(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $           ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c..rlm     entry xy
      if (iseg.eq.16) go to 500
      if (iseg.eq.18) go to 530
      if (iseg.eq.22) go to 570
c
c     ----- xy iseg=3 -----
c
      ixy3=ixy3+1
      jja=ja
      nkl=nklxx(is+1)
      do 420 i=mn(js+1),mx(js+1)
      if (i.lt.jming.or.i.gt.jmaxg) go to 420
      ij=ijxx(arr+i)
c     call dot(t,int(ij+1),c(ia),nkl)
      s(jja)=s(jja)+val1*dot(int(ij+1),c(ia),nkl)
  420 jja=jja+1
c
      iia=ia
      jinf=mn(js+1)
      if (jming.gt.jinf) jinf=jming
      jdiff=jinf-mn(js+1)
      jsup=mx(js+1)
      if (jmaxg.lt.jsup) jsup=jmaxg
      ist=ijxx(arr+jinf)
      do 422 i=1,nkl
      jja=ja+jdiff
      ipt=ist+i
      t=0.0
      do 421 j=jinf,jsup
      t=t+int(ipt)*c(jja)
      ipt=ipt-nkl
  421 jja=jja+1
      s(iia)=s(iia)+t*val1
      iia=iia+1
  422 continue
      return
c
c     ----- xy iseg=16 -----
c
  500 continue
      ixy16=ixy16+1
      a=val1
      b=val1
      intal=1
      intb=3
      go to 520
c
c     ----- xy iseg=18 -----
c
  530 continue
      ixy18=ixy18+1
      a=val1
      b=0.0
      intal=3
      intb=3
      go to 520
c
c     ----- xy iseg=22 -----
c
  570 continue
      ixy22=ixy22+1
      a=val1
      b=val1*val2
      intal=tr1
      intb=tr2
  520 continue
c
c     ----- dm specific mod -----
c
c     a=a+a
c     b=b+b
c
      if (js.ne.0) go to 510
      if (is.ne.0) go to 502
c
c     ----- is.eq.0, js.eq.0 -----
c
      numi=numsym(1)
      if (numi.le.1) return
c     call hi    (h1,0,numi,int,ladd)
c     call square(c(ia),ci,numi)
c     call apbc  (s(ja),h1,ci,1,numi,numi)
c     call ebc   (si,h1,c(ja),numi,1,numi)
c     call fold  (s(ia),si,numi)
      call square(c(ia),ci,numi)
      call ebc   (h1,ci,s(ja),numi,numi,1)
c         for tden-dry
      call square(s(ia),si,numi)
      call apbc   (h1,si,c(ja),numi,numi,1)
c
      call hi    (h1,0,numi,int,ladd)
c
      return
c
c     ----- is.ne.0, js.eq.0 -----
c
  502 continue
      numi=numsym(is+1)
      if (numi.lt.1) return
      numj=numsym(1)
      if (numj.lt.1) return
      ias=ia+wtx(mn(is+1),is+1)
c     call hi    (h1,is,numi,int,ladd)
c     call apbct (s(ja),h1,c(ias),1,numi,numj)
c     call apbtct(s(ias),c(ja),h1,numj,1,numi)
      call ebtc  (h1,c(ias),s(ja),numi,numj,1)
c         for tden-dry
      call apbtc  (h1,s(ias),c(ja),numi,numj,1)
c
      call hi    (h1,is,numi,int,ladd)
c
      return
c
c     ----- js.ne.0 -----
c
  510 continue
      if (is.ne.0) go to 511
c
c     ----- is.eq.0, js.ne.0 -----
c
      numi=numsym(js+1)
      if (numi.le.1) return
      ias=ia+wtx(mn(js+1)+1,is+1)
c     call hi    (h1,js,numi,int,ladd)
c     call square(c(ias),ci,numi)
c     call apbc  (s(ja),h1,ci,1,numi,numi)
c     call ebc   (si,h1,c(ja),numi,1,numi)
c     call fold  (s(ias),si,numi)
      call square(c(ias),ci,numi)
      call ebc   (h1,ci,s(ja),numi,numi,1)
c         for tden-dry
      call square(s(ias),si,numi)
      call apbc   (h1,si,c(ja),numi,numi,1)
c
      call hi    (h1,js,numi,int,ladd)
c
      return
c
c     ----- is.ne.0, js.ne.0 -----
c
  511 continue
      isym=xor(is,js)
      numj=numsym(js+1)
      if (numj.lt.1) return
      numi=numsym(isym+1)
      if (numi.lt.1) return
      if (isym.lt.js) go to 512
      ias=ia+wtx(mn(isym+1),is+1)
c     call hi    (h1,isym,numi,int,ladd)
c     call apbct (s(ja),h1,c(ias),1,numi,numj)
c     call apbtct(s(ias),c(ja),h1,numj,1,numi)
      call ebtc  (h1,c(ias),s(ja),numi,numj,1)
c         for tden
      call apbtc  (h1,s(ias),c(ja),numi,numj,1)
c
      call hi    (h1,isym,numi,int,ladd)
c
      return
  512 continue
      ias=ia+wtx(mn(js+1),is+1)
c     call hi    (h1,isym,numi,int,ladd)
c     call ambc  (s(ja),h1,c(ias),1,numi,numj)
c     call ambc  (s(ias),h1,c(ja),numi,1,numj)
      call embc  (h1,c(ias),s(ja),numi,numj,1)
c         for tden
      call ambc  (h1,s(ias),c(ja),numi,numj,1)
c
      call hi    (h1,isym,numi,int,ladd)
c
      return
c
c******************************** yw ***********************************
c
      entry yw(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $           ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c..rlm     entry yw
c
      iyw19=iyw19+1
      a=val1
      b=val1*val2
c
c     ----- dm specific mod -----
c
c     a=a+a
c     b=b+b
c
      intal=tr1
      intb=tr2
      if (is.ne.0) go to 810
      if (js.ne.0) go to 802
c
c     ----- is.eq.0, js.eq.0 -----
c
      numi=numsym(1)
      if (numi.lt.1) return
c     call hi    (h1,0,numi,int,ladd)
c     call squarw(c(ja),ci,numi)
c     call apbc  (s(ia),h1,ci,1,numi,numi)
c     call ebc   (si,h1,c(ia),numi,1,numi)
c     call foldw (s(ja),si,numi)
      call squarw(s(ja),sj,numi)
      call ebc   (h1,sj,c(ia),numi,numi,1)
c         for tden-dry
      call squarw(c(ja),cj,numi)
      call apbc   (h1,cj,s(ia),numi,numi,1)
c
      call hi    (h1,0,numi,int,ladd)
c
      return
c
c     ----- is.eq.0, js.ne.0 -----
c
  802 continue
      numi=numsym(js+1)
      if (numi.lt.1) return
      numj=numsym(1)
      if (numj.lt.1) return
      jas=ja+wtw(mn(js+1),js+1)
c     call hi    (h1,js,numi,int,ladd)
c     call apbct (s(ia),h1,c(jas),1,numi,numj)
c     call apbtct(s(jas),c(ia),h1,numj,1,numi)
      call ebtc  (h1,s(jas),c(ia),numi,numj,1)
c         for tden-dry
      call apbtc  (h1,c(jas),s(ia),numi,numj,1)
      call hi    (h1,js,numi,int,ladd)
c
      return
c
c     ----- is.ne.0 -----
c
  810 continue
      if (js.ne.0) go to 811
c
c     ----- is.ne.0, js.eq.0 -----
c
      numi=numsym(is+1)
      if (numi.lt.1) return
      jas=ja+wab(mn(is+1))
c     call hi    (h1,is,numi,int,ladd)
c     call squarw(c(jas),ci,numi)
c     call apbc  (s(ia),h1,ci,1,numi,numi)
c     call ebc   (si,h1,c(ia),numi,1,numi)
c     call foldw (s(jas),si,numi)
      call squarw(s(jas),sj,numi)
      call ebc   (h1,sj,c(ia),numi,numi,1)
c         for tden-dry
      call squarw(c(jas),cj,numi)
      call apbc   (h1,cj,s(ia),numi,numi,1)
      call hi    (h1,is,numi,int,ladd)
c
      return
c
c     ----- is.ne.0, js.ne.0 -----
c
  811 continue
      isym=xor(is,js)
      numj=numsym(is+1)
      if (numj.lt.1) return
      numi=numsym(isym+1)
      if (numi.lt.1) return
      if (isym.lt.is) go to 812
      jas=ja+wtw(mn(isym+1),js+1)
c     call hi    (h1,isym,numi,int,ladd)
c     call apbct (s(ia),h1,c(jas),1,numi,numj)
c     call apbtct(s(jas),c(ia),h1,numj,1,numi)
      call ebtc  (h1,s(jas),c(ia),numi,numj,1)
c         for tden
      call apbtc  (h1,c(jas),s(ia),numi,numj,1)
      call hi    (h1,isym,numi,int,ladd)
c
      return
  812 continue
      jas=ja+wtw(mn(is+1),js+1)
c     call hi    (h1,isym,numi,int,ladd)
c     call apbc  (s(ia),h1,c(jas),1,numi,numj)
c     call apbc  (s(jas),h1,c(ia),numi,1,numj)
      call ebc   (h1,s(jas),c(ia),numi,numj,1)
c         for tden
      call apbc   (h1,c(jas),s(ia),numi,numj,1)
      call hi    (h1,isym,numi,int,ladd)
c
      return
c
c****************************** wy *************************************
c
      entry wy(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $           ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c..rlm      entry wy
c
      if (iseg.eq.15) go to 600
      if (iseg.eq.17) go to 630
      if (iseg.eq.21) go to 670
c
c     ----- wy iseg=2 -----
c
      iwy2=iwy2+1
      jja=ja
      nkl=nklww(is+1)
      do 430 i=mn(js+1),mx(js+1)
      if (i.lt.jming.or.i.gt.jmaxg) go to 430
      ij=ijww(arr+i)
c     call dot(t,int(ij+1),c(ia),nkl)
      s(jja)=s(jja)+val1*dot(int(ij+1),c(ia),nkl)
  430 jja=jja+1
c
      njj=numsym(js+1)
      iia=ia
c
      jinf=mn(js+1)
      if (jming.gt.jinf) jinf=jming
      jdiff=jinf-mn(js+1)
      jsup=mx(js+1)
      if (jmaxg.lt.jsup) jsup=jmaxg
      ist=ijww(arr+jinf)
      do 432 i=1,nkl
      ipt=ist+i
      jja=ja+jdiff
      t=0.0
      do 431 j=jinf,jsup
      t=t+int(ipt)*c(jja)
      ipt=ipt-nkl
  431 jja=jja+1
      s(iia)=s(iia)+t*val1
      iia=iia+1
  432 continue
      return
c
c     ----- wy iseg=15 -----
c
  600 continue
      iwy15=iwy15+1
      a=val1
      b=val1
      intal=1
      intb=3
      go to 632
c
c     ----- wy iseg=17 -----
c
  630 continue
      iwy17=iwy17+1
      a=val1
      b=0.0
      intal=3
      intb=3
  632 continue
c
c     ----- dm specific mod -----
c
c     a=a+a
c     b=b+b
c
c
c     ----- special loop for d.o.--s.o. interactions -----
c
      if (is.ne.0) go to 620
      numi=numsym(js+1)
      if (numi.lt.1) return
c     call hs(h1,js,numi,int,ladd)
      iskip=1
      ias=ia+wab(mn(js+1))
      jas=ja
      do 631 i=1,numi
c     s(ias)=s(ias)+h1(i)*c(jas)
c     s(jas)=s(jas)+h1(i)*c(ias)
      h1(i)=c(ias)*s(jas)+c(jas)*s(ias)
c
      jas=jas+1
      iskip=iskip+1
      ias=ias+iskip
  631 continue
c
      call hs(h1,js,numi,int,ladd)
c
      go to 620
c
c     ----- wy iseg=21 -----
c
  670 continue
      iwy21=iwy21+1
      a=val1
      b=val1*val2
      intal=tr1
      intb=tr2
c
c     ----- dm specific mod -----
c
c     a=a+a
c     b=b+b
c
  620 continue
      if (js.ne.0) go to 610
      if (is.ne.0) go to 602
c
c     ----- is.eq.0, js.eq.0 -----
c
      numi=numsym(1)
      if (numi.lt.1) return
c     call hi    (h1,0,numi,int,ladd)
c     call squarw(c(ia),ci,numi)
c     call apbc  (s(ja),h1,ci,1,numi,numi)
c     call ebc   (si,h1,c(ja),numi,1,numi)
c     call foldw (s(ia),si,numi)
      call squarw(c(ia),ci,numi)
      call ebc   (h1,ci,s(ja),numi,numi,1)
c         for tden-dry
      call squarw(s(ia),si,numi)
      call apbc   (h1,si,c(ja),numi,numi,1)
c
      call hi    (h1,0,numi,int,ladd)
c
      return
c
c     ----- is.ne.0, js.eq.0 -----
c
  602 continue
      numi=numsym(is+1)
      if (numi.lt.1) return
      numj=numsym(1)
      if (numj.lt.1) return
      ias=ia+wtw(mn(is+1),is+1)
c     call hi    (h1,is,numi,int,ladd)
c     call apbct (s(ja),h1,c(ias),1,numi,numj)
c     call apbtct(s(ias),c(ja),h1,numj,1,numi)
      call ebtc  (h1,c(ias),s(ja),numi,numj,1)
c         for tden
      call apbtc  (h1,s(ias),c(ja),numi,numj,1)
      call hi    (h1,is,numi,int,ladd)
c
      return
c
c     ----- js.ne.0 -----
c
  610 continue
      if (is.ne.0) go to 611
c
c     ----- is.eq.0, js.ne.0 -----
c
      numi=numsym(js+1)
      if (numi.lt.1) return
      ias=ia+wab(mn(js+1))
c     call hi    (h1,js,numi,int,ladd)
c     call squarw(c(ias),ci,numi)
c     call apbc  (s(ja),h1,ci,1,numi,numi)
c     call ebc   (si,h1,c(ja),numi,1,numi)
c     call foldw (s(ias),si,numi)
      call squarw(c(ias),ci,numi)
      call ebc   (h1,ci,s(ja),numi,numi,1)
c         for tden-dry
      call squarw(s(ias),si,numi)
      call apbc   (h1,si,c(ja),numi,numi,1)
c
      call hi    (h1,js,numi,int,ladd)
c
      return
c
c     ----- is.ne.0, js.ne.0 -----
c
  611 continue
      isym=xor(is,js)
      numj=numsym(js+1)
      if (numj.lt.1) return
      numi=numsym(isym+1)
      if (numi.lt.1) return
      if (isym.lt.js) go to 612
      ias=ia+wtw(mn(isym+1),is+1)
c     call hi    (h1,isym,numi,int,ladd)
c     call apbct (s(ja),h1,c(ias),1,numi,numj)
c     call apbtct(s(ias),c(ja),h1,numj,1,numi)
      call ebtc  (h1,c(ias),s(ja),numi,numj,1)
c         for tden
      call apbtc  (h1,s(ias),c(ja),numi,numj,1)
c
      call hi    (h1,isym,numi,int,ladd)
c
      return
  612 continue
      ias=ia+wtw(mn(js+1),is+1)
c     call hi    (h1,isym,numi,int,ladd)
c     call apbc  (s(ja),h1,c(ias),1,numi,numj)
c     call apbc  (s(ias),h1,c(ja),numi,1,numj)
      call ebc   (h1,c(ias),s(ja),numi,numj,1)
c        for tden-dry
      call apbc   (h1,s(ias),c(ja),numi,numj,1)
      call hi    (h1,isym,numi,int,ladd)
c
      return
c
c********************************** xx *********************************
c
c..rlm     entry xx
      entry xx(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $           ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c
      if (iseg.gt.8) go to 98
      go to (98,98,98,4,5,6,98,8),iseg
c
    4 continue
      ixx4=ixx4+1
      a=-val1*asqrt2
      b=-(a+a)
      go to 10
    5 continue
      ixx5=ixx5+1
      a=-val1*asqrt2+val3
      b=val1*sqrt2
   10 continue
      intal=1
      intb =2
      if (is.ne.0) go to 20
      ias=ia
      do 11 i=1,nsym
      isym=i-1
      numi=numsym(isym+1)
      if (numi.le.1) go to 11
c         no modification required for tden
      call square(c(ias),ci,numi)
      call square(s(ias),sj,numi)
      call ebct  (h1,sj,ci,numi,numi,numi)
      call hiis  (h1,isym,numi,int,kadd)
      ias=ias+numi*(numi-1)/2
   11 continue
      return
c
   20 continue
      ias=ia
      do 21 isym=1,nsym-1
      jsym=xor(isym,is)
      if (jsym.ge.isym) go to 21
      numi=numsym(isym+1)
      numj=numsym(jsym+1)
      numij=numi*numj
      if (numij.le.0) go to 21
c     call hiis  (h1,isym,numi,int,kadd)
c     call hiis  (h2,jsym,numj,int,kadd)
c     call apbct (s(ias),c(ias),h1,numj,numi,numi)
c     call apbc  (s(ias),h2,c(ias),numj,numj,numi)
c         no modification required for tden
      call ebtc  (h1,s(ias),c(ias),numi,numj,numi)
      call hiis  (h1,isym,numi,int,kadd)
      call ebct  (h1,s(ias),c(ias),numj,numi,numj)
      call hiis  (h1,jsym,numj,int,kadd)
c
      ias=ias+numij
   21 continue
      return
c
c     ----- iseg 6 and 8 cases -----
c
    6 continue
      ixx6=ixx6+1
      a=-val1*asqrt2+val3
      b=val1*sqrt2
      intal=3
      intau=1
      intb =2
      intad=1
      intbd=2
      go to 30
    8 continue
      ixx8=ixx8+1
      a=val1
      b=0.0
      intal=tr1
      intau=1
      intb=2
      intad=1
      intbd=2
   30 continue
c
c     ----- dm specific change -----
c
c     a=a+a
c     b=b+b
c
      if (is.ne.js) go to 40
      if (is.ne.0) go to 35
c
      ias=ia
      jas=ja
      do 31 i=1,nsym
      isym=i-1
      numi=numsym(isym+1)
      if (numi.le.1) go to 31
      call square(c(ias),ci,numi)
c     call square(c(jas),cj,numi)
c     call square(s(ias),si,numi)
      call square(s(jas),sj,numi)
      call ebct  (h1,sj,ci,numi,numi,numi)
c         for tden-dry
      call square(s(ias),si,numi)
      call square(c(jas),cj,numi)
      call apbct  (h1,cj,si,numi,numi,numi)
c     call ambc  (h1,ci,sj,numi,numi,numi)
      call hii   (h1,isym,numi,int,kadd)
      ioff=numi*(numi-1)/2
      ias=ias+ioff
      jas=jas+ioff
   31 continue
      return
   35 continue
      ias=ia
      jas=ja
      do 36 i=2,nsym
      isym=i-1
      jsym=xor(is,isym)
      if (jsym.ge.isym) go to 36
      numi=numsym(isym+1)
      numj=numsym(jsym+1)
      ioff=numi*numj
      if (ioff.lt.1) go to 36
c     call hii   (h1,isym,numi,int,kadd)
c     call hii   (h2,jsym,numj,int,kadd)
c     call apbct (s(jas),c(ias),h1,numj,numi,numi)
c     call apbc  (s(jas),h2,c(ias),numj,numj,numi)
c     call apbc  (s(ias),c(jas),h1,numj,numi,numi)
c     call apbtc (s(ias),h2,c(jas),numj,numj,numi)
      call ebct  (h2,s(jas),c(ias),numj,numi,numj)
      call ebtc  (h1,s(jas),c(ias),numi,numj,numi)
c         for tden-dry
      call apbct  (h2,c(jas),s(ias),numj,numi,numj)
      call apbtc  (h1,c(jas),s(ias),numi,numj,numi)
      call hii   (h1,isym,numi,int,kadd)
      call hii   (h2,jsym,numj,int,kadd)
      ias=ias+ioff
      jas=jas+ioff
   36 continue
      return
c
   40 continue
      if (is.ne.0) go to 50
c temp
c]]]]]]]]]]]]]]]]]]undebugged]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
c     write (itape6,'(' xx iseg 6 or 8, is=0 js#0')')
c end
      ias=ia
      do 42 i=1,nsym
      isym=i-1
      ksym=xor(js,isym)
      numi=numsym(isym+1)
      numk=numsym(ksym+1)
      if (numi.le.1.or.numk.le.0) go to 42
      if (ksym.gt.isym) go to 41
      ias=ia+wtx(mn(isym+1)+1,is+1)
      jas=ja+wtx(mn(isym+1),js+1)
      intal=1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call square(c(ias),ci,numi)
c     call apbtct(s(jas),h1,ci,numk,numi,numi)
c     call ebtct (si,c(jas),h1,numi,numk,numi)
c     call fold  (s(ias),si,numi)
      call square(s(ias),si,numi)
      call ebtct (h1,si,c(jas),numi,numi,numk)
c         for tden-dry
      call square(c(ias),ci,numi)
      call apbtct (h1,ci,s(jas),numi,numi,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 42
   41 continue
      ias=ia+wtx(mn(isym+1)+1,is+1)
      jas=ja+wtx(mn(ksym+1),js+1)
      intal=3
c     call hij   (h1,ksym,isym,numk,numi,int,kadd)
c     call square(c(ias),ci,numi)
c     call apbtct(s(jas),ci,h1,numi,numi,numk)
c     call ebtct (si,h1,c(jas),numi,numk,numi)
c     call fold  (s(ias),si,numi)
      call square(c(ias),ci,numi)
      call ebtct (h1,s(jas),ci,numk,numi,numi)
c         for tden-dry
      call square(s(ias),si,numi)
      call apbtct (h1,c(jas),si,numk,numi,numi)
      call hij   (h1,ksym,isym,numk,numi,int,kadd)
c
   42 continue
      return
c
   50 if (js.ne.0) go to 60
c
c]]]]]]]]]]]]]]]]]]]]]]undebugged]]]]]]]]]]]]]]]]]]]]]]]]
c     write (itape6,'(' xx iseg 6 or 8, is#0 js=0')')
c end
      do 52 i=1,nsym
      isym=i-1
      ksym=xor(is,isym)
      numi=numsym(isym+1)
      numk=numsym(ksym+1)
      if (numi.le.1.or.numk.le.0) go to 52
      if (ksym.gt.isym) go to 51
      ias=ia+wtx(mn(isym+1)  ,is+1)
      jas=ja+wtx(mn(isym+1)+1,js+1)
      intal=3
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call square(c(jas),cj,numi)
c     call apbtct(s(ias),h1,cj,numk,numi,numi)
c     call ebtct (sj,c(ias),h1,numi,numk,numi)
c     call fold  (s(jas),sj,numi)
      call square(s(jas),sj,numi)
      call ebtct (h1,sj,c(ias),numi,numi,numk)
c         for tden-dry
      call square(c(jas),cj,numi)
      call apbtct (h1,cj,s(ias),numi,numi,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 52
   51 continue
      ias=ia+wtx(mn(ksym+1)  ,is+1)
      jas=ja+wtx(mn(isym+1)+1,js+1)
      intal=1
c     call hij   (h1,ksym,isym,numk,numi,int,kadd)
c     call square(c(jas),cj,numi)
c     call apbtct(s(ias),cj,h1,numi,numi,numk)
c     call ebtct (sj,h1,c(ias),numi,numk,numi)
c     call fold  (s(jas),sj,numi)
      call square(c(jas),cj,numi)
      call ebtct (h1,s(ias),cj,numk,numi,numi)
c         for tden-dry
      call square(s(jas),sj,numi)
      call apbtct (h1,c(ias),sj,numk,numi,numi)
      call hij   (h1,ksym,isym,numk,numi,int,kadd)
c
   52 continue
      return
c
   60 continue
c temp
c]]]]]]]]]]]]]]]]]]]]]]]]undebugged]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
c     write (6,"(' xx iseg 6 or 8, is#0 js#0')")
c end
      do 66 i=1,nsym
      isym=i-1
      numi=numsym(isym+1)
      if (numi.le.0) go to 66
      jsym=xor(isym,is)
      if (jsym.ge.isym) go to 66
      numj=numsym(jsym+1)
      if (numj.le.0) go to 66
      ias=ia+wtx(mn(isym+1),is+1)
      ksym=xor(jsym,js)
      numk=numsym(ksym+1)
      if (numk.le.0) go to 63
      if (ksym.ge.jsym) go to 61
      jas=ja+wtx(mn(jsym+1),js+1)
      intal=1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call ambtct(s(jas),h1,c(ias),numk,numi,numj)
c     call ambtct(s(ias),c(jas),h1,numj,numk,numi)
      call embtct(h1,s(ias),c(jas),numi,numj,numk)
c         for tden-dry
      call ambtct(h1,c(ias),s(jas),numi,numj,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 63
   61 continue
      if (ksym.ge.isym) go to 62
      jas=ja+wtx(mn(ksym+1),js+1)
      intal=1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call apbc  (s(jas),c(ias),h1,numj,numi,numk)
c     call apbct (s(ias),c(jas),h1,numj,numk,numi)
      call ebtc  (h1,s(ias),c(jas),numi,numj,numk)
c         for tden
      call apbtc  (h1,c(ias),s(jas),numi,numj,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 63
   62 continue
      jas=ja+wtx(mn(ksym+1),js+1)
      intal=3
c     call hij   (h1,ksym,isym,numk,numi,int,kadd)
c     call apbct (s(jas),c(ias),h1,numj,numi,numk)
c     call apbc  (s(ias),c(jas),h1,numj,numk,numi)
      call ebtc  (h1,s(jas),c(ias),numk,numj,numi)
c         for tden-dry
      call apbtc  (h1,c(jas),s(ias),numk,numj,numi)
      call hij   (h1,ksym,isym,numk,numi,int,kadd)
c
c
c     c(i,j) and c(i,k) part
c
   63 continue
      ksym=xor(isym,js)
      numk=numsym(ksym+1)
      if (numk.le.0) go to 66
      if (ksym.ge.jsym) go to 64
      jas=ja+wtx(mn(isym+1),js+1)
      intal=1
c     call hij   (h1,jsym,ksym,numj,numk,int,kadd)
c     call apbtc (s(jas),h1,c(ias),numk,numj,numi)
c     call apbc  (s(ias),h1,c(jas),numj,numk,numi)
      call ebct  (h1,s(ias),c(jas),numj,numi,numk)
c         for tden-dry
      call apbct  (h1,c(ias),s(jas),numj,numi,numk)
      call hij   (h1,jsym,ksym,numj,numk,int,kadd)
c
      go to 66
   64 continue
      if (ksym.ge.isym) go to 65
      jas=ja+wtx(mn(isym+1),js+1)
      intal=3
c     call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c     call apbc  (s(jas),h1,c(ias),numk,numj,numi)
c     call apbtc (s(ias),h1,c(jas),numj,numk,numi)
      call ebct  (h1,s(jas),c(ias),numk,numi,numj)
c         for tden -dry
      call apbct  (h1,c(jas),s(ias),numk,numi,numj)
      call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c
      go to 66
   65 continue
      jas=ja+wtx(mn(ksym+1),js+1)
      intal=3
c     call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c     call ambtct(s(jas),c(ias),h1,numi,numj,numk)
c     call ambtct(s(ias),h1,c(jas),numj,numk,numi)
      call embtct(h1,s(jas),c(ias),numk,numi,numj)
c         for tden-dry
      call ambtct(h1,c(jas),s(ias),numk,numi,numj)
      call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c
   66 continue
      return
c
c
c
   98 continue
      write (itape6,99) iseg
   99 format (//,' ***** error in xx entry, iseg=',i5,//)
      stop
c
c********************************** ww *********************************
c
      entry ww(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $           ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c..rlm     entry ww
c
      if (iseg.gt.8) go to 198
      go to (198,198,198,104,105,106,198,108),iseg
c
  104 continue
      iww4=iww4+1
      a=-val1*asqrt2
      b=-(a+a)
      go to 110
  105 continue
      iww5=iww5+1
      a=-val1*asqrt2
      b=-(a+a)
  110 continue
      intal=1
      intb =2
      if (is.ne.0) go to 120
      ias=ia
      do 111 i=1,nsym
      isym=i-1
      numi=numsym(isym+1)
      if (numi.lt.1) go to 111
c     call hiis  (h1,isym,numi,int,kadd)
c     call squarw(c(ias),ci,numi)
c     call ebc   (sj,h1,ci,numi,numi,numi)
c     call foldw (s(ias),sj,numi)
      call squarw(c(ias),ci,numi)
c         next three lines changed for tden-dry
      call squarw(s(ias),si,numi)
c     call ebct  (h1,ci,ci,numi,numi,numi)
      call ebct  (h1,ci,si,numi,numi,numi)
      call hiis  (h1,isym,numi,int,kadd)
c
      ias=ias+numi*(numi+1)/2
  111 continue
      return
c
  120 continue
      ias=ia
      do 121 isym=1,nsym-1
      jsym=xor(isym,is)
      if (jsym.ge.isym) go to 121
      numi=numsym(isym+1)
      numj=numsym(jsym+1)
      numij=numi*numj
      if (numij.le.0) go to 121
c     call hiis  (h1,isym,numi,int,kadd)
c     call hiis  (h2,jsym,numj,int,kadd)
c     call apbct (s(ias),c(ias),h1,numj,numi,numi)
c     call apbc  (s(ias),h2,c(ias),numj,numj,numi)
      call ebtc  (h1,s(ias),c(ias),numi,numj,numi)
c        no change required for tden  (8/20/84)
cdry  call apbtc  (h1,c(ias),s(ias),numi,numj,numi)
      call hiis  (h1,isym,numi,int,kadd)
      call ebct  (h1,s(ias),c(ias),numj,numi,numj)
cdry  call apbct  (h1,c(ias),s(ias),numj,numi,numj)
      call hiis  (h1,jsym,numj,int,kadd)
c
      ias=ias+numij
  121 continue
      return
c
c     ----- iseg 6 and 8 cases -----
c
  106 continue
      iww6=iww6+1
      a=-val1*asqrt2
      b=-(a+a)
      intal=3
      intau=1
      intb =2
      intad=1
      intbd=2
      go to 130
  108 continue
      iww8=iww8+1
c
c     ----- ww entries iseg=8 analytically zero -----
c
      return
  130 continue
c
c
c     ----- dm specific mod -----
c
c     a=a+a
c     b=b+b
c
      if (is.ne.js) go to 140
      if (is.ne.0) go to 135
      ias=ia
      jas=ja
      do 131 i=1,nsym
      isym=i-1
      numi=numsym(isym+1)
      if (numi.lt.1) go to 131
c     call squarw(c(ias),ci,numi)
c     call squarw(c(jas),cj,numi)
c     call hii   (h1,isym,numi,int,kadd)
c     call ebc   (sj,h1,ci,numi,numi,numi)
c     call ebtc  (si,h1,cj,numi,numi,numi)
c     call foldw (s(ias),si,numi)
c     call foldw (s(jas),sj,numi)
      call squarw(s(jas),sj,numi)
      call squarw(c(ias),ci,numi)
      call ebct  (h1,sj,ci,numi,numi,numi)
c         for tden
      call squarw(c(jas),cj,numi)
      call squarw(s(ias),si,numi)
      call apbct  (h1,cj,si,numi,numi,numi)
      call hii   (h1,isym,numi,int,kadd)
c
      ioff=numi*(numi+1)/2
      ias=ias+ioff
      jas=jas+ioff
  131 continue
      return
  135 continue
      ias=ia
      jas=ja
      do 136 i=2,nsym
      isym=i-1
      jsym=xor(is,isym)
      if (jsym.ge.isym) go to 136
      numi=numsym(isym+1)
      numj=numsym(jsym+1)
      ioff=numi*numj
      if (ioff.lt.1) go to 136
c     call hii   (h1,isym,numi,int,kadd)
c     call hii   (h2,jsym,numj,int,kadd)
c     call apbct (s(jas),c(ias),h1,numj,numi,numi)
c     call apbc  (s(jas),h2,c(ias),numj,numj,numi)
c     call apbc  (s(ias),c(jas),h1,numj,numi,numi)
c     call apbtc (s(ias),h2,c(jas),numj,numj,numi)
      call ebtc  (h1,s(jas),c(ias),numi,numj,numi)
c         for tden-dry
      call apbtc  (h1,c(jas),s(ias),numi,numj,numi)
      call hii   (h1,isym,numi,int,kadd)
      call ebct  (h1,s(jas),c(ias),numj,numi,numj)
      call apbct  (h1,c(jas),s(ias),numj,numi,numj)
      call hii   (h1,jsym,numj,int,kadd)
c
      ias=ias+ioff
      jas=jas+ioff
  136 continue
      return
c
  140 continue
      if (is.ne.0) go to 150
      ias=ia
      do 142 i=1,nsym
      isym=i-1
      ksym=xor(js,isym)
      numi=numsym(isym+1)
      numk=numsym(ksym+1)
      if (numi.lt.1.or.numk.le.0) go to 142
      if (ksym.gt.isym) go to 141
      ias=ia+wab(mn(isym+1))
      jas=ja+wtw(mn(isym+1),js+1)
      intal=1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call squarw(c(ias),ci,numi)
c     call apbtct(s(jas),h1,ci,numk,numi,numi)
c     call ebtct (si,c(jas),h1,numi,numk,numi)
c     call foldw (s(ias),si,numi)
      call squarw(c(ias),ci,numi)
      call ebtct (h1,ci,s(jas),numi,numi,numk)
c         for tden-dry
      call squarw(s(ias),si,numi)
      call apbtct (h1,si,c(jas),numi,numi,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 142
  141 continue
      ias=ia+wab(mn(isym+1))
      jas=ja+wtw(mn(ksym+1),js+1)
      intal=3
c     call hij   (h1,ksym,isym,numk,numi,int,kadd)
c     call squarw(c(ias),ci,numi)
c     call apbtct(s(jas),ci,h1,numi,numi,numk)
c     call ebtct (si,h1,c(jas),numi,numk,numi)
c     call foldw (s(ias),si,numi)
      call squarw(c(ias),ci,numi)
      call ebtct (h1,s(jas),ci,numk,numi,numi)
c        for tden
      call squarw(s(ias),si,numi)
      call apbtct (h1,c(jas),si,numk,numi,numi)
      call hij    (h1,ksym,isym,numk,numi,int,kadd)
c
  142 continue
      return
c
  150 if (js.ne.0) go to 160
      do 152 i=1,nsym
      isym=i-1
      ksym=xor(is,isym)
      numi=numsym(isym+1)
      numk=numsym(ksym+1)
      if (numi.lt.1.or.numk.le.0) go to 152
      if (ksym.gt.isym) go to 151
      ias=ia+wtw(mn(isym+1)  ,is+1)
      jas=ja+wab(mn(isym+1))
      intal=3
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call squarw(c(jas),cj,numi)
c     call apbtct(s(ias),h1,cj,numk,numi,numi)
c     call ebtct (sj,c(ias),h1,numi,numk,numi)
c     call foldw (s(jas),sj,numi)
      call squarw(c(jas),cj,numi)
      call ebtct (h1,cj,s(ias),numi,numi,numk)
c         for tden
      call squarw(s(jas),sj,numi)
      call apbtct (h1,sj,c(ias),numi,numi,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 152
  151 continue
      ias=ia+wtw(mn(ksym+1)  ,is+1)
      jas=ja+wab(mn(isym+1))
      intal=1
c     call hij   (h1,ksym,isym,numk,numi,int,kadd)
c     call squarw(c(jas),cj,numi)
c     call apbtct(s(ias),cj,h1,numi,numi,numk)
c     call ebtct (sj,h1,c(ias),numi,numk,numi)
c     call foldw (s(jas),sj,numi)
      call squarw(s(jas),sj,numi)
      call ebtct (h1,c(ias),sj,numk,numi,numi)
c         for tden-dry
      call squarw(c(jas),cj,numi)
      call apbtct (h1,s(ias),cj,numk,numi,numi)
      call hij   (h1,ksym,isym,numk,numi,int,kadd)
c
  152 continue
      return
c
  160 continue
      do 166 i=1,nsym
      isym=i-1
      numi=numsym(isym+1)
      if (numi.le.0) go to 166
      jsym=xor(isym,is)
      if (jsym.ge.isym) go to 166
      numj=numsym(jsym+1)
      if (numj.le.0) go to 166
      ias=ia+wtw(mn(isym+1),is+1)
      ksym=xor(jsym,js)
      numk=numsym(ksym+1)
      if (numk.le.0) go to 163
      if (ksym.ge.jsym) go to 161
      jas=ja+wtw(mn(jsym+1),js+1)
      intal=1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call apbtct(s(jas),h1,c(ias),numk,numi,numj)
c     call apbtct(s(ias),c(jas),h1,numj,numk,numi)
      call ebtct (h1,c(ias),s(jas),numi,numj,numk)
c         for tden-dry
      call apbtct (h1,s(ias),c(jas),numi,numj,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 163
  161 continue
      if (ksym.ge.isym) go to 162
      jas=ja+wtw(mn(ksym+1),js+1)
      intal=1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call apbc  (s(jas),c(ias),h1,numj,numi,numk)
c     call apbct (s(ias),c(jas),h1,numj,numk,numi)
      call ebtc  (h1,c(ias),s(jas),numi,numj,numk)
c         for tden-dry
      call apbtc  (h1,s(ias),c(jas),numi,numj,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 163
  162 continue
      jas=ja+wtw(mn(ksym+1),js+1)
      intal=3
c     call hij   (h1,ksym,isym,numk,numi,int,kadd)
c     call apbct (s(jas),c(ias),h1,numj,numi,numk)
c     call apbc  (s(ias),c(jas),h1,numj,numk,numi)
      call ebtc  (h1,s(jas),c(ias),numk,numj,numi)
c         for tden-dry
      call apbtc  (h1,c(jas),s(ias),numk,numj,numi)
      call hij   (h1,ksym,isym,numk,numi,int,kadd)
c
c
c     ----- c(i,j) and c(i,k) part -----
c
  163 continue
      ksym=xor(isym,js)
      numk=numsym(ksym+1)
      if (numk.le.0) go to 166
      if (ksym.ge.jsym) go to 164
      jas=ja+wtw(mn(isym+1),js+1)
      intal=1
c     call hij   (h1,jsym,ksym,numj,numk,int,kadd)
c     call apbtc (s(jas),h1,c(ias),numk,numj,numi)
c     call apbc  (s(ias),h1,c(jas),numj,numk,numi)
      call ebct  (h1,c(ias),s(jas),numj,numi,numk)
c         for tden
      call apbct  (h1,s(ias),c(jas),numj,numi,numk)
      call hij   (h1,jsym,ksym,numj,numk,int,kadd)
c
      go to 166
  164 continue
      if (ksym.ge.isym) go to 165
      jas=ja+wtw(mn(isym+1),js+1)
      intal=3
c     call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c     call apbc  (s(jas),h1,c(ias),numk,numj,numi)
c     call apbtc (s(ias),h1,c(jas),numj,numk,numi)
      call ebct  (h1,s(jas),c(ias),numk,numi,numj)
c         for tden-dry
      call apbct  (h1,c(jas),s(ias),numk,numi,numj)
      call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c
      go to 166
  165 continue
      jas=ja+wtw(mn(ksym+1),js+1)
      intal=3
c     call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c     call apbtct(s(jas),c(ias),h1,numi,numj,numk)
c     call apbtct(s(ias),h1,c(jas),numj,numk,numi)
      call ebtct (h1,s(jas),c(ias),numk,numi,numj)
c         for tden-dry
      call apbtct (h1,c(jas),s(ias),numk,numi,numj)
      call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c
  166 continue
      return
c
c
c
  198 continue
      write (itape6,199) iseg
  199 format (//,' ***** error in ww entry, iseg=',i5,//)
      stop
c
c******************************* wx ************************************
c
      entry wx(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $           ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c..rlm     entry wx
c
      if (iseg.ne.7) go to 298
c
c     ----- iseg 7 case -----
c
      iwx7=iwx7+1
      a=val1*hsqrt3
c
c     ----- dm specific mod -----
c
c     a=a+a
c
      b=0.0d+00
      intal=tr1
      intau=1
      intb =2
      intad=1
      intbd=2
c
      if (is.ne.js) go to 240
      if (is.ne.0) go to 235
c
c     ----- isym=jsym=symmetric -----
c
      do 231 i=1,nsym
      isym=i-1
      numi=numsym(isym+1)
      if (numi.le.1) go to 231
      ias=ia+wab(mn(isym+1))
      jas=ja+wtx(mn(isym+1)+1,js+1)
c     call squarw(c(ias),ci,numi)
c     call square(c(jas),cj,numi)
c     call hii   (h1,isym,numi,int,kadd)
c     call ebc   (sj,h1,ci,numi,numi,numi)
c     call ebtc  (si,h1,cj,numi,numi,numi)
c     call foldw (s(ias),si,numi)
c     call fold  (s(jas),sj,numi)
      call square(s(jas),sj,numi)
      call squarw(c(ias),ci,numi)
      call ebct  (h1,sj,ci,numi,numi,numi)
c         for tden-dry
      call square(c(jas),cj,numi)
      call squarw(s(ias),si,numi)
      call apbct  (h1,cj,si,numi,numi,numi)
      call hii   (h1,isym,numi,int,kadd)
c
  231 continue
      return
c
c     ----- isym=jsym#symmetric -----
c
  235 continue
      do 236 i=2,nsym
      isym=i-1
      jsym=xor(is,isym)
      if (jsym.ge.isym) go to 236
      numi=numsym(isym+1)
      numj=numsym(jsym+1)
      if (numi.lt.1.or.numj.lt.1) go to 236
      ias=ia+wtw(mn(isym+1),is+1)
      jas=ja+wtx(mn(isym+1),js+1)
c     call hii   (h1,isym,numi,int,kadd)
c     call hii   (h2,jsym,numj,int,kadd)
c     call apbct (s(jas),c(ias),h1,numj,numi,numi)
c     call ambc  (s(jas),h2,c(ias),numj,numj,numi)
c     call apbc  (s(ias),c(jas),h1,numj,numi,numi)
c     call ambtc (s(ias),h2,c(jas),numj,numj,numi)
      call ebtc  (h1,s(jas),c(ias),numi,numj,numi)
c         for tden
      call apbtc  (h1,c(jas),s(ias),numi,numj,numi)
      call hii   (h1,isym,numi,int,kadd)
      call embct (h1,s(jas),c(ias),numj,numi,numj)
      call ambct (h1,c(jas),s(ias),numj,numi,numj)
      call hii   (h1,jsym,numj,int,kadd)
c
  236 continue
      return
c
  240 continue
      if (is.ne.0) go to 250
c
c     ----- isym=symmetric, jsym#symmetric -----
c
      ias=ia
      do 242 i=1,nsym
      isym=i-1
      ksym=xor(js,isym)
      numi=numsym(isym+1)
      numk=numsym(ksym+1)
      if (numi.lt.1.or.numk.le.0) go to 242
      if (ksym.gt.isym) go to 241
      ias=ia+wab(mn(isym+1))
      jas=ja+wtx(mn(isym+1),js+1)
      intal=1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call squarw(c(ias),ci,numi)
c     call ambtct(s(jas),h1,ci,numk,numi,numi)
c     call embtct(si,c(jas),h1,numi,numk,numi)
c     call foldw (s(ias),si,numi)
      call squarw(c(ias),ci,numi)
c         following line needed modification for tden
c     call embtct(h1,ci,c(jas),numi,numi,numk)
      call embtct(h1,ci,s(jas),numi,numi,numk)
c         for tden-dry
      call squarw(s(ias),si,numi)
      call ambtct(h1,si,c(jas),numi,numi,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 242
  241 continue
      ias=ia+wab(mn(isym+1))
      jas=ja+wtx(mn(ksym+1),js+1)
      intal=tr1
c     call hij   (h1,ksym,isym,numk,numi,int,kadd)
c     call squarw(c(ias),ci,numi)
c     call apbtct(s(jas),ci,h1,numi,numi,numk)
c     call ebtct (si,h1,c(jas),numi,numk,numi)
c     call foldw (s(ias),si,numi)
      call squarw(c(ias),ci,numi)
      call ebtct (h1,s(jas),ci,numk,numi,numi)
c         for tden-dry
      call squarw(s(ias),si,numi)
      call apbtct (h1,c(jas),si,numk,numi,numi)
      call hij   (h1,ksym,isym,numk,numi,int,kadd)
c
  242 continue
      return
c
c     ----- jsym=symmetric, isym#symmetric -----
c
  250 if (js.ne.0) go to 260
      do 252 i=1,nsym
      isym=i-1
      ksym=xor(is,isym)
      numi=numsym(isym+1)
      numk=numsym(ksym+1)
      if (numi.le.1.or.numk.le.0) go to 252
      if (ksym.gt.isym) go to 251
      ias=ia+wtw(mn(isym+1)  ,is+1)
      jas=ja+wtx(mn(isym+1)+1,js+1)
      intal=tr1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call square(c(jas),cj,numi)
c     call ambtct(s(ias),h1,cj,numk,numi,numi)
c     call embtct(sj,c(ias),h1,numi,numk,numi)
c     call fold  (s(jas),sj,numi)
      call square(s(jas),sj,numi)
      call embtct(h1,sj,c(ias),numi,numi,numk)
c         for tden
      call square(c(jas),cj,numi)
      call ambtct(h1,cj,s(ias),numi,numi,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 252
  251 continue
      ias=ia+wtw(mn(ksym+1)  ,is+1)
      jas=ja+wtx(mn(isym+1)+1,js+1)
      intal=1
c     call hij   (h1,ksym,isym,numk,numi,int,kadd)
c     call square(c(jas),cj,numi)
c     call apbtct(s(ias),cj,h1,numi,numi,numk)
c     call ebtct (sj,h1,c(ias),numi,numk,numi)
c     call fold  (s(jas),sj,numi)
      call square(s(jas),sj,numi)
      call ebtct (h1,c(ias),sj,numk,numi,numi)
c         for tden-dry
      call square(c(jas),cj,numi)
      call apbtct (h1,s(ias),cj,numk,numi,numi)
      call hij   (h1,ksym,isym,numk,numi,int,kadd)
c
  252 continue
      return
c
c     ----- isym#jsym, neither is symmetric -----
c
  260 continue
      do 266 i=1,nsym
      isym=i-1
      numi=numsym(isym+1)
      if (numi.le.0) go to 266
      jsym=xor(isym,is)
      if (jsym.ge.isym) go to 266
      numj=numsym(jsym+1)
      if (numj.le.0) go to 266
      ias=ia+wtw(mn(isym+1),is+1)
      ksym=xor(jsym,js)
      numk=numsym(ksym+1)
      if (numk.le.0) go to 263
c
c     ----- c(i,j) and c(j,k) part -----
c
      if (ksym.ge.jsym) go to 261
      jas=ja+wtx(mn(jsym+1),js+1)
      intal=1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call ambtct(s(jas),h1,c(ias),numk,numi,numj)
c     call ambtct(s(ias),c(jas),h1,numj,numk,numi)
      call embtct(h1,c(ias),s(jas),numi,numj,numk)
c         for tden
      call ambtct(h1,s(ias),c(jas),numi,numj,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 263
  261 continue
      if (ksym.ge.isym) go to 262
      jas=ja+wtx(mn(ksym+1),js+1)
      intal=1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call apbc  (s(jas),c(ias),h1,numj,numi,numk)
c     call apbct (s(ias),c(jas),h1,numj,numk,numi)
      call ebtc  (h1,c(ias),s(jas),numi,numj,numk)
c         for tden-dry
      call apbtc  (h1,s(ias),c(jas),numi,numj,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 263
  262 continue
      jas=ja+wtx(mn(ksym+1),js+1)
      intal=tr1
c     call hij   (h1,ksym,isym,numk,numi,int,kadd)
c     call apbct (s(jas),c(ias),h1,numj,numi,numk)
c     call apbc  (s(ias),c(jas),h1,numj,numk,numi)
      call ebtc  (h1,s(jas),c(ias),numk,numj,numi)
c         for tden-dry
      call apbtc  (h1,c(jas),s(ias),numk,numj,numi)
      call hij   (h1,ksym,isym,numk,numi,int,kadd)
c
c
c     ----- c(i,j) and c(i,k) part -----
c
  263 continue
      ksym=xor(isym,js)
      numk=numsym(ksym+1)
      if (numk.le.0) go to 266
      if (ksym.ge.jsym) go to 264
      jas=ja+wtx(mn(isym+1),js+1)
      intal=1
c     call hij   (h1,jsym,ksym,numj,numk,int,kadd)
c     call ambtc (s(jas),h1,c(ias),numk,numj,numi)
c     call ambc  (s(ias),h1,c(jas),numj,numk,numi)
      call embct (h1,c(ias),s(jas),numj,numi,numk)
c         for tden-dry
      call ambct (h1,s(ias),c(jas),numj,numi,numk)
      call hij   (h1,jsym,ksym,numj,numk,int,kadd)
c
      go to 266
  264 continue
      if (ksym.ge.isym) go to 265
      jas=ja+wtx(mn(isym+1),js+1)
      intal=tr1
c     call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c     call ambc  (s(jas),h1,c(ias),numk,numj,numi)
c     call ambtc (s(ias),h1,c(jas),numj,numk,numi)
      call embct (h1,s(jas),c(ias),numk,numi,numj)
c         for tden
      call ambct (h1,c(jas),s(ias),numk,numi,numj)
      call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c
      go to 266
  265 continue
      jas=ja+wtx(mn(ksym+1),js+1)
      intal=tr1
c     call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c     call apbtct(s(jas),c(ias),h1,numi,numj,numk)
c     call apbtct(s(ias),h1,c(jas),numj,numk,numi)
      call ebtct (h1,s(jas),c(ias),numk,numi,numj)
c         for tden-dry
      call apbtct (h1,c(jas),s(ias),numk,numi,numj)
      call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c
  266 continue
      return
c
c
c
  298 continue
      write (itape6,299) iseg
  299 format (//,' ***** error in wx entry, iseg=',i5,//)
      stop
c
c******************************* xw ************************************
c
      entry xw(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $           ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c
      if (iseg.ne.9) go to 398
c
c     ----- iseg 9 case -----
c
      ixw9=ixw9+1
      a=val1*hsqrt3
c
c     ----- dm specific mod -----
c
c     a=a+a
c
      b=0.0d+00
      intal=tr1
      intau=1
      intb =2
      intad=1
      intbd=2
c
      if (is.ne.js) go to 340
      if (is.ne.0) go to 335
c
c     ----- isym=jsym=symmetric -----
c
      do 331 i=1,nsym
      isym=i-1
      numi=numsym(isym+1)
      if (numi.le.1) go to 331
      ias=ia+wtx(mn(isym+1)+1,is+1)
      jas=ja+wab(mn(isym+1))
c     call square(c(ias),ci,numi)
c     call squarw(c(jas),cj,numi)
c     call hii   (h1,isym,numi,int,kadd)
c     call ebc   (sj,h1,ci,numi,numi,numi)
c     call ebtc  (si,h1,cj,numi,numi,numi)
c     call fold  (s(ias),si,numi)
c     call foldw (s(jas),sj,numi)
      call squarw(s(jas),sj,numi)
      call square(c(ias),ci,numi)
      call ebct  (h1,sj,ci,numi,numi,numi)
c         for tden-dry
      call squarw(c(jas),cj,numi)
      call square(s(ias),si,numi)
      call apbct  (h1,cj,si,numi,numi,numi)
      call hii   (h1,isym,numi,int,kadd)
c
  331 continue
      return
c
c     ----- isym=jsym#symmetric -----
c
  335 continue
      do 336 i=2,nsym
      isym=i-1
      jsym=xor(is,isym)
      if (jsym.ge.isym) go to 336
      numi=numsym(isym+1)
      numj=numsym(jsym+1)
      if (numi.lt.1.or.numj.lt.1) go to 336
      ias=ia+wtx(mn(isym+1),is+1)
      jas=ja+wtw(mn(isym+1),js+1)
c     call hii   (h1,isym,numi,int,kadd)
c     call hii   (h2,jsym,numj,int,kadd)
c     call apbct (s(jas),c(ias),h1,numj,numi,numi)
c     call ambc  (s(jas),h2,c(ias),numj,numj,numi)
c     call apbc  (s(ias),c(jas),h1,numj,numi,numi)
c     call ambtc (s(ias),h2,c(jas),numj,numj,numi)
      call ebtc  (h1,s(jas),c(ias),numi,numj,numi)
c         for tden
      call apbtc  (h1,c(jas),s(ias),numi,numj,numi)
      call hii   (h1,isym,numi,int,kadd)
      call embct (h1,s(jas),c(ias),numj,numi,numj)
      call ambct (h1,c(jas),s(ias),numj,numi,numj)
      call hii   (h1,jsym,numj,int,kadd)
c
  336 continue
      return
c
  340 continue
      if (is.ne.0) go to 350
c
c     ----- isym=symmetric, jsym#symmetric -----
c
      ias=ia
      do 342 i=1,nsym
      isym=i-1
      ksym=xor(js,isym)
      numi=numsym(isym+1)
      numk=numsym(ksym+1)
      if (numi.le.1.or.numk.le.0) go to 342
      if (ksym.gt.isym) go to 341
      ias=ia+wtx(mn(isym+1)+1,is+1)
      jas=ja+wtw(mn(isym+1)  ,js+1)
      intal=1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call square(c(ias),ci,numi)
c     call ambtct(s(jas),h1,ci,numk,numi,numi)
c     call embtct(si,c(jas),h1,numi,numk,numi)
c     call fold  (s(ias),si,numi)
      call square(c(ias),ci,numi)
      call embtct(h1,ci,s(jas),numi,numi,numk)
c         tden-dry
      call square(s(ias),si,numi)
      call ambtct(h1,si,c(jas),numi,numi,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 342
  341 continue
      ias=ia+wtx(mn(isym+1)+1,is+1)
      jas=ja+wtw(mn(ksym+1)  ,js+1)
      intal=tr1
c     call hij   (h1,ksym,isym,numk,numi,int,kadd)
c     call square(c(ias),ci,numi)
c     call apbtct(s(jas),ci,h1,numi,numi,numk)
c     call ebtct (si,h1,c(jas),numi,numk,numi)
c     call fold  (s(ias),si,numi)
      call square(c(ias),ci,numi)
      call ebtct (h1,s(jas),ci,numk,numi,numi)
c         for tden-dry
      call square(s(ias),si,numi)
      call apbtct (h1,c(jas),si,numk,numi,numi)
      call hij   (h1,ksym,isym,numk,numi,int,kadd)
c
  342 continue
      return
c
c     ----- jsym=symmetric, isym#symmetric -----
c
  350 if (js.ne.0) go to 360
      do 352 i=1,nsym
      isym=i-1
      ksym=xor(is,isym)
      numi=numsym(isym+1)
      numk=numsym(ksym+1)
      if (numi.lt.1.or.numk.le.0) go to 352
      if (ksym.gt.isym) go to 351
      ias=ia+wtx(mn(isym+1)  ,is+1)
      jas=ja+wab(mn(isym+1))
      intal=tr1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call squarw(c(jas),cj,numi)
c     call ambtct(s(ias),h1,cj,numk,numi,numi)
c     call embtct(sj,c(ias),h1,numi,numk,numi)
c     call foldw (s(jas),sj,numi)
      call squarw(s(jas),sj,numi)
      call embtct(h1,sj,c(ias),numi,numi,numk)
c         for tden-dry
      call squarw(c(jas),cj,numi)
      call ambtct(h1,cj,s(ias),numi,numi,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 352
  351 continue
      ias=ia+wtx(mn(ksym+1)  ,is+1)
      jas=ja+wab(mn(isym+1))
      intal=1
c     call hij   (h1,ksym,isym,numk,numi,int,kadd)
c     call squarw(c(jas),cj,numi)
c     call apbtct(s(ias),cj,h1,numi,numi,numk)
c     call ebtct (sj,h1,c(ias),numi,numk,numi)
c     call foldw (s(jas),sj,numi)
      call squarw(s(jas),sj,numi)
      call ebtct (h1,c(ias),sj,numk,numi,numi)
c         tden-dry
      call squarw(c(jas),cj,numi)
      call apbtct (h1,s(ias),cj,numk,numi,numi)
      call hij   (h1,ksym,isym,numk,numi,int,kadd)
c
  352 continue
      return
c
c     ----- isym#jsym, neither is symmetric -----
c
  360 continue
      do 366 i=1,nsym
      isym=i-1
      numi=numsym(isym+1)
      if (numi.le.0) go to 366
      jsym=xor(isym,is)
      if (jsym.ge.isym) go to 366
      numj=numsym(jsym+1)
      if (numj.le.0) go to 366
      ias=ia+wtx(mn(isym+1),is+1)
      ksym=xor(jsym,js)
      numk=numsym(ksym+1)
      if (numk.le.0) go to 363
c
c     ----- c(i,j) and c(j,k) part -----
c
      if (ksym.ge.jsym) go to 361
      jas=ja+wtw(mn(jsym+1),js+1)
      intal=1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call apbtct(s(jas),h1,c(ias),numk,numi,numj)
c     call apbtct(s(ias),c(jas),h1,numj,numk,numi)
      call ebtct (h1,c(ias),s(jas),numi,numj,numk)
c         for tden
      call apbtct (h1,s(ias),c(jas),numi,numj,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 363
  361 continue
      if (ksym.ge.isym) go to 362
      jas=ja+wtw(mn(ksym+1),js+1)
      intal=1
c     call hij   (h1,isym,ksym,numi,numk,int,kadd)
c     call apbc  (s(jas),c(ias),h1,numj,numi,numk)
c     call apbct (s(ias),c(jas),h1,numj,numk,numi)
      call ebtc  (h1,c(ias),s(jas),numi,numj,numk)
c         for tden-dry
      call apbtc  (h1,s(ias),c(jas),numi,numj,numk)
      call hij   (h1,isym,ksym,numi,numk,int,kadd)
c
      go to 363
  362 continue
      jas=ja+wtw(mn(ksym+1),js+1)
      intal=tr1
c     call hij   (h1,ksym,isym,numk,numi,int,kadd)
c     call apbct (s(jas),c(ias),h1,numj,numi,numk)
c     call apbc  (s(ias),c(jas),h1,numj,numk,numi)
      call ebtc  (h1,s(jas),c(ias),numk,numj,numi)
c         for tden-dry
      call apbtc  (h1,c(jas),s(ias),numk,numj,numi)
      call hij   (h1,ksym,isym,numk,numi,int,kadd)
c
c
c     ----- c(i,j) and c(i,k) part -----
c
  363 continue
      ksym=xor(isym,js)
      numk=numsym(ksym+1)
      if (numk.le.0) go to 366
      if (ksym.ge.jsym) go to 364
      jas=ja+wtw(mn(isym+1),js+1)
      intal=1
c     call hij   (h1,jsym,ksym,numj,numk,int,kadd)
c     call ambtc (s(jas),h1,c(ias),numk,numj,numi)
c     call ambc  (s(ias),h1,c(jas),numj,numk,numi)
      call embct (h1,c(ias),s(jas),numj,numi,numk)
c         for tden-dry
      call ambct (h1,s(ias),c(jas),numj,numi,numk)
      call hij   (h1,jsym,ksym,numj,numk,int,kadd)
c
      go to 366
  364 continue
      if (ksym.ge.isym) go to 365
      jas=ja+wtw(mn(isym+1),js+1)
      intal=tr1
c     call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c     call ambc  (s(jas),h1,c(ias),numk,numj,numi)
c     call ambtc (s(ias),h1,c(jas),numj,numk,numi)
      call embct (h1,s(jas),c(ias),numk,numi,numj)
c         for tden-dry
      call ambct (h1,c(jas),s(ias),numk,numi,numj)
      call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c
      go to 366
  365 continue
      jas=ja+wtw(mn(ksym+1),js+1)
      intal=tr1
c     call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c     call ambtct(s(jas),c(ias),h1,numi,numj,numk)
c     call ambtct(s(ias),h1,c(jas),numj,numk,numi)
      call embtct(h1,s(jas),c(ias),numk,numi,numj)
c         for tden-dry
      call ambtct(h1,c(jas),s(ias),numk,numi,numj)
      call hij   (h1,ksym,jsym,numk,numj,int,kadd)
c
  366 continue
      return
c
c
c
  398 continue
      write (itape6,399) iseg
  399 format (//,' ***** error in xw entry, iseg=',i5,//)
      stop
c
      end
