*deck @(#)initst.f	1.2  7/30/91
      subroutine initst(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
c
c
c
      implicit real*8 (a-h,o-z)
c
      integer xor
      integer arr,tr1,tr2,asm,os,wtw,wtx,wty,wab,ss,symorb
      integer bmax,orbfrm
      real*8 int(nmax),c(nwks,mxvc),s(nwks,mxvc)
c
      dimension ci(1),cj(1),si(1),sj(1),h1(1),h2(1)
      dimension kadd(symorb),ladd(symorb),ijadd(numij),wtw(orbfrm,nsym)
      dimension wtx(orbfrm,nsym),wty(orbfrm),wab(orbfrm),ss(norbs)
      dimension ijxx(numij ),nklxx(nsym,orbfrm),ijww(numij )
      dimension nklww(nsym,orbfrm),klxx(1),klww(1)
c
      common /coefs/  a,b,intal,intau,intb,intad,intbd
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
c  universal identity of the objects in these common
c     common /all/ val1,val2,val3,arr,itr1,itr2,ia,ja,itype,isegt
c    *,lvfrm1,nlwki,nlwkj,imax,imin
c     common /symq/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
c    #,             numsym(8)
      common /symq/ asm,js,is,mx(8),mn(8),os(8),numsym(8)
      common /io/     itape5,itape6
      common /all/val1,val2,val3,arr,tr1,tr2,ia,ja,m,iseg,n,n1,n2
     *,           imax,imin
      common /count/  icount,ixx4,ixx5,ixx6,ixx8,iww4,iww5,iww6,iww8
     #,               iwx7,ixw9,ixy3,ixy16,ixy18,ixy22,iwy2,iwy15,iwy17
     #,               iwy21,iyx20,iyw19
      common /x4x901/ nijvir,nrefs
      common /minmax/ iming,imaxg,jming,jmaxg
      common /nvectr/ nvc,mxvc
      real*8 val1,val2,val3,t1(64)
c
c
      sqrt2=sqrt(2.0d+00)
      sqrt3=sqrt(3.0d+00)
      sqt1p5=sqrt(1.5d+00)
      asqrt2=1.0d+00/sqrt2
      hsqrt3=0.5d+00*sqrt(3.0d+00)
      return
c
c***************************** xx 4x ***********************************
c
      entry xx4x
c
      jminh=jming
      if (jming.eq.1) jminh=0
c
      do 402 i=imin,imax
         ip=i*(i-1)/2
         jsy=xor(is,ss(i))
         jmin=mn(jsy+1)
         if (jmin.ge.i) go to 401
         nkl=nklxx(is+1,n)
         jmax=mx(jsy+1)
         if (jmax.ge.i) jmax=i-1
         iia=ia+wtx(i,is+1)
         do 400 j=jmin,jmax
            if (j.lt.jming.or.j.gt.jmaxg) go to 400
            ijxxpt=ijxx(ip+j)+1
            do 1001 ivc=1,nvc
               ciia=c(iia,ivc)
               if (ciia.ne.0.0d+00) then
                  call saxpy(nkl,ciia,int(ijxxpt),1,s(ia,ivc),1)
               end if
 1001       continue
  400    iia=iia+1
  401    continue
  402 continue
      return
c
c***************************** ww 4x ***********************************
c
      entry ww4x
c
      jminh=jming
      if (jming.eq.1) jminh=0
c
      do 414 i=imin,imax
         ip=i*(i-1)/2
         jsy=xor(is,ss(i))
         jmin=mn(jsy+1)
         if (jmin.gt.i) go to 413
         nkl=nklww(is+1,n)
         jmax=mx(jsy+1)
         if (jmax.gt.i) jmax=i
         if (i.ne.jmin) go to 410
         iia=ia+wab(i)
         go to 411
  410    continue
         iia=ia+wtw(i,is+1)
  411    continue
         do 412 j=jmin,jmax
            if (j.lt.jming.or.j.gt.jmaxg) go to 412
            ijwwpt=ijww(ip+j)+1
            do 1002 ivc=1,nvc
               ciia=c(iia,ivc)
               if (ciia.ne.0.0d+00) then
                  call saxpy(nkl,ciia,int(ijwwpt),1,s(ia,ivc),1)
               end if
 1002       continue
  412    iia=iia+1
  413    continue
  414 continue
      return
c
c******************************** yx ***********************************
c
      entry yx
c
      iyx20=iyx20+1
      a=val1
      b=val1*val2
      intal=tr1
      intb=tr2
      if (is.ne.0) go to 710
      if (js.ne.0) go to 702
c
c     ----- is.eq.0, js.eq.0 -----
c
      numi=numsym(1)
      if (numi.le.1) return
      call hi    (h1(1),0,numi)
      do 1003 ivc=1,nvc
         call square(c(ja,ivc),ci(1),numi)
         call apbc  (s(ia,ivc),h1(1),ci(1),1,numi,numi)
         call ebc   (si(1),h1(1),c(ia,ivc),numi,1,numi)
         call fold  (s(ja,ivc),si(1),numi)
 1003 continue
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
      call hi    (h1(1),js,numi)
      do 1004 ivc=1,nvc
         call apbct (s(ia,ivc),h1(1),c(jas,ivc),1,numi,numj)
         call atpbc(s(jas,ivc),h1(1),c(ia,ivc),numi,1,numj)
 1004 continue
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
      call hi    (h1(1),is,numi)
      do 1005 ivc=1,nvc
         call square(c(jas,ivc),ci(1),numi)
         call apbc  (s(ia,ivc),h1(1),ci(1),1,numi,numi)
         call ebc   (si(1),h1(1),c(ia,ivc),numi,1,numi)
         call fold  (s(jas,ivc),si(1),numi)
 1005 continue
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
      call hi    (h1(1),isym,numi)
      do 1006 ivc=1,nvc
         call apbct (s(ia,ivc),h1(1),c(jas,ivc),1,numi,numj)
         call atpbc(s(jas,ivc),h1(1),c(ia,ivc),numi,1,numj)
 1006 continue
      return
c
  712 continue
      jas=ja+wtx(mn(is+1),js+1)
      call hi    (h1(1),isym,numi)
      do 1007 ivc=1,nvc
         call ambc  (s(ia,ivc),h1(1),c(jas,ivc),1,numi,numj)
         call ambc  (s(jas,ivc),h1(1),c(ia,ivc),numi,1,numj)
 1007 continue
      return
c
c****************************** xy *************************************
c
      entry xy
      if (iseg.eq.16) go to 500
      if (iseg.eq.18) go to 530
      if (iseg.eq.22) go to 570
c
c     ----- xy iseg=3 -----
c
      ixy3=ixy3+1
      jja=ja
      nkl=nklxx(is+1,n)
      do 420 i=mn(js+1),mn(js+1)+numsym(js+1)-1
         if (i.lt.jming.or.i.gt.jmaxg) go to 420
         ij=ijxx(arr+i)
         do 1008 ivc=1,nvc
            s(jja,ivc)=s(jja,ivc)+
     $           val1*sdot(nkl,int(ij+1),1,c(ia,ivc),1)
 1008    continue
  420 jja=jja+1
c
      nkloff=nklxx(is+1,orbfrm)
      iia=ia
      jinf=mn(js+1)
      if (jming.gt.jinf) jinf=jming
      jdiff=jinf-mn(js+1)
      jsup=mn(js+1)+numsym(js+1)-1
      if (jmaxg.lt.jsup) jsup=jmaxg
      ist=ijxx(arr+jinf)
      do 422 i=1,nkl
         jja=ja+jdiff
         ipt=ist+i
         call rzero(t1,nvc)
         do 421 j=jinf,jsup
            do 1009 ivc=1,nvc
               t1(ivc)=t1(ivc)+int(ipt)*c(jja,ivc)
 1009       continue
            ipt=ipt-nkloff
  421    jja=jja+1
         do 1010 ivc=1,nvc
            s(iia,ivc)=s(iia,ivc)+t1(ivc)*val1
 1010    continue
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
      if (js.ne.0) go to 510
      if (is.ne.0) go to 502
c
c     ----- is.eq.0, js.eq.0 -----
c
      numi=numsym(1)
      if (numi.le.1) return
      call hi    (h1(1),0,numi)
      do 1011 ivc=1,nvc
         call square(c(ia,ivc),ci(1),numi)
         call apbc  (s(ja,ivc),h1(1),ci(1),1,numi,numi)
         call ebc   (si(1),h1(1),c(ja,ivc),numi,1,numi)
         call fold  (s(ia,ivc),si(1),numi)
 1011 continue
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
      call hi    (h1(1),is,numi)
      do 1012 ivc=1,nvc
         call apbct (s(ja,ivc),h1(1),c(ias,ivc),1,numi,numj)
         call atpbc(s(ias,ivc),h1(1),c(ja,ivc),numi,1,numj)
 1012 continue
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
      call hi    (h1(1),js,numi)
      do 1013 ivc=1,nvc
         call square(c(ias,ivc),ci(1),numi)
         call apbc  (s(ja,ivc),h1(1),ci(1),1,numi,numi)
         call ebc   (si(1),h1(1),c(ja,ivc),numi,1,numi)
         call fold  (s(ias,ivc),si(1),numi)
 1013 continue
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
      call hi    (h1(1),isym,numi)
      do 1014 ivc=1,nvc
         call apbct (s(ja,ivc),h1(1),c(ias,ivc),1,numi,numj)
         call atpbc(s(ias,ivc),h1(1),c(ja,ivc),numi,1,numj)
 1014 continue
      return
c
  512 continue
      ias=ia+wtx(mn(js+1),is+1)
      call hi    (h1(1),isym,numi)
      do 1015 ivc=1,nvc
         call ambc  (s(ja,ivc),h1(1),c(ias,ivc),1,numi,numj)
         call ambc  (s(ias,ivc),h1(1),c(ja,ivc),numi,1,numj)
 1015 continue
      return
c
c******************************** yw ***********************************
c
      entry yw
c
      iyw19=iyw19+1
      a=val1
      b=val1*val2
      intal=tr1
      intb=tr2
      if (is.ne.0) go to 810
      if (js.ne.0) go to 802
c
c     ----- is.eq.0, js.eq.0 -----
c
      numi=numsym(1)
      if (numi.lt.1) return
      call hi    (h1(1),0,numi)
      do 1016 ivc=1,nvc
         call squarw(c(ja,ivc),ci(1),numi)
         call apbc  (s(ia,ivc),h1(1),ci(1),1,numi,numi)
         call ebc   (si(1),h1(1),c(ia,ivc),numi,1,numi)
         call foldw (s(ja,ivc),si(1),numi)
 1016 continue
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
      call hi    (h1(1),js,numi)
      do 1017 ivc=1,nvc
         call apbct (s(ia,ivc),h1(1),c(jas,ivc),1,numi,numj)
         call atpbc(s(jas,ivc),h1(1),c(ia,ivc),numi,1,numj)
 1017 continue
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
      call hi    (h1(1),is,numi)
      do 1018 ivc=1,nvc
         call squarw(c(jas,ivc),ci(1),numi)
         call apbc  (s(ia,ivc),h1(1),ci(1),1,numi,numi)
         call ebc   (si(1),h1(1),c(ia,ivc),numi,1,numi)
         call foldw (s(jas,ivc),si(1),numi)
 1018 continue
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
      call hi    (h1(1),isym,numi)
      do 1019 ivc=1,nvc
         call apbct (s(ia,ivc),h1(1),c(jas,ivc),1,numi,numj)
         call atpbc(s(jas,ivc),h1(1),c(ia,ivc),numi,1,numj)
 1019 continue
      return
c
  812 continue
      jas=ja+wtw(mn(is+1),js+1)
      call hi    (h1(1),isym,numi)
      do 1020 ivc=1,nvc
         call apbc  (s(ia,ivc),h1(1),c(jas,ivc),1,numi,numj)
         call apbc  (s(jas,ivc),h1(1),c(ia,ivc),numi,1,numj)
 1020 continue
      return
c
c****************************** wy *************************************
c
      entry wy
c
      if (iseg.eq.15) go to 600
      if (iseg.eq.17) go to 630
      if (iseg.eq.21) go to 670
c
c     ----- wy iseg=2 -----
c
      iwy2=iwy2+1
      jja=ja
      nkl=nklww(is+1,n)
      do 430 i=mn(js+1),mn(js+1)+numsym(js+1)-1
         if (i.lt.jming.or.i.gt.jmaxg) go to 430
         ij=ijww(arr+i)
         do 1021 ivc=1,nvc
            s(jja,ivc)=s(jja,ivc)+
     $           val1*sdot(nkl,int(ij+1),1,c(ia,ivc),1)
 1021    continue
  430 jja=jja+1
c
      nkloff=nklww(is+1,orbfrm)
      njj=numsym(js+1)
      iia=ia
c
      jinf=mn(js+1)
      if (jming.gt.jinf) jinf=jming
      jdiff=jinf-mn(js+1)
      jsup=mn(js+1)+numsym(js+1)-1
      if (jmaxg.lt.jsup) jsup=jmaxg
      ist=ijww(arr+jinf)
      do 432 i=1,nkl
         ipt=ist+i
         jja=ja+jdiff
         call rzero(t1,nvc)
         do 431 j=jinf,jsup
            do 1022 ivc=1,nvc
               t1(ivc)=t1(ivc)+int(ipt)*c(jja,ivc)
 1022       continue
            ipt=ipt-nkloff
  431    jja=jja+1
         do 1023 ivc=1,nvc
            s(iia,ivc)=s(iia,ivc)+t1(ivc)*val1
 1023    continue
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
c     ----- special loop for d.o.--s.o. interactions -----
c
      if (is.ne.0) go to 620
      numi=numsym(js+1)
      if (numi.lt.1) return
      call hs(h1(1),js,numi)
      iskip=1
      ias=ia+wab(mn(js+1))
      jas=ja
      do 631 i=1,numi
         do 1024 ivc=1,nvc
            s(ias,ivc)=s(ias,ivc)+h1(i)*c(jas,ivc)
            s(jas,ivc)=s(jas,ivc)+h1(i)*c(ias,ivc)
 1024    continue
         jas=jas+1
         iskip=iskip+1
         ias=ias+iskip
  631 continue
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
  620 continue
      if (js.ne.0) go to 610
      if (is.ne.0) go to 602
c
c     ----- is.eq.0, js.eq.0 -----
c
      numi=numsym(1)
      if (numi.lt.1) return
      call hi    (h1(1),0,numi)
      do 1025 ivc=1,nvc
         call squarw(c(ia,ivc),ci(1),numi)
         call apbc  (s(ja,ivc),h1(1),ci(1),1,numi,numi)
         call ebc   (si(1),h1(1),c(ja,ivc),numi,1,numi)
         call foldw (s(ia,ivc),si(1),numi)
 1025 continue
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
      call hi    (h1(1),is,numi)
      do 1026 ivc=1,nvc
         call apbct (s(ja,ivc),h1(1),c(ias,ivc),1,numi,numj)
         call atpbc(s(ias,ivc),h1(1),c(ja,ivc),numi,1,numj)
 1026 continue
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
      call hi    (h1(1),js,numi)
      do 1027 ivc=1,nvc
         call squarw(c(ias,ivc),ci(1),numi)
         call apbc  (s(ja,ivc),h1(1),ci(1),1,numi,numi)
         call ebc   (si(1),h1(1),c(ja,ivc),numi,1,numi)
         call foldw (s(ias,ivc),si(1),numi)
 1027 continue
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
      call hi    (h1(1),isym,numi)
      do 1028 ivc=1,nvc
         call apbct (s(ja,ivc),h1(1),c(ias,ivc),1,numi,numj)
         call atpbc(s(ias,ivc),h1(1),c(ja,ivc),numi,1,numj)
 1028 continue
      return
c
  612 continue
      ias=ia+wtw(mn(js+1),is+1)
      call hi    (h1(1),isym,numi)
      do 1029 ivc=1,nvc
         call apbc  (s(ja,ivc),h1(1),c(ias,ivc),1,numi,numj)
         call apbc  (s(ias,ivc),h1(1),c(ja,ivc),numi,1,numj)
 1029 continue
      return
c
c********************************** xx *********************************
c
      entry xx
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
         call hiis  (h1(1),isym,numi)
         do 1030 ivc=1,nvc
            call square(c(ias,ivc),ci(1),numi)
            call ebc   (sj(1),h1(1),ci(1),numi,numi,numi)
            call fold  (s(ias,ivc),sj(1),numi)
 1030    continue
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
         call hiis  (h1(1),isym,numi)
         call hiis  (h2(1),jsym,numj)
         do 1031 ivc=1,nvc
            call atpbct(s(ias,ivc),h1(1),c(ias,ivc),numi,numi,numj)
            call apbc  (s(ias,ivc),h2(1),c(ias,ivc),numj,numj,numi)
 1031    continue
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
      if (is.ne.js) go to 40
      if (is.ne.0) go to 35
      ias=ia
      jas=ja
      do 31 i=1,nsym
         isym=i-1
         numi=numsym(isym+1)
         if (numi.le.1) go to 31
         call hii   (h1(1),isym,numi)
         do 1032 ivc=1,nvc
            call square(c(ias,ivc),ci(1),numi)
            call square(c(jas,ivc),cj(1),numi)
            call ebc   (sj(1),h1(1),ci(1),numi,numi,numi)
            call ebtc  (si(1),h1(1),cj(1),numi,numi,numi)
            call fold  (s(ias,ivc),si(1),numi)
            call fold  (s(jas,ivc),sj(1),numi)
 1032    continue
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
         call hii   (h1(1),isym,numi)
         call hii   (h2(1),jsym,numj)
         do 1033 ivc=1,nvc
            call atpbct(s(jas,ivc),h1(1),c(ias,ivc),numi,numi,numj)
            call apbc  (s(jas,ivc),h2(1),c(ias,ivc),numj,numj,numi)
            call apbctt(s(ias,ivc),h1(1),c(jas,ivc),numi,numi,numj)
            call apbtc (s(ias,ivc),h2(1),c(jas,ivc),numj,numj,numi)
 1033    continue
         ias=ias+ioff
         jas=jas+ioff
   36 continue
      return
c
   40 continue
      if (is.ne.0) go to 50
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
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1035 ivc=1,nvc
            call square(c(ias,ivc),ci(1),numi)
            call apbtct(s(jas,ivc),h1(1),ci(1),numk,numi,numi)
            call atebc(si(1),h1(1),c(jas,ivc),numi,numk,numi)
            call fold  (s(ias,ivc),si(1),numi)
 1035    continue
         go to 42
   41    continue
         ias=ia+wtx(mn(isym+1)+1,is+1)
         jas=ja+wtx(mn(ksym+1),js+1)
         intal=3
         call hij   (h1(1),ksym,isym,numk,numi)
         do 1036 ivc=1,nvc
            call square(c(ias,ivc),ci(1),numi)
            call atpbc (s(jas,ivc),h1(1),ci(1),numk,numi,numi)
            call ebtct (si(1),h1(1),c(jas,ivc),numi,numk,numi)
            call fold  (s(ias,ivc),si(1),numi)
 1036    continue
   42 continue
      return
c
   50 if (js.ne.0) go to 60
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
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1038 ivc=1,nvc
            call square(c(jas,ivc),cj(1),numi)
            call apbtct(s(ias,ivc),h1(1),cj(1),numk,numi,numi)
            call atebc (sj(1),h1(1),c(ias,ivc),numi,numk,numi)
            call fold  (s(jas,ivc),sj(1),numi)
 1038    continue
         go to 52
   51    continue
         ias=ia+wtx(mn(ksym+1)  ,is+1)
         jas=ja+wtx(mn(isym+1)+1,js+1)
         intal=1
         call hij   (h1(1),ksym,isym,numk,numi)
         do 1039 ivc=1,nvc
            call square(c(jas,ivc),cj(1),numi)
            call atpbc (s(ias,ivc),h1(1),cj(1),numk,numi,numi)
            call ebtct (sj(1),h1(1),c(ias,ivc),numi,numk,numi)
            call fold  (s(jas,ivc),sj(1),numi)
 1039    continue
   52 continue
      return
c
   60 continue
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
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1040 ivc=1,nvc
            call ambtct(s(jas,ivc),h1(1),c(ias,ivc),numk,numi,numj)
            call atmbc (s(ias,ivc),h1(1),c(jas,ivc),numi,numk,numj)
 1040    continue
         go to 63
   61    continue
         if (ksym.ge.isym) go to 62
         jas=ja+wtx(mn(ksym+1),js+1)
         intal=1
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1041 ivc=1,nvc
            call apbctt(s(jas,ivc),h1(1),c(ias,ivc),numk,numi,numj)
            call atpbct(s(ias,ivc),h1(1),c(jas,ivc),numi,numk,numj)
 1041    continue
         go to 63
   62    continue
         jas=ja+wtx(mn(ksym+1),js+1)
         intal=3
         call hij   (h1(1),ksym,isym,numk,numi)
         do 1042 ivc=1,nvc
            call atpbct(s(jas,ivc),h1(1),c(ias,ivc),numk,numi,numj)
            call apbctt(s(ias,ivc),h1(1),c(jas,ivc),numi,numk,numj)
 1042    continue
c
c     c(i,j) and c(i,k) part
c
   63    continue
         ksym=xor(isym,js)
         numk=numsym(ksym+1)
         if (numk.le.0) go to 66
         if (ksym.ge.jsym) go to 64
         jas=ja+wtx(mn(isym+1),js+1)
         intal=1
         call hij   (h1(1),jsym,ksym,numj,numk)
         do 1043 ivc=1,nvc
            call apbtc (s(jas,ivc),h1(1),c(ias,ivc),numk,numj,numi)
            call apbc  (s(ias,ivc),h1(1),c(jas,ivc),numj,numk,numi)
 1043    continue
         go to 66
   64    continue
         if (ksym.ge.isym) go to 65
         jas=ja+wtx(mn(isym+1),js+1)
         intal=3
         call hij   (h1(1),ksym,jsym,numk,numj)
         do 1044 ivc=1,nvc
            call apbc  (s(jas,ivc),h1(1),c(ias,ivc),numk,numj,numi)
            call apbtc (s(ias,ivc),h1(1),c(jas,ivc),numj,numk,numi)
 1044    continue
         go to 66
   65    continue
         jas=ja+wtx(mn(ksym+1),js+1)
         intal=3
         call hij   (h1(1),ksym,jsym,numk,numj)
         do 1045 ivc=1,nvc
            call atmbc (s(jas,ivc),h1(1),c(ias,ivc),numk,numj,numi)
            call ambtct(s(ias,ivc),h1(1),c(jas,ivc),numj,numk,numi)
 1045    continue
   66 continue
      return
c
c
c
   98 continue
      write (6,99) iseg
   99 format (//,' ***** error in xx entry, iseg=',i5,//)
      stop
c
c********************************** ww *********************************
c
      entry ww
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
         call hiis  (h1(1),isym,numi)
         do 1046 ivc=1,nvc
            call squarw(c(ias,ivc),ci(1),numi)
            call ebc   (sj(1),h1(1),ci(1),numi,numi,numi)
            call foldw (s(ias,ivc),sj(1),numi)
 1046    continue
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
         call hiis  (h1(1),isym,numi)
         call hiis  (h2(1),jsym,numj)
         do 1047 ivc=1,nvc
            call atpbct(s(ias,ivc),h1(1),c(ias,ivc),numi,numi,numj)
            call apbc  (s(ias,ivc),h2(1),c(ias,ivc),numj,numj,numi)
 1047    continue
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
      if (is.ne.js) go to 140
      if (is.ne.0) go to 135
      ias=ia
      jas=ja
      do 131 i=1,nsym
         isym=i-1
         numi=numsym(isym+1)
         if (numi.lt.1) go to 131
         call hii   (h1(1),isym,numi)
         do 1048 ivc=1,nvc
            call squarw(c(ias,ivc),ci(1),numi)
            call squarw(c(jas,ivc),cj(1),numi)
            call ebc   (sj(1),h1(1),ci(1),numi,numi,numi)
            call ebtc  (si(1),h1(1),cj(1),numi,numi,numi)
            call foldw (s(ias,ivc),si(1),numi)
            call foldw (s(jas,ivc),sj(1),numi)
 1048    continue
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
         call hii   (h1(1),isym,numi)
         call hii   (h2(1),jsym,numj)
         do 1049 ivc=1,nvc
            call apbct (s(jas,ivc),c(ias,ivc),h1(1),numj,numi,numi)
            call apbc  (s(jas,ivc),h2(1),c(ias,ivc),numj,numj,numi)
            call apbctt(s(ias,ivc),h1(1),c(jas,ivc),numi,numi,numj)
            call apbtc (s(ias,ivc),h2(1),c(jas,ivc),numj,numj,numi)
 1049    continue
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
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1050 ivc=1,nvc
            call squarw(c(ias,ivc),ci(1),numi)
            call apbtct(s(jas,ivc),h1(1),ci(1),numk,numi,numi)
            call atebc (si(1),h1(1),c(jas,ivc),numi,numk,numi)
            call foldw (s(ias,ivc),si(1),numi)
 1050    continue
         go to 142
  141    continue
         ias=ia+wab(mn(isym+1))
         jas=ja+wtw(mn(ksym+1),js+1)
         intal=3
         call hij   (h1(1),ksym,isym,numk,numi)
         do 1051 ivc=1,nvc
            call squarw(c(ias,ivc),ci(1),numi)
            call atpbc (s(jas,ivc),h1(1),ci(1),numk,numi,numi)
            call ebtct (si(1),h1(1),c(jas,ivc),numi,numk,numi)
            call foldw (s(ias,ivc),si(1),numi)
 1051    continue
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
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1052 ivc=1,nvc
            call squarw(c(jas,ivc),cj(1),numi)
            call apbtct(s(ias,ivc),h1(1),cj(1),numk,numi,numi)
            call atebc (sj(1),h1(1),c(ias,ivc),numi,numk,numi)
            call foldw (s(jas,ivc),sj(1),numi)
 1052    continue
         go to 152
  151    continue
         ias=ia+wtw(mn(ksym+1)  ,is+1)
         jas=ja+wab(mn(isym+1))
         intal=1
         call hij   (h1(1),ksym,isym,numk,numi)
         do 1053 ivc=1,nvc
            call squarw(c(jas,ivc),cj(1),numi)
            call atpbc (s(ias,ivc),h1(1),cj(1),numk,numi,numi)
            call ebtct (sj(1),h1(1),c(ias,ivc),numi,numk,numi)
            call foldw (s(jas,ivc),sj(1),numi)
 1053    continue
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
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1054 ivc=1,nvc
            call apbtct(s(jas,ivc),h1(1),c(ias,ivc),numk,numi,numj)
            call atpbc (s(ias,ivc),h1(1),c(jas,ivc),numi,numk,numj)
 1054    continue
         go to 163
  161    continue
         if (ksym.ge.isym) go to 162
         jas=ja+wtw(mn(ksym+1),js+1)
         intal=1
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1055 ivc=1,nvc
            call apbctt(s(jas,ivc),h1(1),c(ias,ivc),numk,numi,numj)
            call atpbct(s(ias,ivc),h1(1),c(jas,ivc),numi,numk,numj)
 1055    continue
         go to 163
  162    continue
         jas=ja+wtw(mn(ksym+1),js+1)
         intal=3
         call hij   (h1(1),ksym,isym,numk,numi)
         do 1056 ivc=1,nvc
            call atpbct(s(jas,ivc),h1(1),c(ias,ivc),numk,numi,numj)
            call apbctt(s(ias,ivc),h1(1),c(jas,ivc),numi,numk,numj)
 1056    continue
c
c     ----- c(i,j) and c(i,k) part -----
c
  163    continue
         ksym=xor(isym,js)
         numk=numsym(ksym+1)
         if (numk.le.0) go to 166
         if (ksym.ge.jsym) go to 164
         jas=ja+wtw(mn(isym+1),js+1)
         intal=1
         call hij   (h1(1),jsym,ksym,numj,numk)
         do 1057 ivc=1,nvc
            call apbtc (s(jas,ivc),h1(1),c(ias,ivc),numk,numj,numi)
            call apbc  (s(ias,ivc),h1(1),c(jas,ivc),numj,numk,numi)
 1057    continue
         go to 166
  164    continue
         if (ksym.ge.isym) go to 165
         jas=ja+wtw(mn(isym+1),js+1)
         intal=3
         call hij   (h1(1),ksym,jsym,numk,numj)
         do 1058 ivc=1,nvc
            call apbc  (s(jas,ivc),h1(1),c(ias,ivc),numk,numj,numi)
            call apbtc (s(ias,ivc),h1(1),c(jas,ivc),numj,numk,numi)
 1058    continue
         go to 166
  165    continue
         jas=ja+wtw(mn(ksym+1),js+1)
         intal=3
         call hij   (h1(1),ksym,jsym,numk,numj)
         do 1059 ivc=1,nvc
            call atpbc (s(jas,ivc),h1(1),c(ias,ivc),numk,numj,numi)
            call apbtct(s(ias,ivc),h1(1),c(jas,ivc),numj,numk,numi)
 1059    continue
  166 continue
      return
c
c
c
  198 continue
      write (6,199) iseg
  199 format (//,' ***** error in ww entry, iseg=',i5,//)
      stop
c
c******************************* wx ************************************
c
      entry wx
c
      if (iseg.ne.7) go to 298
c
c     ----- iseg 7 case -----
c
      iwx7=iwx7+1
      a=val1*hsqrt3
      b=0.0
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
         call hii   (h1(1),isym,numi)
         do 1060 ivc=1,nvc
            call squarw(c(ias,ivc),ci(1),numi)
            call square(c(jas,ivc),cj(1),numi)
            call ebc   (sj(1),h1(1),ci(1),numi,numi,numi)
            call ebtc  (si(1),h1(1),cj(1),numi,numi,numi)
            call foldw (s(ias,ivc),si(1),numi)
            call fold  (s(jas,ivc),sj(1),numi)
 1060    continue
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
         call hii   (h1(1),isym,numi)
         call hii   (h2(1),jsym,numj)
         do 1061 ivc=1,nvc
            call atpbct(s(jas,ivc),h1(1),c(ias,ivc),numi,numi,numj)
            call ambc  (s(jas,ivc),h2(1),c(ias,ivc),numj,numj,numi)
            call apbctt(s(ias,ivc),h1(1),c(jas,ivc),numi,numi,numj)
            call ambtc (s(ias,ivc),h2(1),c(jas,ivc),numj,numj,numi)
 1061    continue
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
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1062 ivc=1,nvc
            call squarw(c(ias,ivc),ci(1),numi)
            call ambtct(s(jas,ivc),h1(1),ci(1),numk,numi,numi)
            call atembc (si(1),h1(1),c(jas,ivc),numi,numk,numi)
            call foldw (s(ias,ivc),si(1),numi)
 1062    continue
         go to 242
  241    continue
         ias=ia+wab(mn(isym+1))
         jas=ja+wtx(mn(ksym+1),js+1)
         intal=tr1
         call hij   (h1(1),ksym,isym,numk,numi)
         do 1063 ivc=1,nvc
            call squarw(c(ias,ivc),ci(1),numi)
            call atpbc (s(jas,ivc),h1(1),ci(1),numk,numi,numi)
            call ebtct (si(1),h1(1),c(jas,ivc),numi,numk,numi)
            call foldw (s(ias,ivc),si(1),numi)
 1063    continue
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
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1064 ivc=1,nvc
            call square(c(jas,ivc),cj,numi)
            call ambtct(s(ias,ivc),h1(1),cj(1),numk,numi,numi)
            call atembc (sj(1),h1(1),c(ias,ivc),numi,numk,numi)
            call fold  (s(jas,ivc),sj(1),numi)
 1064    continue
         go to 252
  251    continue
         ias=ia+wtw(mn(ksym+1)  ,is+1)
         jas=ja+wtx(mn(isym+1)+1,js+1)
         intal=1
         call hij   (h1(1),ksym,isym,numk,numi)
         do 1065 ivc=1,nvc
            call square(c(jas,ivc),cj(1),numi)
            call atpbc (s(ias,ivc),h1(1),cj(1),numk,numi,numi)
            call ebtct (sj(1),h1(1),c(ias,ivc),numi,numk,numi)
            call fold  (s(jas,ivc),sj(1),numi)
 1065    continue
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
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1066 ivc=1,nvc
            call ambtct(s(jas,ivc),h1(1),c(ias,ivc),numk,numi,numj)
            call atmbc (s(ias,ivc),h1(1),c(jas,ivc),numi,numk,numj)
 1066    continue
         go to 263
  261    continue
         if (ksym.ge.isym) go to 262
         jas=ja+wtx(mn(ksym+1),js+1)
         intal=1
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1067 ivc=1,nvc
            call apbctt(s(jas,ivc),h1(1),c(ias,ivc),numk,numi,numj)
            call atpbct(s(ias,ivc),h1(1),c(jas,ivc),numi,numk,numj)
 1067    continue
         go to 263
  262    continue
         jas=ja+wtx(mn(ksym+1),js+1)
         intal=tr1
         call hij   (h1(1),ksym,isym,numk,numi)
         do 1068 ivc=1,nvc
            call atpbct(s(jas,ivc),h1(1),c(ias,ivc),numk,numi,numj)
            call apbctt(s(ias,ivc),h1(1),c(jas,ivc),numi,numk,numj)
 1068    continue
c
c     ----- c(i,j) and c(i,k) part -----
c
  263    continue
         ksym=xor(isym,js)
         numk=numsym(ksym+1)
         if (numk.le.0) go to 266
         if (ksym.ge.jsym) go to 264
         jas=ja+wtx(mn(isym+1),js+1)
         intal=1
         call hij   (h1(1),jsym,ksym,numj,numk)
         do 1069 ivc=1,nvc
            call ambtc (s(jas,ivc),h1(1),c(ias,ivc),numk,numj,numi)
            call ambc  (s(ias,ivc),h1(1),c(jas,ivc),numj,numk,numi)
 1069    continue
         go to 266
  264    continue
         if (ksym.ge.isym) go to 265
         jas=ja+wtx(mn(isym+1),js+1)
         intal=tr1
         call hij   (h1(1),ksym,jsym,numk,numj)
         do 1070 ivc=1,nvc
            call ambc  (s(jas,ivc),h1(1),c(ias,ivc),numk,numj,numi)
            call ambtc (s(ias,ivc),h1(1),c(jas,ivc),numj,numk,numi)
 1070    continue
         go to 266
  265    continue
         jas=ja+wtx(mn(ksym+1),js+1)
         intal=tr1
         call hij   (h1(1),ksym,jsym,numk,numj)
         do 1071 ivc=1,nvc
            call atpbc (s(jas,ivc),h1(1),c(ias,ivc),numk,numj,numi)
            call apbtct(s(ias,ivc),h1(1),c(jas,ivc),numj,numk,numi)
 1071    continue
  266 continue
      return
c
c
c
  298 continue
      write (6,299) iseg
  299 format (//,' ***** error in wx entry, iseg=',i5,//)
      stop
c
c******************************* xw ************************************
c
      entry xw
c
      if (iseg.ne.9) go to 398
c
c     ----- iseg 9 case -----
c
      ixw9=ixw9+1
      a=val1*hsqrt3
      b=0.0
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
         call hii   (h1(1),isym,numi)
         do 1072 ivc=1,nvc
            call square(c(ias,ivc),ci(1),numi)
            call squarw(c(jas,ivc),cj(1),numi)
            call ebc   (sj(1),h1(1),ci(1),numi,numi,numi)
            call ebtc  (si(1),h1(1),cj(1),numi,numi,numi)
            call fold  (s(ias,ivc),si(1),numi)
            call foldw (s(jas,ivc),sj(1),numi)
 1072    continue
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
         call hii   (h1(1),isym,numi)
         call hii   (h2(1),jsym,numj)
         do 1073 ivc=1,nvc
            call atpbct(s(jas,ivc),h1(1),c(ias,ivc),numi,numi,numj)
            call ambc  (s(jas,ivc),h2(1),c(ias,ivc),numj,numj,numi)
            call apbctt(s(ias,ivc),h1(1),c(jas,ivc),numi,numi,numj)
            call ambtc (s(ias,ivc),h2(1),c(jas,ivc),numj,numj,numi)
 1073    continue
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
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1074 ivc=1,nvc
            call square(c(ias,ivc),ci(1),numi)
            call ambtct(s(jas,ivc),h1(1),ci(1),numk,numi,numi)
            call atembc(si(1),h1(1),c(jas,ivc),numi,numk,numi)
            call fold  (s(ias,ivc),si(1),numi)
 1074    continue
         go to 342
  341    continue
         ias=ia+wtx(mn(isym+1)+1,is+1)
         jas=ja+wtw(mn(ksym+1)  ,js+1)
         intal=tr1
         call hij   (h1(1),ksym,isym,numk,numi)
         do 1075 ivc=1,nvc
            call square(c(ias,ivc),ci(1),numi)
            call atpbc (s(jas,ivc),h1(1),ci(1),numk,numi,numi)
            call ebtct (si(1),h1(1),c(jas,ivc),numi,numk,numi)
            call fold  (s(ias,ivc),si(1),numi)
 1075     continue
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
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1076 ivc=1,nvc
            call squarw(c(jas,ivc),cj(1),numi)
            call ambtct(s(ias,ivc),h1(1),cj(1),numk,numi,numi)
            call atembc(sj(1),h1(1),c(ias,ivc),numi,numk,numi)
            call foldw (s(jas,ivc),sj(1),numi)
 1076    continue
         go to 352
  351    continue
         ias=ia+wtx(mn(ksym+1)  ,is+1)
         jas=ja+wab(mn(isym+1))
         intal=1
         call hij   (h1(1),ksym,isym,numk,numi)
         do 1077 ivc=1,nvc
            call squarw(c(jas,ivc),cj(1),numi)
            call atpbc (s(ias,ivc),h1(1),cj(1),numk,numi,numi)
            call ebtct (sj(1),h1(1),c(ias,ivc),numi,numk,numi)
            call foldw (s(jas,ivc),sj(1),numi)
 1077    continue
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
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1078 ivc=1,nvc
            call apbtct(s(jas,ivc),h1(1),c(ias,ivc),numk,numi,numj)
            call atpbc (s(ias,ivc),h1(1),c(jas,ivc),numi,numk,numj)
 1078    continue
         go to 363
  361    continue
         if (ksym.ge.isym) go to 362
         jas=ja+wtw(mn(ksym+1),js+1)
         intal=1
         call hij   (h1(1),isym,ksym,numi,numk)
         do 1079 ivc=1,nvc
            call apbctt(s(jas,ivc),h1(1),c(ias,ivc),numk,numi,numj)
            call atpbct(s(ias,ivc),h1(1),c(jas,ivc),numi,numk,numj)
 1079    continue
         go to 363
  362    continue
         jas=ja+wtw(mn(ksym+1),js+1)
         intal=tr1
         call hij   (h1(1),ksym,isym,numk,numi)
         do 1080 ivc=1,nvc
            call atpbct(s(jas,ivc),h1(1),c(ias,ivc),numk,numi,numj)
            call apbctt(s(ias,ivc),h1(1),c(jas,ivc),numi,numk,numj)
 1080    continue
c
c     ----- c(i,j) and c(i,k) part -----
c
  363    continue
         ksym=xor(isym,js)
         numk=numsym(ksym+1)
         if (numk.le.0) go to 366
         if (ksym.ge.jsym) go to 364
         jas=ja+wtw(mn(isym+1),js+1)
         intal=1
         call hij   (h1(1),jsym,ksym,numj,numk)
         do 1081 ivc=1,nvc
            call ambtc (s(jas,ivc),h1(1),c(ias,ivc),numk,numj,numi)
            call ambc  (s(ias,ivc),h1(1),c(jas,ivc),numj,numk,numi)
 1081    continue
         go to 366
  364    continue
         if (ksym.ge.isym) go to 365
         jas=ja+wtw(mn(isym+1),js+1)
         intal=tr1
         call hij   (h1(1),ksym,jsym,numk,numj)
         do 1082 ivc=1,nvc
            call ambc  (s(jas,ivc),h1(1),c(ias,ivc),numk,numj,numi)
            call ambtc (s(ias,ivc),h1(1),c(jas,ivc),numj,numk,numi)
 1082    continue
         go to 366
  365    continue
         jas=ja+wtw(mn(ksym+1),js+1)
         intal=tr1
         call hij   (h1(1),ksym,jsym,numk,numj)
         do 1083 ivc=1,nvc
            call atmbc(s(jas,ivc),h1(1),c(ias,ivc),numk,numj,numi)
            call ambtct(s(ias,ivc),h1(1),c(jas,ivc),numj,numk,numi)
 1083    continue
  366 continue
      return
c
c
c
  398 continue
      write (6,399) iseg
  399 format (//,' ***** error in xw entry, iseg=',i5,//)
      stop
c
      end
