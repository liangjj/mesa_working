*deck @(#)initdg.f	5.1  11/6/94
      subroutine initdg(int,d,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      implicit real*8 (a-h,o-z)
c
      integer xor
      integer arr,tr1,tr2,os,asm,wtw,wtx,wty,wab,ss,symorb,aos
      integer bmax,orbfrm,ssj
      dimension ijadd(numij),kadd(symorb),ladd(symorb),wtw(orbfrm,nsym)
      dimension wtx(orbfrm,nsym),wty(orbfrm),wab(orbfrm),ss(norbs)
c
      real*8 int(nmax),d(nwks)
c
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
c     universal designation for these commons local variant follows
c     common /all/ val1,val2,val3,arr,itr1,itr2,ia,ja,itype,isegt
c    *,lvfrm1,nlwki,nlwkj,imax,imin
c     common /symq/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
c    #,             numsym(8)
      common /io/     itape5,itape6
      common /all/val1,val2,val3,arr,tr1,tr2,ia,ja,m,iseg,n,n1,n2
     *,           imax,imin
      common /symq/ asm,js,is,mx(8),mn(8),os(8),numsym(8)
c:
      common /minmax/ iming,imaxg,jming,jmaxg
c:
      save sqrt2,sqt1p5
c
      sqrt2=sqrt(2.0d+00)
      sqt1p5=sqrt(1.5d+00)
cps      do 5 i=1,nsym
cps         mn(i)=1000
cps    5 mx(i)=0
cps      ism=ss(1)
cps      mn(ism)=1
cps      do 20 i=2,n
cps         if(ss(i).eq.ism) go to 20
cps         mx(ism)=i-1
cps         ism=ss(i)
cps         mn(ism)=i
cps   20 continue
cps      mx(ism)=n
      return
c
      entry diagonal(int,d,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
c
      aos=os(asm+1)
      ia=ja
      go to (314,305,304,315,316,317,318,321),m
      write(itape6,10) m
   10 format(1x,'unknown type value--',i5)
      stop
c    yy   no   10   iseg=4
  304 if (iseg.eq.5) go to 3043
      if (iseg.eq.6) go to 3043
      if (iseg.eq.8) go to 3042
      if (iseg.ne.4) go to 999
      val4=-val1/sqrt2
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.eq.1000) return
      if (m1.gt.n) m1=n
      if (m1.lt.m2) return
      ii=ia+wty(m2)
      jj=ii
      do 460 i=m2,m1
         lsm=xor(asm,ss(i))
         los=os(lsm+1)
         lad=arr+ladd(i+los)+kadd(i+aos)
         z=val4*(int(lad+1)-2.*int(lad+2))
         d(jj)=d(jj)+z
  460 jj=jj+1
      return
c  yy   no 10 iseg=8
 3042 val4=val1*sqt1p5
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.eq.1000) return
      if (m1.gt.n) m1=n
      if (m1.lt.m2) return
      jj=ja+wty(m2)
      do 430 i=m2,m1
         lsm=xor(asm,ss(i))
         los=os(lsm+1)
         lad=arr+ladd(i+los)+kadd(i+aos)
         z=val4*int(lad+1)
         d(jj)=d(jj)+z
  430 jj=jj+1
      return
c  yy  no 10  iseg=5,6
 3043 continue
      xx=-val1/sqrt2+val3*sqt1p5
      if (abs(xx).lt.1.0d-06) go to 490
      val5=sqrt2*val1
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.eq.1000) return
      if (m1.gt.n) m1=n
      if (m1.lt.m2) return
      jj=ja+wty(m2)
      do 493 i=m2,m1
         lsm=xor(asm,ss(i))
         los=os(lsm+1)
         lad=arr+ladd(i+los)+kadd(i+aos)
         z=xx*int(lad+1)+val5*int(lad+2)
         d(jj)=d(jj)+z
  493 jj=jj+1
      return
  490 val5=sqrt2*val1
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.eq.1000) return
      if (m1.gt.n) m1=n
      if (m1.gt.m2) return
      jj=ja+wty(m2)
      do 494 i=m2,m1
         lsm=xor(asm,ss(i))
         los=os(lsm+1)
         lad=arr+ladd(i+los)+kadd(i+aos)
         z=val5*int(lad+2)
         d(jj)=d(jj)+z
  494 jj=jj+1
      return
c   ww    no 68  iseg=4,5,6
  314 if (iseg.eq.8) return
      if (iseg.eq.6) go to 299
      if (iseg.eq.4.or.iseg.eq.5) go to 299
      go to 999
  299 if (js.ne.0) go to 6600
      val4=-val1*sqrt2
      do 1401 i=1,n
         jj=ja+wab(i)
         lsm=xor(asm,ss(i))
         los=os(lsm+1)
         lad=kadd(i+aos)+ladd(i+los)+arr
         z=val4*(int(lad+1)-2.*int(lad+2))
         d(jj)=d(jj)+z
 1401 continue
c   no 81
 6600 val4=-val1/sqrt2
      do 1406 i=2,n
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i)) go to 1406
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.eq.1000) go to 1406
         if (ssj.eq.ss(i)) m1=i-1
         if (m1.lt.m2) go to 1406
         lsm=xor(asm,ss(i))
         los=os(lsm+1)
         lad1=arr+kadd(i+aos)+ladd(i+los)
         z1=val4*(int(lad1+1)-2.0d+00*int(lad1+2))
         jj=ja+wtw(i,1+js)+wty(m2)
         do 1407 j=m2,m1
            lsm=xor(asm,ss(j))
            los=os(lsm+1)
            lad=arr+kadd(j+aos)+ladd(j+los)
            z=val4*(int(lad+1)-2.*int(lad+2))+z1
            d(jj)=d(jj)+z
 1407    jj=jj+1
 1406 continue
      return
c    xx   no   11   iseg=4
 3051 continue
      val4=-val1/sqrt2
      do 501 i=2,n
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i)) go to 501
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.eq.1000) go to 501
         if (ssj.eq.ss(i)) m1=i-1
         if (m1.lt.m2) go to 501
         jj=ja+wtx(i,1+js)+wty(m2)
         lsm=xor(asm,ss(i))
         los=os(lsm+1)
         lad1=arr+kadd(i+aos)+ladd(i+los)
         z1=val4*(int(lad1+1)-2.0d+00*int(lad1+2))
         do 502 j=m2,m1
            lsm=xor(asm,ss(j))
            los=os(lsm+1)
            lad=arr+kadd(j+aos)+ladd(j+los)
            z=val4*(int(lad+1)-2.*int(lad+2))+z1
            d(jj)=d(jj)+z
  502    jj=jj+1
  501 continue
      return
  305 if (iseg.eq.4) go to 3051
      if (iseg.eq.5) go to 3052
      if (iseg.eq.8) go to 3054
      if (iseg.ne.6) go to 999
c xx   no  11   iseg=6,5
 3052 xx=-val1/sqrt2+val3
      val5=sqrt2*val1
      if (abs(xx).lt.1.0d-06) go to 550
      do 560 i=2,n
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i)) go to 560
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.eq.1000) go to 560
         if (ssj.eq.ss(i)) m1=i-1
         if (m1.lt.m2) go to 560
         jj=ja+wtx(i,1+js)+wty(m2)
         lsm=xor(asm,ss(i))
         los=os(lsm+1)
         lad1=arr+kadd(i+aos)+ladd(i+los)
         z1=xx*int(lad1+1)+val5*int(lad1+2)
         do 561 j=m2,m1
            lsm=xor(asm,ss(j))
            los=os(lsm+1)
            lad=arr+kadd(j+aos)+ladd(j+los)
            z=xx*int(lad+1)+val5*int(lad+2)+z1
            d(jj)=d(jj)+z
  561    jj=jj+1
  560 continue
      return
  550 do 562 i=2,n
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i)) go to 562
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.eq.1000) go to 562
         if (ssj.eq.ss(i)) m1=i-1
         if (m1.lt.m2) go to 562
         jj=ja+wtx(i,1+js)+wty(m2)
         lsm=xor(asm,ss(i))
         los=os(lsm+1)
         lad1=arr+kadd(i+aos)+ladd(i+los)
         z1=val5*int(lad1+2)
         do 563 j=m2,m1
            lsm=xor(asm,ss(j))
            los=os(lsm+1)
            lad=arr+kadd(j+aos)+ladd(j+los)
            z=val5*int(lad+2)+z1
            d(jj)=d(jj)+z
  563    jj=jj+1
  562 continue
      return
c   xx   no 11  four externals
  316 jja=ja
      m5=imax
      m6=imin
      if (m6.lt.2) m6=2
      do 2002 i=m6,m5
         if (i.gt.imaxg.or.i.lt.iming.or.i.gt.jmaxg.or.i.lt.jming)goto
     *   2002
         lad=ijadd(i*(i+1)/2)
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i)) go to 2002
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.eq.1000) go to 2002
         if (ssj.eq.ss(i)) m1=i-1
         if (m1.lt.m2) go to 2002
         jj=ja+wtx(i,1+js)+wty(m2)
         kos=os(1)
         los=os(ss(i)+1)
         lad2=lad+kadd(i+kos)+ladd(i+los)
         z1=int(lad2+2)
         do 2003 j=m2,m1
            los=os(ss(j)+1)
c:    lad1=lad+kadd(j+kos)+ladd(j+los)
            d(jj)=d(jj)+                           z1
 2003    jj=jj+1
 2002 continue
      jja=ja
      m5=imax
      m6=imin
      if (m6.lt.2) m6=2
      do 8002 i=m6,m5
         if (i.gt.imaxg.or.i.lt.iming.or.i.gt.jmaxg.or.i.lt.jming)goto
     *   8002
         lad=ijadd(i*(i+1)/2)
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i)) go to 8002
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.eq.1000) go to 8002
         if (ssj.eq.ss(i)) m1=i-1
         if (m1.lt.m2) go to 8002
         jj=ja+wtx(i,1+js)+wty(m2)
         kos=os(1)
         los=os(ss(i)+1)
c:    lad2=lad+kadd(i+kos)+ladd(i+los)
c:    z1=int(lad2+2)
         do 8003 j=m2,m1
            los=os(ss(j)+1)
            lad1=lad+kadd(j+kos)+ladd(j+los)
            d(jj)=d(jj)+(-int(lad1+1)+int(lad1+2))
 8003    jj=jj+1
 8002 continue
      do 2082 i=2,n
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i))go to 2082
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.eq.1000) go to 2082
         if (ssj.eq.ss(i)) m1=i-1
         if (m1.lt.m2) go to 2082
         m5=imax
         m6=imin
         if (m5.gt.m1) m5=m1
         if (m6.lt.m2) m6=m2
         jj=ja+wtx(i,1+js)+wty(m6)
         kos=os(1)
         do 2083 j=m6,m5
            if (j.gt.imaxg.or.j.lt.iming.or.j.gt.jmaxg.or.j.lt.jming)
     *      goto 2083
            los=os(ss(j)+1)
            lad=ijadd(j*(j+1)/2)+kadd(j+kos)+ladd(j+los)
            z=int(lad+2)
            d(jj)=d(jj)+z
 2083    jj=jj+1
 2082 continue
      return
c   ww   no   68  four externals
  317 continue
      m5=imax
      m6=imin
      if(is.ne.0) go to 6601
      do 2030 i=m6,m5
         if (i.gt.imaxg.or.i.lt.iming.or.i.gt.jmaxg.or.i.lt.jming)goto
     *   2030
         jj=ja+wab(i)
         kos=os(1)
         los=os(ss(i)+1)
         lad=ijadd(i*(i+1)/2)+kadd(i+kos)+ladd(i+los)
         z=int(lad+2)+int(lad+2)+int(lad+1)
         d(jj)=d(jj)+z
 2030 continue
c   no  81
 6601 m6=imin
      if (m6.lt.2) m6=2
      do 2046 i=m6,m5
         if (i.gt.imaxg.or.i.lt.iming.or.i.gt.jmaxg.or.i.lt.jming)goto
     *   2046
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i)) go to 2046
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.eq.1000) go to 2046
         if (ssj.eq.ss(i)) m1=i-1
         if (m1.lt.m2) go to 2046
         jj=ja+wtw(i,1+js)+wty(m2)
         lad=ijadd(i*(i+1)/2)
         kos=os(1)
         los=os(ss(i)+1)
         lad2=lad+kadd(i+kos)+ladd(i+los)
         z1=int(lad2+2)
         do 2047 j=m2,m1
c:    los=os(ss(j)+1)
c:    lad1=lad+kadd(j+kos)+ladd(j+los)
c:    z=int(lad1+1)+int(lad1+2)+z1
            d(jj)=d(jj)+z1
 2047    jj=jj+1
 2046 continue
      m6=imin
      if (m6.lt.2) m6=2
      do 8046 i=m6,m5
         if (i.gt.imaxg.or.i.lt.iming.or.i.gt.jmaxg.or.i.lt.jming)goto
     *   8046
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i)) go to 8046
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.eq.1000) go to 8046
         if (ssj.eq.ss(i)) m1=i-1
         if (m1.lt.m2) go to 8046
         jj=ja+wtw(i,1+js)+wty(m2)
         lad=ijadd(i*(i+1)/2)
         kos=os(1)
         los=os(ss(i)+1)
c:    lad2=lad+kadd(i+kos)+ladd(i+los)
c:    z1=int(lad2+2)
         do 8047 j=m2,m1
            los=os(ss(j)+1)
            lad1=lad+kadd(j+kos)+ladd(j+los)
            z=int(lad1+1)+int(lad1+2)
            d(jj)=d(jj)+z
 8047    jj=jj+1
 8046 continue
      do 2084 i=2,n
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i)) go to 2084
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.eq.1000) go to 2084
         if (ssj.eq.ss(i)) m1=i-1
         if (m1.lt.m2) go to 2084
         m5=imax
         m6=imin
         if (m5.gt.m1) m5=m1
         if (m6.lt.m2) m6=m2
         jj=ja+wtw(i,1+js)+wty(m6)
         kos=os(1)
         do 2085 j=m6,m5
            if (j.gt.imaxg.or.j.lt.iming.or.j.gt.jmaxg.or.j.lt.jming)
     *      goto 2085
            los=os(ss(j)+1)
            lad=ijadd(j*(j+1)/2)+kadd(j+kos)+ladd(j+los)
            z=int(lad+2)
            d(jj)=d(jj)+z
 2085    jj=jj+1
 2084 continue
      return
c   yy   no 10 four externals
  315 continue
      m1=mx(js+1)
      m2=mn(js+1)
      if (m2.eq.1000) return
      if (m1.gt.n) m1=n
      if (m1.lt.m2) return
      m5=imax
      m6=imin
      if (m5.gt.m1) m5=m1
      if (m6.lt.m2) m6=m2
      jj=ja+wty(m6)
      do 2080 i=m6,m5
         if (i.gt.imaxg.or.i.lt.iming.or.i.gt.jmaxg.or.i.lt.jming)goto
     *   2080
         kos=os(1)
         los=os(ss(i)+1)
         lad=ijadd(i*(i+1)/2)+kadd(i+kos)+ladd(i+los)+2
         d(jj)=d(jj)+int(lad)
 2080 jj=jj+1
      return
c   xx    no  11  iseg=8
 3054 continue
      do 510 i=2,n
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i)) go to 510
         m1=mx(ssj+1)
         m2=mn(ssj+1)
         if (m2.eq.1000) go to 510
         if (ssj.eq.ss(i)) m1=i-1
         if (m1.lt.m2) go to 510
         jj=ja+wtx(i,1+js)+wty(m2)
         lsm=xor(asm,ss(i))
         los=os(lsm+1)
         lad1=arr+kadd(i+aos)+ladd(i+los)+1
         z1=val1*int(lad1)
         do 511 j=m2,m1
            lsm=xor(asm,ss(j))
            los=os(lsm+1)
            lad=arr+kadd(j+aos)+ladd(j+los)
            z=val1*int(lad+1)+z1
            d(jj)=d(jj)+z
  511    jj=jj+1
  510 continue
      return
c   closed internal loop
  318 z=val1*(int(arr+tr1)+val2*int(arr+tr2))
      jj=ja
      do 4000 i=1,n1
         d(jj)=d(jj)+z
 4000 jj=jj+1
      return
c zz   arrival is an error
  321 continue
  999 write (6,998) m,iseg
  998 format (2i4)
      return
      end
