*deck @(#)initdg.f	5.1  11/6/94
      subroutine initdg(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,
     1                  ss)
c
c
c
      implicit real*8 (a-h,o-z)
      integer xor
      integer arr,tr1,tr2,os,asm,wtw,wtx,wty,wab,ss,symorb,aos
      integer bmax,orbfrm,ssj
      real*8 int(nmax),c(nwks),s(nwks),int1(2)
      dimension ijadd(numij),kadd(symorb),ladd(symorb),wtw(orbfrm,nsym)
      dimension wtx(orbfrm,nsym),wty(orbfrm),wab(orbfrm),ss(norbs)
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,nrow4o,nwksoc,nlevoc,norboc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
c  universal designation for these commons local variant follows
c     common /all/ val1,val2,val3,arr,itr1,itr2,ia,ja,itype,isegt
c    *,lvfrm1,nlwki,nlwkj,imax,imin
c     common /symm/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
c    #,             numsym(8)
      common /tapes/itap20,itape5,itape6,itape8,itap12,itap03,itap04
     *,             itape3,itap05,itap06
      common /all/val1,val2,val3,arr,tr1,tr2,ia,ja,m,iseg,n,n1,n2
     *,           imax,imin
      common /symm/ asm,js,is,mx(8),mn(8),os(8),numsym(8)
c:
      common /minmax/ iming,imaxg,jming,jmaxg
      common/dry1/idry,jdry,ksegsv,indxs(100)
c:
c
      save sqrt2,sqt1p5
c
      sqrt2=sqrt(2.0d+00)
      sqt1p5=sqrt(1.5d+00)
c
      do 5 i=1,nsym
      mn(i)=1000
    5 mx(i)=0
      mn(1)=1
      ism=ss(1)
      do 20 i=2,n
      if(ss(i).eq.ism) go to 20
      mx(ism)=i-1
      ism=ss(i)
      mn(ism)=i
   20 continue
      mx(ism)=n
c     write(itape6,30) (mn(j),j=1,nsym)
   30 format(' minsym=',8i5)
c     write(itape6,35) (mx(j),j=1,nsym)
   35 format(' maxsym=',8i5)
      return
c
c
c
      entry diagonal(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,
     $                  ss)
c..rlm     entry diagonal
c
      aos=os(asm+1)
      ia=ja
      go to (314,305,304,315,316,317,318,321),m
      write(itape6,10) m
   10 format(1x,'unknown type value--',i5)
      call lnkerr(' m911: error ')
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
c     z=val4*(int(lad+1)-2.*int(lad+2))
c     d(jj)=d(jj)+z
      z=val4*c(jj)*s(jj)
      int(lad+1)=int(lad+1)+z
      int(lad+2)=int(lad+2)-2*z
c
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
c     z=val4*int(lad+1)
c     d(jj)=d(jj)+z
      int(lad+1)=int(lad+1)+val4*c(jj)*s(jj)
c
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
c     z=xx*int(lad+1)+val5*int(lad+2)
c     d(jj)=d(jj)+z
      z=c(jj)*s(jj)
      int(lad+1)=int(lad+1)+xx*z
      int(lad+2)=int(lad+2)+val5*z
c
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
c     z=val5*int(lad+2)
c     d(jj)=d(jj)+z
      int(lad+2)=int(lad+2)+val5*c(jj)*s(jj)
c
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
c     z=val4*(int(lad+1)-2.*int(lad+2))
c     d(jj)=d(jj)+z
      z=val4*c(jj)*s(jj)
      int(lad+1)=int(lad+1)+z
      int(lad+2)=int(lad+2)-z*2
c
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
c     z1=val4*(int(lad1+1)-2.0d+00*int(lad1+2))
      z1=0.0d+00
      jj=ja+wtw(i,1+js)+wty(m2)
      do 1407 j=m2,m1
      lsm=xor(asm,ss(j))
      los=os(lsm+1)
      lad=arr+kadd(j+aos)+ladd(j+los)
c     z=val4*(int(lad+1)-2.*int(lad+2))+z1
c     d(jj)=d(jj)+z
      z=val4*c(jj)*s(jj)
      z1=z1+z
      int(lad+1)=int(lad+1)+z
      int(lad+2)=int(lad+2)-z*2
c
 1407 jj=jj+1
c
      int(lad1+1)=int(lad1+1)+z1
      int(lad1+2)=int(lad1+2)-z1*2
c
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
c     z1=val4*(int(lad1+1)-2.0d+00*int(lad1+2))
      z1=0.0d+00
      do 502 j=m2,m1
      lsm=xor(asm,ss(j))
      los=os(lsm+1)
      lad=arr+kadd(j+aos)+ladd(j+los)
c     z=val4*(int(lad+1)-2.*int(lad+2))+z1
c     d(jj)=d(jj)+z
      z=val4*c(jj)*s(jj)
      z1=z1+z
      int(lad+1)=int(lad+1)+z
      int(lad+2)=int(lad+2)-z*2
c
  502 jj=jj+1
c
      int(lad1+1)=int(lad1+1)+z1
      int(lad1+2)=int(lad1+2)-z1*2
c
  501 continue
      return
  305 if (iseg.eq.4) go to 3051
      if (iseg.eq.5) go to 3052
      if (iseg.eq.8) go to 3054
      if (iseg.ne.6) go to 999
c   xx   no  11   iseg=6,5
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
c     z1=xx*int(lad1+1)+val5*int(lad1+2)
      z1=0.0d+00
      do 561 j=m2,m1
      lsm=xor(asm,ss(j))
      los=os(lsm+1)
      lad=arr+kadd(j+aos)+ladd(j+los)
c     z=xx*int(lad+1)+val5*int(lad+2)+z1
c     d(jj)=d(jj)+z
      z=c(jj)*s(jj)
      z1=z1+z
      int(lad+1)=int(lad+1)+xx*z
      int(lad+2)=int(lad+2)+val5*z
c
  561 jj=jj+1
c
      int(lad1+1)=int(lad1+1)+xx*z1
      int(lad1+2)=int(lad1+2)+val5*z1
c
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
c     z1=val5*int(lad1+2)
      z1=0.0d+00
      do 563 j=m2,m1
      lsm=xor(asm,ss(j))
      los=os(lsm+1)
      lad=arr+kadd(j+aos)+ladd(j+los)
c     z=val5*int(lad+2)+z1
c     d(jj)=d(jj)+z
      z=val5*c(jj)*s(jj)
      z1=z1+z
      int(lad+2)=int(lad+2)+z
c
  563 jj=jj+1
c
      int(lad1+2)=int(lad1+2)+z1
c
  562 continue
      return
c   xx   no 11  four externals
  316 jja=ja
      m5=imax
      m6=imin
      if (m6.lt.2) m6=2
      do 2002 i=m6,m5
      if (i.gt.imaxg.or.i.lt.iming.or.i.gt.jmaxg.or.i.lt.jming)goto 2002
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
c     z1=int(lad2+2)
      z1=0.0d+00
      do 2003 j=m2,m1
c     los=os(ss(j)+1)
c:c   lad1=lad+kadd(j+kos)+ladd(j+los)
c     d(jj)=d(jj)+                           z1
      z1=z1+c(jj)*s(jj)
 2003 jj=jj+1
      int(lad2+2)=int(lad2+2)+z1
      iix=indxs(i)+i
      int1(iix)=int1(iix)+z1
c
 2002 continue
      jja=ja
      m5=imax
      m6=imin
      if (m6.lt.2) m6=2
      do 8002 i=m6,m5
      if (i.gt.imaxg.or.i.lt.iming.or.i.gt.jmaxg.or.i.lt.jming)goto 8002
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
c     d(jj)=d(jj)+(-int(lad1+1)+int(lad1+2))
      z=c(jj)*s(jj)
      int(lad1+1)=int(lad1+1)-z
      int(lad1+2)=int(lad1+2)+z
c
 8003 jj=jj+1
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
      if (j.gt.imaxg.or.j.lt.iming.or.j.gt.jmaxg.or.j.lt.jming)goto 2083
      los=os(ss(j)+1)
      lad=ijadd(j*(j+1)/2)+kadd(j+kos)+ladd(j+los)
c     z=int(lad+2)
c     d(jj)=d(jj)+z
      int(lad+2)=int(lad+2)+c(jj)*s(jj)
 2083 jj=jj+1
 2082 continue
      return
c   ww   no   68  four externals
  317 continue
      m5=imax
      m6=imin
      if(is.ne.0) go to 6601
      do 2030 i=m6,m5
      if (i.gt.imaxg.or.i.lt.iming.or.i.gt.jmaxg.or.i.lt.jming)goto 2030
      jj=ja+wab(i)
      kos=os(1)
      los=os(ss(i)+1)
      lad=ijadd(i*(i+1)/2)+kadd(i+kos)+ladd(i+los)
c     z=int(lad+2)+int(lad+2)+int(lad+1)
c     d(jj)=d(jj)+z
      int(lad+2)=int(lad+2)+2*c(jj)*s(jj)
      int(lad+1)=int(lad+1)+  c(jj)*s(jj)
cdry      7/25/84
      iix=indxs(i)+i
      int1(iix)=int1(iix)+2*c(jj)*s(jj)
c
 2030 continue
c   no  81
 6601 m6=imin
      if (m6.lt.2) m6=2
      do 2046 i=m6,m5
      if (i.gt.imaxg.or.i.lt.iming.or.i.gt.jmaxg.or.i.lt.jming)goto 2046
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
c     z1=int(lad2+2)
      z1=0.0d+00
      do 2047 j=m2,m1
c:    los=os(ss(j)+1)
c:    lad1=lad+kadd(j+kos)+ladd(j+los)
c:    z=int(lad1+1)+int(lad1+2)+z1
c     d(jj)=d(jj)+z1
      z1=z1+c(jj)*s(jj)
 2047 jj=jj+1
c
      int(lad2+2)=int(lad2+2)+z1
      iix=indxs(i)+i
      int1(iix)=int1(iix)+z1
c
 2046 continue
      m6=imin
      if (m6.lt.2) m6=2
      do 8046 i=m6,m5
      if (i.gt.imaxg.or.i.lt.iming.or.i.gt.jmaxg.or.i.lt.jming)goto 8046
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
c     z=int(lad1+1)+int(lad1+2)
c     d(jj)=d(jj)+z
      z=c(jj)*s(jj)
      int(lad1+1)=int(lad1+1)+z
      int(lad1+2)=int(lad1+2)+z
c
 8047 jj=jj+1
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
      if (j.gt.imaxg.or.j.lt.iming.or.j.gt.jmaxg.or.j.lt.jming)goto 2085
      los=os(ss(j)+1)
      lad=ijadd(j*(j+1)/2)+kadd(j+kos)+ladd(j+los)
c     z=int(lad+2)
c     d(jj)=d(jj)+z
      int(lad+2)=int(lad+2)+c(jj)*s(jj)
      iix=indxs(j)+j
      int1(iix)=int1(iix)+c(jj)*s(jj)
c
 2085 jj=jj+1
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
      if (i.gt.imaxg.or.i.lt.iming.or.i.gt.jmaxg.or.i.lt.jming)goto 2080
      kos=os(1)
      los=os(ss(i)+1)
      lad=ijadd(i*(i+1)/2)+kadd(i+kos)+ladd(i+los)+2
c     d(jj)=d(jj)+int(lad)
      int(lad)=int(lad)+c(jj)*s(jj)
      iix=indxs(i)+i
      int1(iix)=int1(iix)+c(jj)*s(jj)
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
c     z1=val1*int(lad1)
      z1=0.0d+00
      do 511 j=m2,m1
      lsm=xor(asm,ss(j))
      los=os(lsm+1)
      lad=arr+kadd(j+aos)+ladd(j+los)
c     z=val1*int(lad+1)+z1
c     d(jj)=d(jj)+z
      z=val1*c(jj)*s(jj)
      z1=z1+z
      int(lad+1)=int(lad+1)+z
c
  511 jj=jj+1
c
      int(lad1)=int(lad1)+z1
c
  510 continue
      return
c
c     ----- closed internal loop -----
c
  318 continue
c     z=dot(d(ja),d(ja),n1)*val1
      z=dot(c(ja),s(ja),n1)*val1
      int(arr+tr1)=int(arr+tr1)+z
      int(arr+tr2)=int(arr+tr2)+z*val2
      if(ksegsv.gt.15)return
      if((idry.ne.jdry).or.(tr1.ne.2))then
      write(itape6,50501)
50501 format(' idry jdry tr2 ',3i6)
      stop ' idry ne jdry   or  tr1 ne 2  in initdg '
      else
      iix=indxs(idry)+jdry
      int1(iix)=int1(iix)+z
      endif
      return
c
c     ----- zz arrival is an error -----
c
  321 continue
  999 write (itape6,998) m,iseg
  998 format (2i4)
      return
      end
