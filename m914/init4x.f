*deck @(#)init4x.f	1.1  11/30/90
      subroutine init4x(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
c
c
c
      implicit real*8 (a-h,o-z)
c
cvax  extended dummy int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss
c
      integer xor
      integer arr,tr1,tr2,asm,aos,os,wtw,wtx,wty,wab,ss,ssi,ssj,symorb
      integer bmax,orbfrm
      real*8 int(nmax),c(nwks,mxvc),s(nwks,mxvc)
      dimension kadd(symorb),ladd(symorb),ijadd(numij),wtw(orbfrm,nsym)
      dimension wtx(orbfrm,nsym),wty(orbfrm),wab(orbfrm),ss(norbs)
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
      common /minmax/ iming,imaxg,jming,jmaxg
      common /nvectr/ nvc,mxvc
c
      real*8 z,val1,val2,val3
      sqrt2=sqrt(2.0d+00)
      sqrt3=sqrt(3.0d+00)
      sqt1p5=sqrt(1.5d+00)
      return
c
      entry shape4
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
      do 9000 ivc=1,nvc
         jj=ja+wty(m6)
         do 2000 i=m6,m5
c
            if (i.gt.jmaxg.or.i.lt.jming) go to 2000
c
            if (is.eq.js) m3=i-1
            if (m3.lt.m4) go to 2000
            kos=os(1)
            los=os(ss(i)+1)
            lad=ijadd(i*(i+1)/2)+kadd(i+kos)+3+ladd(m4+los)
            ii=ia+wty(m4)
            cjj=c(jj,ivc)
            sjj=0.0
            do 2001 j=m4,m3
               z=int(lad)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               ii=ii+1
               lad=lad+3
 2001       continue
            s(jj,ivc)=s(jj,ivc)+sjj
 2000    jj=jj+1
 9000 continue
c
      return
c
c   xx   no 13 & 161 four externals
c***********************************************************************
  315 continue
      do 9001 ivc=1,nvc
         m5=imax
         m6=imin
         if (m6.lt.3) m6=3
         if (m5.lt.m6) go to 2003
         do 2004 i=m6,m5
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
               lad1=lad+kadd(j+kos)+ladd(m4+los)
               ii=iia+wty(m4)
               cjj=c(jj,ivc)
               sjj=0.0
               do 2006 k=m4,m3
                  z=-int(lad1+1)+int(lad1+2)
                  s(ii,ivc)=s(ii,ivc)+z*cjj
                  sjj=sjj+z*c(ii,ivc)
                  ii=ii+1
                  lad1=lad1+3
 2006          continue
               s(jj,ivc)=s(jj,ivc)+sjj
 2005       jj=jj+1
 2004    continue
 2003    continue
c
c     no 161
c
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
               if (ssi.eq.ssj) m3=j-1
               if (m3.lt.m4) go to 7085
               los=os(ss(j)+1)
               lad3=ijadd(j*(j+1)/2)+kadd(j+kos)+3+ladd(m4+los)
               ii=iia+wty(m4)
               cjj=c(jj,ivc)
               sjj=0.0
               do 7086 k=m4,m3
                  z=int(lad3)
                  s(ii,ivc)=s(ii,ivc)+z*cjj
                  sjj=sjj+z*c(ii,ivc)
                  ii=ii+1
                  lad3=lad3+3
 7086          continue
               s(jj,ivc)=s(jj,ivc)+sjj
 7085       jj=jj+1
 7084    continue
 8084    continue
c
c   no  14
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
            lad3=ijadd(i*(i+1)/2)+kadd(i+kos)+3
            do 2008 j=m2,m1
               if (ssi.eq.ssj) m3=j-1
               if (m3.lt.m4) go to 2008
               z1=int(lad3+ladd(j+los1))
               ksm=xor(ss(i),ss(j))
               kos=os(ksm+1)
               lad1=ijadd(i*(i-1)/2+j)
               ii=ia+wtx(j,is+1)+wty(m4)
               jj=jja+wty(m4)
cdir$ ivdep
               do 2009 k=m4,m3
                  lsm=xor(ksm,ss(k))
                  los=os(lsm+1)
                  lad2=lad1+kadd(k+kos)+ladd(k+los)
                  z=-int(lad2+1)+int(lad2+2)+z1
                  s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
                  s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
                  ii=ii+1
                  jj=jj+1
 2009          continue
 2008       continue
 2007    continue
 8007    continue
c
c    no 16   both
c
         if (m5.lt.m6) go to 8010
         do 2010 i=m6,m5
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
               cjj=c(jj,ivc)
               sjj=0.0
               do 2012 k=m4,m3
                  lad4=lad2+ladd(k+los)
                  lad5=lad3+ladd(k+los)
                  z=-int(lad4+3)-int(lad5+2)+int(lad5+1)
                  s(ii,ivc)=s(ii,ivc)+z*cjj
                  sjj=sjj+z*c(ii,ivc)
                  ii=ii+1
 2012          continue
               s(jj,ivc)=s(jj,ivc)+sjj
               jj=jj+1
 2011       continue
 2010    continue
 8010    continue
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
               ksm=xor(ss(i),ss(j))
               kos=os(ksm+1)
               lad1=ijadd(i*(i-1)/2+j)
               cjj=c(jj,ivc)
               sjj=0.0
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
                  lad2=lad1+kadd(k+kos)+ladd(m4+los)
                  do 2016 l=m4,m3
                     z=int(lad2+1)-int(lad2+3)
                     s(ii,ivc)=s(ii,ivc)+z*cjj
                     sjj=sjj+z*c(ii,ivc)
                     ii=ii+1
                     lad2=lad2+3
 2016             continue
 2015          continue
               s(jj,ivc)=s(jj,ivc)+sjj
               jj=jj+1
 2014       continue
 2013    continue
 8013    continue
c
c      no  21
c
         if(m5.lt.m6) go to 8017
         do 2017 i=m6,m5
            ssj=xor(ss(i),js)
            if (ssj.gt.ss(i)) go to 2017
            jja=ja+wtx(i,js+1)
            do 2018 j=3,i-1
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
                  lsm=xor(ksm,ss(k))
                  los=os(lsm+1)
                  lad2=lad1+kadd(k+kos)+ladd(m4+los)
                  cii=c(ii,ivc)
                  sii=0.0
                  do 2020 l=m4,m3
                     z=-int(lad2+1)+int(lad2+2)
                     sii=sii+z*c(jj,ivc)
                     s(jj,ivc)=s(jj,ivc)+z*cii
                     jj=jj+1
                     lad2=lad2+3
 2020             continue
                  s(ii,ivc)=s(ii,ivc)+sii
 2019          ii=ii+1
 2018       continue
 2017    continue
 8017    continue
c
c   no  22
c
         if(m5.lt.m6) go to 9001
            do 2021 i=m6,m5
            ssj=xor(ss(i),js)
            if (ssj.gt.ss(i)) go to 2021
            jja=ja+wtx(i,js+1)
            do 2022 j=3,i-1
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
                  lad2=lad1+kadd(k+kos)+ladd(m4+los)
                  cjj=c(jj,ivc)
                  sjj=0.0
                  do 2024 l=m4,m3
                     z=-int(lad2+3)+int(lad2+2)
                     s(ii,ivc)=s(ii,ivc)+z*cjj
                     sjj=sjj+z*c(ii,ivc)
                     ii=ii+1
                     lad2=lad2+3
 2024             continue
                  s(jj,ivc)=s(jj,ivc)+sjj
 2023          jj=jj+1
 2022       continue
 2021    continue
 9001 continue
c
      return
c
c   ww   no 69a four externals
c***********************************************************************
  316 continue
      do 9002 ivc=1,nvc
      m5=imax
      m6=imin
      if (m6.lt.2) m6=2
      if (m5.lt.m6) go to 6586
      if (is.ne.0.or.js.ne.0) go to 6586
      do 2031 i=m6,m5
         jj=ja+wab(i)
         lad=ijadd(i*(i+1)/2)
         cjj=c(jj,ivc)
         kos=os(1)
         sjj=0.0
         do 2032 j=1,i-1
            los=os(ss(j)+1)
            lad1=kadd(j+kos)+ladd(j+los)+lad
            ii=ia+wab(j)
            z=int(lad1+1)
            s(ii,ivc)=s(ii,ivc)+z*cjj
            sjj=sjj+z*c(ii,ivc)
 2032    continue
         s(jj,ivc)=s(jj,ivc)+sjj
 2031 continue
c
c    no   72a
c
 6586 if (js.ne.0) go to  6635
      if (m5.lt.m6) go to 6635
      do 2033 i=m6,m5
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
         lad=ijadd(i*(i+1)/2)+kadd(i+kos)+ladd(m2+los)
         cjj=c(jj,ivc)
         sjj=0.0
         do 2034 j=m2,m1
            z=sqrt2*(int(lad+3)+int(lad+1))
            s(ii,ivc)=s(ii,ivc)+z*cjj
            sjj=sjj+z*c(ii,ivc)
            lad=lad+3
 2034    ii=ii+1
         s(jj,ivc)=s(jj,ivc)+sjj
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
            if (ssj.gt.ss(j)) go to 2036
            if (ssj.eq.ss(j)) m3=j-1
            if (m3.lt.m4) go to 2036
            jj=jja+wty(m4)
            ii=ia+wab(j)
            ksm=xor(ss(i),ss(j))
            kos=os(ksm+1)
            los=os(ss(i)+1)
            lad1=ijadd(i*(i-1)/2+j)+kadd(j+kos)+ladd(m4+los)
            cii=c(ii,ivc)
            sii=0.0
            do 2037 k=m4,m3
               z=sqrt2*int(lad1+1)
               sii=sii+z*c(jj,ivc)
               s(jj,ivc)=s(jj,ivc)+z*cii
               lad1=lad1+3
 2037       jj=jj+1
            s(ii,ivc)=s(ii,ivc)+sii
 2036    continue
 2035 continue
c
c   no  74a
c
 6636 if (js.ne.0) go to 6637
      if (m5.lt.m6) go to 6637
      do 2038 i=m6,m5
         jj=ja+wab(i)
         lad=ijadd(i*(i+1)/2)
         cjj=c(jj,ivc)
         kos=os(1)
         sjj=0.0
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
            lad1=lad+kadd(j+kos)+ladd(m2+los)
            do 2040 k=m2,m1
               z=sqrt2*int(lad1+1)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               lad1=lad1+3
 2040       ii=ii+1
 2039    continue
         s(jj,ivc)=s(jj,ivc)+sjj
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
            ksm=xor(ss(i),ss(j))
            kos=os(ksm+1)
            lad1=ijadd(i*(i-1)/2+j)
            cjj=c(jj,ivc)
            sjj=0.0
            do 2043 k=1,j-1
               ii=ia+wab(k)
               lsm=xor(ksm,ss(k))
               los=os(lsm+1)
               lad2=lad1+kadd(k+kos)+ladd(k+los)
               z=sqrt2*int(lad2+1)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
 2043       continue
            s(jj,ivc)=s(jj,ivc)+sjj
            jj=jj+1
 2042    continue
 2041 continue
 8041 continue
c
c    no  78
c
      m6=imin
      if (m6.lt.2) m6=2
      if (m5.lt.m6) go to 6638
      do 2044 i=m6,m5
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
         lad=ijadd(i*(i+1)/2)+kadd(i+kos)+ladd(m2+los)
         do 2045 j=m2,m1
            ii=ia+wab(j)
            z=sqrt2*(int(lad+3)+int(lad+2))
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
            lad=lad+3
 2045    jj=jj+1
 2044 continue
c
c    no   84a
c
 6638 m6=imin
      if (m6.lt.3) m6=3
      if (m5.lt.m6) goto 8054
      do 2048 i=m6,m5
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
            z1=int(lad3+ladd(j+los1))
            ii=ia+wtw(j,is+1)+wty(m4)
            jj=jja+wty(m4)
            ksm=xor(ss(i),ss(j))
            kos=os(ksm+1)
            lad1=ijadd(i*(i-1)/2+j)
cdir$ ivdep
            do 2050 k=m4,m3
               lsm=xor(ksm,ss(k))
               los=os(lsm+1)
               lad2=lad1+kadd(k+kos)+ladd(k+los)
               z=int(lad2+1)+int(lad2+2)+z1
               s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
               s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
               ii=ii+1
 2050       jj=jj+1
 2049    continue
 2048 continue
c
c   no  88a
c
      do 2051 i=m6,m5
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
            lad2=ijadd(i*(i-1)/2+j)+kadd(j+kos)
            cjj=c(jj,ivc)
            sjj=0.0
            do 2053 k=m4,m3
               z=int(lad1+ladd(k+los1)+3)+int(lad2+ladd(k+los1)+2)
               z=z+int(lad2+ladd(k+los1)+1)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
 2053       ii=ii+1
            s(jj,ivc)=s(jj,ivc)+sjj
 2052    jj=jj+1
 2051 continue
c
c   no   89 & 160
c
      do 2054 i=m6,m5
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
            los=os(ss(j)+1)
            lad1=lad+kadd(j+kos)+ladd(m4+los)
            cjj=c(jj,ivc)
            sjj=0.0
            do 2056 k=m4,m3
               z=int(lad1+1)+int(lad1+2)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               lad1=lad1+3
 2056       ii=ii+1
            s(jj,ivc)=s(jj,ivc)+sjj
 2055    jj=jj+1
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
            los=os(ss(j)+1)
            lad3=ijadd(j*(j+1)/2)+kadd(j+kos)+3+ladd(m4+los)
            if (ssi.eq.ssj) m3=j-1
            if (m3.lt.m4) go to 7088
            ii=iia+wty(m4)
            cjj=c(jj,ivc)
            sjj=0.0
            do 7089 k=m4,m3
               z=int(lad3)
               s(ii,ivc)=s(ii,ivc)+z*cjj
               sjj=sjj+z*c(ii,ivc)
               lad3=lad3+3
 7089       ii=ii+1
            s(jj,ivc)=s(jj,ivc)+sjj
 7088    jj=jj+1
 7087 continue
 7080 continue
c     no    99a
c
      m5=imax
      m6=imin
      if (m6.lt.4) m6=4
      if (m5.lt.m6) go to 9002
      do 2057 i=m6,m5
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i)) go to 2057
         jja=ja+wtw(i,js+1)
         do 2058 j=3,i-1
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
               lad2=lad1+kadd(k+kos)+ladd(m4+los)
               cjj=c(jj,ivc)
               sjj=0.0
               do 2060 l=m4,m3
                  z=int(lad2+3)+int(lad2+2)
                  s(ii,ivc)=s(ii,ivc)+z*cjj
                  sjj=sjj+z*c(ii,ivc)
                  lad2=lad2+3
 2060          ii=ii+1
               s(jj,ivc)=s(jj,ivc)+sjj
 2059       jj=jj+1
 2058    continue
 2057 continue
c
c     no   100a
c
      do 2061 i=m6,m5
         ssj=xor(ss(i),js)
         if (ssj.gt.ss(i)) go to 2061
         jja=ja+wtw(i,js+1)
         do 2062 j=3,i-1
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
               lad2=lad1+kadd(k+kos)+ladd(m4+los)
               cii=c(ii,ivc)
               sii=0.0
               do 2064 l=m4,m3
                  z=int(lad2+1)+int(lad2+2)
                  sii=sii+z*c(jj,ivc)
                  s(jj,ivc)=s(jj,ivc)+z*cii
                  lad2=lad2+3
 2064          jj=jj+1
               s(ii,ivc)=s(ii,ivc)+sii
 2063       ii=ii+1
 2062    continue
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
            ksm=xor(ss(i),ss(j))
            kos=os(ksm+1)
            lad1=ijadd(i*(i-1)/2+j)
            cjj=c(jj,ivc)
            sjj=0.0
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
               lad2=lad1+kadd(k+kos)+ladd(m4+los)
               do 2068 l=m4,m3
                  z=int(lad2+1)+int(lad2+3)
                  s(ii,ivc)=s(ii,ivc)+z*cjj
                  sjj=sjj+z*c(ii,ivc)
                  lad2=lad2+3
 2068          ii=ii+1
 2067       continue
            s(jj,ivc)=s(jj,ivc)+sjj
            jj=jj+1
 2066    continue
 2065 continue
 9002 continue
      return
c
c     closed internal loop
c***********************************************************************
  318 z=val1*(int(arr+tr1)+val2*int(arr+tr2))
      ii=ia
      jj=ja
      do 4000 i=1,n1
cdir$ ivdep
         do 9003 ivc=1,nvc
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
 9003    continue
         ii=ii+1
 4000 jj=jj+1
      return
c
c     ----- internal case for tracks of (3,2,1) -----
c***********************************************************************
  319 continue
      z=val1*(int(arr+1)+int(arr+2)+int(arr+3))
      ii=ia
      jj=ja
      do 3999 i=1,n1
cdir$ ivdep
         do 9004 ivc=1,nvc
            s(jj,ivc)=s(jj,ivc)+z*c(ii,ivc)
            s(ii,ivc)=s(ii,ivc)+z*c(jj,ivc)
 9004    continue
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
