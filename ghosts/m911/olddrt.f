*deck @(#)olddrt.f	1.1  11/30/90
      subroutine olddrt(levpt,      aval,bval,rowsym,arc
     #,                 nlwks,              weight,wab,wtw,wtx
     #,                 wty,kadd,ladd,ijadd,ijgrp,inint,inext,jmnnxt
     #,                 jmxnxt,ningrp,orbsym,imngrp,imxgrp,jmngrp,jmxgrp
     #,                 ijxx,klxx,nklxx,ijww,klww,nklww)
c
c
c
      implicit integer (a-z)
      real*8 civer,drtver,rx
c
      common /lbls/ lblint(26),lbldrt(26),civer,drtver
      common /tapes/itap20,itape5,itape6,itape8,itap12,itap03,itap04
     *,             itape3,itap05,itap06
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,nrow4o,nwksoc,nlevoc,norboc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /all/  rx(3),ix(7),lvfrm1,nlwki,nlwkj,imax,imin
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /symm/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
     #,             numsym(8)
      common /xx4xbl/ nijvir
c
      dimension levpt(nlevs),aval(nrows),bval(nrows)
      dimension arc(nrows4),nlwks(nrows),rowsym(nrows)
      dimension                weight(nrows4),wab(lvfrm1)
      dimension wtw(lvfrm1,nsym),wtx(lvfrm1,nsym),wty(lvfrm1)
      dimension ijadd(numij),kadd(symorb),ladd(symorb)
      dimension ijgrp(numij),inint(norbs),inext(norbs)
      dimension jmnnxt(norbs),jmxnxt(norbs),orbsym(norbs)
      dimension ningrp(ngroup),imngrp(ngroup),imxgrp(ngroup)
      dimension jmngrp(ngroup),jmxgrp(ngroup)
      dimension ijxx(numij ),klxx(nijvir),nklxx(nsym)
      dimension ijww(numij ),klww(nijvir),nklww(nsym)
c
      read (itape8) kadd
      read (itape8) ladd
      read (itape8) ijadd
      read (itape8) ijgrp
      read (itape8) inint
      read (itape8) inext
      read (itape8) jmnnxt
      read (itape8) jmxnxt
      read (itape8) ningrp
      read (itape8)
      read (itape8) orbsym
      read (itape8)
      read (itape8)
      read (itape8)
      read (itape8) levpt
      read (itape8)
      read (itape8) aval
      read (itape8) bval
      read (itape8) rowsym
      read (itape8) arc
      read (itape8) nlwks
      if(drtver.lt.2.0d+00) read (itape8)
      if(drtver.lt.2.0d+00) read (itape8)
      if(drtver.lt.2.0d+00) read (itape8)
      read (itape8) weight
      read (itape8) wab
      read (itape8) wtw
      read (itape8) wtx
      read (itape8) wty
      if (drtver.lt.3.0d+00) go to 9
      read (itape8) ijxx
      read (itape8) klxx
      read (itape8) nklxx
      read (itape8) ijww
      read (itape8) klww
      read (itape8) nklww
    9 continue
      nfrmin=levpt(levfrm)+1
      nlwkmx=0
      do 10 i=nfrmin,nrowoc
      if(nlwks(i).gt.nlwkmx) nlwkmx=nlwks(i)
   10 continue
      iblk=1
      imxgrp(iblk)=norbs
      jmxgrp(iblk)=norbs
      do 25 i=norbs,1,-1
      do 20 j=i,1,-1
      ij=i*(i-1)/2+j
c     print '(' i,j',2i3,' iblk and ijgrp(ij)',2i3)',i,j,iblk,ijgrp(ij)
      if(ijgrp(ij).eq.iblk) go to 20
      if (i.eq.j.and.i.ne.imxgrp(iblk)) go to 22
      imngrp(iblk)=i
      jmngrp(iblk)=j+1
      go to 24
  22  continue
      imngrp(iblk)=i+1
      jmngrp(iblk)=1
  24  continue
      iblk=iblk+1
      imxgrp(iblk)=i
      jmxgrp(iblk)=j
  20  continue
  25  continue
      imngrp(ngroup)=1
      jmngrp(ngroup)=1
c     print '(' imx ',(20i3))',imxgrp
c     print '(' imn ',(20i3))',imngrp
c     print '(' jmx ',(20i3))',jmxgrp
c     print '(' jmn ',(20i3))',jmngrp
      if(iblk.eq.ngroup) go to 35
      write(itape6,30)iblk,ngroup
  30  format(1h ,'error in number of blocks,iblk=',i3,' ngroup=',i3)
      call lnkerr(' m911: error ')
  35  continue
      do 40 i=1,nsym
      minsym(i)=1000
   40 maxsym(i)=0
      minsym(1)=1
      ism=orbsym(1)
      do 45 i=2,lvfrm1
      if(orbsym(i).eq.ism) go to 45
      maxsym(ism)=i-1
      ism=orbsym(i)
      minsym(ism)=i
   45 continue
      maxsym(ism)=lvfrm1
      write(itape6,50) (minsym(j),j=1,nsym)
   50 format(' minsym=',8i5)
      write(itape6,55) (maxsym(j),j=1,nsym)
   55 format(' maxsym=',8i5)
      return
      end
