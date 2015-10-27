*deck @(#)getdrt.f	2.1  10/10/91
      subroutine getdrt(levpt,      aval,bval,rowsym,arc
     #,                 nlwks,              weight,wab,wtw,wtx
     #,                 wty,kadd,ladd,ijadd,ijgrp,inint,inext,jmnnxt
     #,                 jmxnxt,ningrp,orbsym,imngrp,imxgrp,jmngrp,jmxgrp
     #,                 ijxx,klxx,nklxx,ijww,klww,nklww)
c
c
c
      implicit integer (a-z)
c
cvax  extended dummy    levpt,      aval,bval,rowsym,arc
cvax  extended dummy    nlwks,              weight,wab,wtw,wtx
cvax  extended dummy    wty,kadd,ladd,ijadd,ijgrp,inint,inext,jmnnxt
cvax  extended dummy    jmxnxt,ningrp,orbsym,imngrp,imxgrp,jmngrp,jmxgrp
cvax  extended dummy    ijxx,klxx,nklxx,ijww,klww,nklww
c
      real*8 civer,drtver,q
c
      common /lbls/ lblint(26),lbldrt(26),civer,drtver
      common /io/     itape5,itape6
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /all/  q(3),iq,ix(6),lvfrm1,nlwki,nlwkj,imax,imin
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /symq/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
     #,             numsym(8)
      common /x4x901/ nijvir,nrefs
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
      dimension ijxx(numij ),klxx(nijvir),nklxx(nsym,orbfrm)
      dimension ijww(numij ),klww(nijvir),nklww(nsym,orbfrm)
c
      call iosys('read integer kadd from rwf',-1,kadd,0,' ')
      call iosys('read integer ladd from rwf',-1,ladd,0,' ')
      call iosys('read integer ijadd from rwf',-1,ijadd,0,' ')
      call iosys('read integer ijgrp from rwf',-1,ijgrp,0,' ')
      call iosys('read integer inint from rwf',-1,inint,0,' ')
      call iosys('read integer inext from rwf',-1,inext,0,' ')
      call iosys('read integer jmnnxt from rwf',-1,jmnnxt,0,' ')
      call iosys('read integer jmxnxt from rwf',-1,jmxnxt,0,' ')
      call iosys('read integer ningrp from rwf',-1,ningrp,0,' ')
      call iosys('read integer orbsym from rwf',-1,orbsym,0,' ')
      call iosys('read integer levpt from rwf',-1,levpt,0,' ')
      call iosys('read integer a from rwf',-1,aval,0,' ')
      call iosys('read integer b from rwf',-1,bval,0,' ')
      call iosys('read integer s from rwf',-1,rowsym,0,' ')
      call iosys('read integer arc from rwf',-1,arc,0,' ')
      call iosys('read integer nlwks from rwf',-1,nlwks,0,' ')
      call iosys('read integer weight from rwf',-1,weight,0,' ')
      call iosys('read integer wtab from rwf',-1,wab,0,' ')
      call iosys('read integer wtw from rwf',-1,wtw,0,' ')
      call iosys('read integer wtx from rwf',-1,wtx,0,' ')
      call iosys('read integer wty from rwf',-1,wty,0,' ')
      call iosys('read integer ijxx from rwf',-1,ijxx,0,' ')
      call iosys('read integer klxx from rwf',-1,klxx,0,' ')
      call iosys('read integer nklxx from rwf',-1,nklxx,0,' ')
      call iosys('read integer ijww from rwf',-1,ijww,0,' ')
      call iosys('read integer klww from rwf',-1,klww,0,' ')
      call iosys('read integer nklww from rwf',-1,nklww,0,' ')
c
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
  22        continue
            imngrp(iblk)=i+1
            jmngrp(iblk)=1
  24        continue
            iblk=iblk+1
            imxgrp(iblk)=i
            jmxgrp(iblk)=j
  20     continue
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
      stop
  35  continue
      do 40 i=1,nsym
         minsym(i)=1000
   40 maxsym(i)=0
      ism=orbsym(1)
      minsym(ism)=1
      do 45 i=2,lvfrm1
         if(orbsym(i).eq.ism) go to 45
         maxsym(ism)=i-1
         ism=orbsym(i)
         minsym(ism)=i
   45 continue
      maxsym(ism)=lvfrm1
      return
      end
