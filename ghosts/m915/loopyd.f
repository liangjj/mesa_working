*deck @(#)loopyd.f	2.1  10/10/91
      subroutine loopyd(nabc,nlwks,           iwght,nabca,ijadd
     *,                 ijgrp,kadd,ladd,inext,iarc,jmnnxt,jmxnxt
     *,                 nabcs,isym,            isegm,jsegm,imain
     *,                 isub,iuwkmn,iuwksb,itrack,lmin,acoef,bcoef
     *,                 hdwgt,hdseg,hdpnt,imxgrp,imngrp,jmxgrp,jmngrp
     *,                 aint,file,unit)
c
c
c
      implicit real*8 (a-h,o-z)
c
      integer xor
      character*16 file,unit
c
      character*16 path
      dimension aint(nmax)
      dimension nlcsmn(22),lcond(8),coeffs(20,21),cfs(420)
      integer puwkt,refwlk,symorb,bmax,orbfrm
      integer hpt,hseg,hsegmx,hdwgt,hdseg,hdpnt,segwgt,harpt
      logical pagein
c   level characteristics
      dimension isegm(nlevs),jsegm(nlevs),imain(nlevs),isub(nlevs)
      dimension iuwkmn(nlevs),iuwksb(nlevs),itrack(nlevs),isym(norbs)
      dimension acoef(nlevs),bcoef(nlevs),lmin(nlevs)
      dimension hdseg(nlevs),hdpnt(nlevs),hdwgt(nlevs)
c     dimension levpt(nlevs),levnr(nlevs)
c   graph description arrays
      dimension nabc(nrows ),nlwks(nrows)              ,nabcs(nrows )
      dimension              iwght(nrows4),iarc(nrows4),nabca(nrows )
c   integral addressing arrays
      dimension inext(norbs),jmnnxt(norbs),jmxnxt(norbs)
      dimension ijadd(numij),ijgrp(numij),kadd(symorb),ladd(symorb)
      dimension imngrp(ngroup),imxgrp(ngroup),jmngrp(ngroup)
      dimension jmxgrp(ngroup)
c
      dimension ibuf(30,1000),xbuf(30,1000)
      equivalence (ibuf(1,1),xbuf(1,1))
c
      common /level/ levi
c
      common /symq/ jsmt,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
     #,             numsym(8)
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /d901/ rep,fzcore,eguess,eci,refwlk,mxiter,cnverg,icnvg
     *,             iter,sqcdif,czero,nroot
      common /io/     itape5,itape6
      common /optns/  iguess,irstrt,irooti,irootf,ikl34
      common/loops/ nuwk,puwkt,iuwk,juwk,itrak,ipt1,ipt2
      common/all/   acf,d,ccf,ladt,itr1,itr2,ia,ja,itype,isegt,lvfrm1
     *,             nlwki,nlwkj,imax,imin
      common /page/ iword3,iword4,pagein,noffi,noffj,npass,nwordl,iword6
      common /rstrt/iblck1,inxt,lowtri,ndvdit,iblock
      common /tables/ jsegnr(22),jsegpt(22),iarcmn(228),iarcsb(228)
     *,itrk(228),jcond(228),kcond(228),nxtseg(228),jsegpx(3)
      common /minmax/ iming,imaxg,jmin,jmax
c
      common /nbhl/ ibhl
c
      equivalence (coeffs,cfs)
      data nlcsmn/6*1,0,1,14*0/
c
      crite = 0.00001d+00
      root2 = sqrt(2.0d+00)
      rootn2 = -root2
      toor2 = 1.0d+00 / root2
      toorn2 = -toor2
      jsegpt(1)=0
      hsegmx=4
      do 125 i=1,norbs
  125 isym(i)=isym(i)-1
      do 130 j=1,21
         jsegpt(j+1)=jsegnr(j)
         do 130 i=1,2
 130  coeffs(i,j)=0.0
      do 140 i=3,20
         a = dble(i-2)
         coeffs(i,1) = sqrt(a/(a+1.0d+00))
         coeffs(i,2) = -coeffs(i,1)
         coeffs(i,3) = coeffs(i,1)/sqrt(2.0d+00)
         coeffs(i,4) = -coeffs(i,3)
         coeffs(i,5) = sqrt((a+1.0d+00)/a)
         coeffs(i,6) = -coeffs(i,5)
         coeffs(i,7) = coeffs(i,5)/sqrt(2.0d+00)
         coeffs(i,8) = -coeffs(i,7)
         coeffs(i,9) = sqrt((a+2.0d+00)/(a*2.0d+00))
         coeffs(i,10) = -coeffs(i,9)
         coeffs(i,11) = sqrt(a/(2.0d+00*(a+2.0d+00)))
         coeffs(i,12) = -coeffs(i,11)
         coeffs(i,13) = sqrt(2.0d+00/(a*(a+1.0d+00)))
         coeffs(i,14) = sqrt(a*(a+2.0d+00))/(a+1.0d+00)
         coeffs(i,15) = -sqrt(a*(a+2.0d+00))/(a+1.0d+00)
         coeffs(i,16) = sqrt((a-1.0d+00)*(a+2.0d+00)/(a*(a+1.0d+00)))
         coeffs(i,17) = -coeffs(i,16)
         coeffs(i,18)=-sqrt(2.0d+00/(a*(a+2.0d+00)))
         coeffs(i,19) = 1.0d+00/a
         coeffs(i,20) = -1.0d+00/a
         coeffs(i,21) = -sqrt(2.0d+00)/a
 140  continue
      do 160 i=1,nsym
         ismoff(i)=(i-1)*norbs
 160  lcond(i)=0
      i=isym(1)+1
      lcond(i)=1
      lcond(1)=1
      nsm=0
      do 170 iorb=2,norbs
         do 165 i=1,nsym
            if(lcond(i).eq.0) go to 165
            ism=i-1
            j=xor(ism,isym(iorb))+1
            if(lcond(j).gt.0) go to 165
            lcond(j)=iorb
            nsm=nsm+1
            if(nsm.eq.nsym) go to 175
 165     continue
 170  continue
 175  continue
      do 180 i=1,nsym
         if(lcond(i).eq.0) lcond(i)=norbs+1
 180  continue
c..bhl
      call iosys('read integer ndci from rwf',1,ndci,0,' ')
      ldci=30*maxbuf*ndci
      call iosys('read integer nfci from rwf',1,nfci,0,' ')
      lfci=30*maxbuf*nfci
      ntci=ldci+lfci+maxbuf
c..bhl
      ibhl=0
c
c--------------------------------------------------------------------loopy
c
      entry loopy
c
      call iosys('rewind "'//file//'" on '//unit,0,0,0,' ')
c
      ibhl=ibhl+1
      if(ibhl.eq.1) then
       path='diagonal form'
      call iosys('read integer ndci from rwf',1,ndci,0,' ')
      ldci=30*maxbuf*ndci
      call iosys('create real "'//path//'" on ciform',
     $ ldci,0,0,' ')
      else
       path='ci form'
      call iosys('read integer nfci from rwf',1,nfci,0,' ')
      lfci=30*maxbuf*nfci
      call iosys('create real "'//path//'" on ciform',
     $ lfci,0,0,' ')
      end if
c
      noutf=0
      ixf=1
      maxbuf=1000
c
c      call iosys('rewind "'//file//'" on '//unit,0,0,0,' ')
c
      do 500 iblock=1,ngroup
c         call iosys('read real "'//file//'" from '//unit//
c     #              ' without rewinding',nmax,aint,0,' ')
c         if(ndvdit.gt.0) goto 188
c         write(itape6,187)iblock,imxgrp(iblock),imngrp(iblock),
c     #   jmxgrp(iblock),jmngrp(iblock)
c 187     format(' process integrals from group',i4,' first level',i4
c     #   ,      ' last',i4,' jmax',i4,' jmin',i4)
c 188     continue
         imaxg=imxgrp(iblock)
         iming=imngrp(iblock)
         i=imaxg
c
 190     continue
         levi=i+1
         npt=1
         puwkt=1
         jmax=i
         jmin=1
         if (i.eq.imaxg) jmax=jmxgrp(iblock)
         if (i.eq.iming) jmin=jmngrp(iblock)
         jsm=isym(i)
         iad=(i*(i-1))/2
         ij=iad+jmax
         noffi=0
         noffj=0
         levv=nlevs
         hseg=4
         levh=nlevs
         pagein=.false.
         if(nwks.le.nwksmx) pagein=.true.
         if(levi.eq.nlevs) go to 1070
c
c     ----- generate head segments to loops -----
c
         levh=nlevs+1
         npt=1
         hdwgt(nlevs)=0
         hdpnt(nlevs)=1
 1000    levh=levh-1
         segwgt=0
         hseg=0
         hpt=npt
 1010    hseg=hseg+1
         if(hseg.le.hsegmx) go to 1030
         levh=levh+1
         if(levh.gt.nlevs) go to 480
         puwkt=puwkt-hdwgt(levh)
         hseg=hdseg(levh)
         hpt=hdpnt(levh)
         go to 1010
 1030    continue
         harpt=(hpt-1)*4+hseg
         npt=iarc(harpt)
         if(npt.eq.0) go to 1010
         puwkt=puwkt-segwgt
         segwgt=iwght(harpt)
         puwkt=puwkt+segwgt
         hdwgt(levh)=segwgt
         hdseg(levh)=hseg
         lev=levh-1
         hdpnt(lev)=npt
         iuwk=0
         juwk=0
         nlwki=nlwks(npt)
         if(pagein) go to 1035
         call paged
         if(pagein) levv=lev
         go to 1040
 1035    continue
         if(levv.gt.lev) go to 1040
         call pageout
         call paged
         levv=nlevs
         if(pagein) levv=lev
 1040    continue
         if(lev.gt.levfrm) go to 1050
         if((nabca(npt)*2+nabc(npt)).le.2) go to 1060
 1050    if(lev.gt.levi) go to 1000
         go to 1070
 1060    continue
c
c     ----- loops with four external indices -----
c
c 10/20/89
c..bhl         if(ndvdit.eq.1.and.iguess.eq.0) go to 480
c 10/20/89
c
         itype=3*nabca(npt)+nabc(npt)
         if(itype.eq.0) go to 1010
         ifsym=nabcs(npt)
         jfsym=ifsym
         lvfrm1=lev-1
         imax=levi-1
         imin=levi-1
         if(.not.pagein) stop 'paging malfunction'
c..bhl
         ixf=ixf+1
         ibuf(1,ixf)=ladt
         xbuf(2,ixf)=acf
         xbuf(3,ixf)=d
         xbuf(4,ixf)=ccf
         ibuf(5,ixf)=itr1
         ibuf(6,ixf)=itr2
         ibuf(7,ixf)=ia
         ibuf(8,ixf)=ja
         ibuf(9,ixf)=itype
         ibuf(10,ixf)=isegt
         ibuf(11,ixf)=lvfrm1
         ibuf(12,ixf)=nlwki
         ibuf(13,ixf)=nlwkj
         ibuf(14,ixf)=imax
         ibuf(15,ixf)=imin
         ibuf(16,ixf)=nuwk
         ibuf(17,ixf)=puwkt
         ibuf(18,ixf)=iuwk
         ibuf(19,ixf)=juwk
         ibuf(20,ixf)=itrak
         ibuf(21,ixf)=ipt1
         ibuf(22,ixf)=ipt2
         ibuf(23,ixf)=iming
         ibuf(24,ixf)=imaxg
         ibuf(25,ixf)=jmin
         ibuf(26,ixf)=jmax
         ibuf(27,ixf)=jsmt
         ibuf(28,ixf)=ifsym
         ibuf(29,ixf)=jfsym
         ibuf(30,ixf)=3
c
      if(ixf.eq.maxbuf) then
       ibuf(1,1)=ixf
       ibuf(2,1)=iblock
       ibuf(3,1)=0
       call iosys('write real "'//path//'" to ciform without rewinding',
     $ 30*maxbuf,xbuf,0,' ')
       ixf=1
       noutf=noutf+1
      end if
c..bhl
c         call allext
c..bhl
         go to 1010
 1070    lev=levi
         levm=lev-1
         isegm(lev)=1
         iseg=1
         imn=npt
         isb=npt
         kseg=0
         ksegmx=jsegnr(iseg)
         lmin(lev)=lcond(jsm+1)
         iuwkmn(lev)=0
         iuwksb(lev)=0
         imain(lev)=npt
         isub(lev)=npt
         acoef(lev)=1.0d+00
c
c     ----- test next segment of group -----
c
200      kseg=kseg+1
         if(kseg.gt.ksegmx) go to 440
         ksb=iarcsb(kseg)
         jarpt=4*(isb-1)+ksb
         ksb=iarc(jarpt)
         if(ksb.eq.0) go to 200
         kmn=iarcmn(kseg)
         iarpt=4*(imn-1)+kmn
         kmn=iarc(iarpt)
         if(kmn.eq.0) go to 200
         if(ndvdit.gt.0) go to 210
         if(kmn.ne.ksb) go to 200
 210     continue
         jsegm(lev)=kseg
         iuwkmn(levm)=iuwkmn(lev)+iwght(iarpt)
         iuwksb(levm)=iuwksb(lev)+iwght(jarpt)
         lmin(levm)=lmin(lev)
         if(jcond(kseg))220,240,230
 220     continue
         if (levm.le.jmin) go to 200
         goto 240
 230     continue
         if (levm.gt.jmax) go to 200
         ij=iad+levm
         jad=ijadd(ij)
         ksm=xor(jsm,isym(levm))
         lmin(levm)=lcond(ksm+1)
         ksmptx=ismoff(ksm+1)
 240     continue
         if(kcond(kseg).eq.0)goto 260
         ksmpt=levm+ksmptx
         kad=jad+kadd(ksmpt)
         lsm=xor(ksm,isym(levm))
         lmin(levm)=lcond(lsm+1)
         lsmptx=ismoff(lsm+1)
 260     continue
         if(itrk(kseg))280,280,270
270      itrack(levm)=itrk(kseg)
         goto 290
280      itrack(levm)=itrack(lev)
290      continue
         go to( 1, 1, 1, 1, 1, 1,44, 1, 3, 3, 1,45, 4, 4,50,51, 1,40, 1,
     *    1,     6, 6, 7,46, 9,54, 2,52,41, 5,47, 8, 2,53, 1, 1,42,48, 2
     *   ,52,    36,55,43, 1,11,11,12,10,13,49, 2,53, 1,77,77,77, 1,79,
     *   77, 1,    80,78, 1,71,67,68,87,75,83,69,68,76,70,82,71,71,67,68
     *   ,67,87,    75,83,69,68,83,68,76,70,69,70,82,71, 1, 6,16, 6, 6,
     *   17,16,74,     8, 1, 1,18,19,18,18,22,24,20,19,24,19,23,21,20,21
     *   , 1, 1,11,    11,27,11,28,81,27,13, 1, 1, 3, 2, 4, 2, 1, 2, 2,
     *   1,71,63,72,    84,65,85,73,29,66,64,71, 1,30,56,86, 1,57, 1, 1,
     *   56, 1,32,31,    58, 1, 1,59,33,34,61,35,88,62,60, 1, 1, 6, 1, 9
     *   , 2, 5, 2, 1,     1, 2,36,11,10, 2, 1, 6, 1, 9, 2, 5, 2, 1, 1,
     *   2,36,11,10, 2,     1, 1,37,38, 4, 2, 2, 1, 1, 2,39,37, 3, 2, 1,
     *    6, 1, 9, 2, 5,     2, 1, 1, 2,36,11,10, 2) kseg
  1      acoef(levm)=acoef(lev)
         goto 120
  2      acoef(levm)=-acoef(lev)
         goto 120
   3     ia = nabc(imn) + 2
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
   4     ia = nabc(imn) + 83
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
   5     ia = nabc(imn) + 82
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
   6     ia = nabc(imn) + 261
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
   7     ia = nabc(imn) + 1
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
   8     ia = nabc(imn) + 102
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
   9     ia = nabc(imn) + 362
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  10     ia = nabc(imn) + 3
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  11     ia = nabc(imn) + 263
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  12     ia = nabc(imn) + 84
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  13     ia = nabc(imn) + 23
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  16     ia = nabc(imn) + 281
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  17     ia = nabc(imn) + 402
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  18     ia = nabc(imn) + 162
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  19     ia = nabc(imn) + 222
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  20     ia = nabc(imn) + 143
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  21     ia = nabc(imn) + 42
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  22     ia = nabc(imn) + 302
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  23     ia = nabc(imn) + 303
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  24     ia = nabc(imn) + 342
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  27     ia = nabc(imn) + 283
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  28     ia=nabc(imn)+404
         acoef(levm)=acoef(lev)*cfs(ia)
         go to 120
  29     acoef(levm) = acoef(lev) * root2
         go to 120
  30     ia = nabc(imn) + 301
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  31     ia = nabc(imn) + 304
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  32     ia = nabc(imn) + 244
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  33     ia = nabc(imn) + 322
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  34     ia = nabc(imn) + 243
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  35     ia = nabc(imn) + 242
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  36     ia = nabc(imn) + 384
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  37     ia = nabc(imn) + 262
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  38     ia = nabc(imn) + 363
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  39     ia = nabc(imn) + 383
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  86     ia = nabc(imn) + 241
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  40     ia = nabc(imn) + 122
         ib=ia-61
         acoef(levm) = acoef(lev) * cfs(ia)
         bcoef(levm) = acoef(lev) * cfs(ib)
         go to 120
  41     ib = nabc(imn) + 162
         acoef(levm) = acoef(lev) * toorn2
         bcoef(levm) = acoef(lev) * cfs(ib)
         go to 120
  42     ia = nabc(imn) + 43
         ib = ia + 81
         acoef(levm) = acoef(lev) * cfs(ia)
         bcoef(levm) = acoef(lev) * cfs(ib)
         go to 120
  43     ib = nabc(imn) + 222
         acoef(levm) = acoef(lev) * toorn2
         bcoef(levm) = acoef(lev) * cfs(ib)
         go to 120
  44     ib=nabc(imn)+221
         acoef(levm) = acoef(lev) * toor2
         bcoef(levm) = acoef(lev) * cfs(ib)
         go to 120
  45     ib = nabc(imn) + 163
         acoef(levm) = acoef(lev) * toor2
         bcoef(levm) = acoef(lev) * cfs(ib)
         go to 120
  46     ib = nabc(imn) + 162
         acoef(levm) = acoef(lev) * toor2
         bcoef(levm) = acoef(lev) * cfs(ib)
         go to 120
  47     ia = nabc(imn) + 122
         ib = ia - 81
         acoef(levm) = acoef(lev) * cfs(ia)
         bcoef(levm) = acoef(lev) * cfs(ib)
         go to 120
  48     ib = nabc(imn) + 222
         acoef(levm) = acoef(lev) * toor2
         bcoef(levm) = acoef(lev) * cfs(ib)
         go to 120
  49     ia = nabc(imn) + 43
         ib = ia + 101
         acoef(levm) = acoef(lev) * cfs(ia)
         bcoef(levm) = acoef(lev) * cfs(ib)
         go to 120
  50     acoef(levm) = acoef(lev) + acoef(lev)
         d=0.5d+00
         go to 120
  51     acoef(levm)=acoef(lev)*root2
         go to 120
  52     acoef(levm) = -acoef(lev)
         d= -1.0d+00
         go to 120
  53     acoef(levm) = -acoef(lev) - acoef(lev)
         d = -0.5d+00
         go to 120
  54     ia=nabc(imn)+362
         d=1.0d+00/cfs(ia)
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  55     ia = nabc(imn) + 384
         d=1.0d+00/cfs(ia)
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
  56     acoef(levm) = acoef(lev)
         d = -1.0d+00
         go to 120
  57     ia = nabc(imn) + 82
         acoef(levm) = acoef(lev) * cfs(ia)
         d=-1.0d+00
         go to 120
  58     ia = nabc(imn) + 3
         acoef(levm) = acoef(lev) * cfs(ia)
         d=-1.0d+00
         go to 120
  59     ia = nabc(imn) + 123
         acoef(levm) = acoef(lev) * cfs(ia)
         d=-1.0d+00
         go to 120
  60     ia = nabc(imn) + 222
         acoef(levm) = acoef(lev) * cfs(ia)
         d=-1.0d+00
         go to 120
  61     ia = nabc(imn) + 62
         acoef(levm) = acoef(lev) * cfs(ia)
         d=-1.0d+00
         go to 120
  62     ia = nabc(imn) + 162
         acoef(levm) = acoef(lev) * cfs(ia)
         d=-1.0d+00
         go to 120
  63     ia = nabc(imn) + 42
         ib = ia + 81
         acof = acoef(lev) * cfs(ia)
         bcof = bcoef(lev) * cfs(ib)
         d = acof + bcof
         if(abs(d).lt.crite) go to 110
         acoef(levm) = d
         d = (acof - bcof) / d
         go to 120
  64     ib = nabc(imn) + 222
         acof = acoef(lev) * toorn2
         bcof = bcoef(lev) * cfs(ib)
         d = acof + bcof
         if(abs(d).lt.crite) go to 110
         acoef(levm) = d
         d = (acof - bcof) / d
         go to 120
  65     ia = nabc(imn) + 123
         ib = ia - 61
         acof = acoef(lev) * cfs(ia)
         bcof = bcoef(lev) * cfs(ib)
         d = acof + bcof
         if(abs(d).lt.crite) go to 110
         acoef(levm) = d
         d = (acof - bcof) / d
         go to 120
  66     ib = nabc(imn) + 162
         acof = acoef(lev) * toorn2
         bcof = bcoef(lev) * cfs(ib)
         d = acof + bcof
         if(abs(d).lt.crite) go to 110
         acoef(levm) = d
         d = (acof - bcof) / d
         go to 120
  67     ib = nabc(imn) + 162
         dx=acoef(lev)*toorn2
         d=dx+bcoef(lev)*cfs(ib)
         if(abs(d).lt.crite) go to 111
         acoef(levm) = d
         d=-(dx+dx)/d
         go to 120
  68     ib = nabc(imn) + 222
         dx=acoef(lev)*toorn2
         d=dx+bcoef(lev)*cfs(ib)
         if(abs(d).lt.crite) go to 111
         acoef(levm) = d
         d=-(dx+dx)/d
         go to 120
  69     ia = nabc(imn) + 62
         ib = ia + 81
         dx=acoef(lev)*cfs(ia)
         d=dx+bcoef(lev)*cfs(ib)
         if(abs(d).lt.crite) go to 111
         acoef(levm) = d
         d=-(dx+dx)/d
         go to 120
  70     ia = nabc(imn) + 143
         ib = ia - 101
         dx=acoef(lev)*cfs(ia)
         d=dx+bcoef(lev)*cfs(ib)
         if(abs(d).lt.crite) go to 111
         acoef(levm) = d
         d=-(dx+dx)/d
         go to 120
  87     ib = nabc(imn) + 162
         dx=acoef(lev)*toorn2
         d=dx+bcoef(lev)*cfs(ib)
         if(abs(d).lt.crite) go to 111
         acoef(levm)=d
         d=-(dx+dx)/d
         go to 120
  71     acoef(levm) = acoef(lev)
         bcoef(levm) = bcoef(lev)
         go to 120
  72     ib = nabc(imn) + 322
         acoef(levm) = -acoef(lev)
         bcoef(levm) = bcoef(lev) * cfs(ib)
         go to 120
  73     ib = nabc(imn) + 323
         acoef(levm) = -acoef(lev)
         bcoef(levm) = bcoef(lev) * cfs(ib)
         go to 120
  74     ia=nabc(imn)+21
         acoef(levm)=acoef(lev)*cfs(ia)
         go to 120
  75     ib = nabc(imn) + 302
         acoef(levm) = acoef(lev)
         bcoef(levm) = bcoef(lev) * cfs(ib)
         go to 120
  76     ib = nabc(imn) + 303
         acoef(levm) = acoef(lev)
         bcoef(levm) = bcoef(lev) * cfs(ib)
         go to 120
  77     acoef(levm)=acoef(lev)*toorn2
         d=-2.0d+00
         go to 120
  78     acoef(levm)=acoef(lev)*rootn2
         d=-2.0d+00
         go to 120
  79     ia=nabc(imn)+62
         acoef(levm)=acoef(lev)*cfs(ia)
         d=-2.0d+00
         go to 120
  80     ia=nabc(imn)+143
         acoef(levm)=acoef(lev)*cfs(ia)
         d=-2.0d+00
         go to 120
  81     ia=nabc(imn)+104
         acoef(levm)=acoef(lev)*cfs(ia)
         go to 120
  82     acoef(levm) = acoef(lev) * rootn2
         d=-2.0d+00
         go to 120
  83     ia = nabc(imn) + 342
         acoef(levm) = bcoef(lev) * cfs(ia)
         go to 120
  84     ia = nabc(imn) + 243
         acoef(levm) = bcoef(lev) * cfs(ia)
         go to 120
  85     ia = nabc(imn) + 242
         acoef(levm) = bcoef(lev) * cfs(ia)
         go to 120
  88     ia = nabc(imn) + 323
         acoef(levm) = acoef(lev) * cfs(ia)
         go to 120
 110     itrack(levm)=3
         acoef(levm)=acof-bcof
         go to 120
111      itrack(levm) = 2
         acoef(levm)=-(dx+dx)
  120    continue
c
c     ----- check to see if page of vector is correct -----
c
 1100    if(pagein) go to 1110
         iuwk=iuwksb(levm)
         juwk=iuwkmn(levm)
         nlwki=nlwks(ksb)
         nlwkj=nlwks(kmn)
         call pages
         if(pagein) levv=levm
         go to 1130
 1110    continue
         if(levm.lt.levv) go to 1130
         iuwk=iuwksb(levm)
         juwk=iuwkmn(levm)
         nlwki=nlwks(ksb)
         nlwkj=nlwks(kmn)
         if(levm.gt.levv) go to 1120
         call pageij
         levv=nlevs
         if(pagein) levv=levm
         go to 1130
 1120    continue
         call pageout
         go to 1100
 1130    continue
         if(nxtseg(kseg).gt.0) go to 400
         if(isym(levm).ne.lsm) go to 200
         lsmpt=levm+lsmptx
         lad=kad+ladd(lsmpt)
         if(kmn-ksb.eq.0) go to 380
         levl=levm
         ksegmx=4
310      lev=levm
         levm=lev-1
         if(levm.gt.0) go to 315
         write(itape6,313)
 313     format(' problems with partial space')
         stop
 315     continue
         kseg=0
         imain(lev)=kmn
         imn=kmn
         isub(lev)=ksb
         isb=ksb
320      kseg=kseg+1
         if(kseg.gt.ksegmx)goto 360
         jarpt=4*(isb-1)+kseg
         ksb=iarc(jarpt)
         if(ksb.le.0)goto 320
         iarpt=4*(imn-1)+kseg
         kmn=iarc(iarpt)
         if(kmn.le.0)goto 320
         jsegm(lev)=kseg
         iuwkmn(levm)=iuwkmn(lev)+iwght(iarpt)
         iuwksb(levm)=iuwksb(lev)+iwght(jarpt)
         if(levm.gt.levfrm) goto 310
         if(kmn-ksb.ne.0) stop
         nlwki = nlwks(kmn)
         iuwk = iuwksb(levm)
         juwk = iuwkmn(levm)
         itrak = itrack(levl)
         acf = acoef(levl)
         ladt=lad
c..bhl
         ixf=ixf+1
         ibuf(1,ixf)=ladt
         xbuf(2,ixf)=acf
         xbuf(3,ixf)=d
         xbuf(4,ixf)=ccf
         ibuf(5,ixf)=itr1
         ibuf(6,ixf)=itr2
         ibuf(7,ixf)=ia
         ibuf(8,ixf)=ja
         ibuf(9,ixf)=itype
         ibuf(10,ixf)=isegt
         ibuf(11,ixf)=lvfrm1
         ibuf(12,ixf)=nlwki
         ibuf(13,ixf)=nlwkj
         ibuf(14,ixf)=imax
         ibuf(15,ixf)=imin
         ibuf(16,ixf)=nuwk
         ibuf(17,ixf)=puwkt
         ibuf(18,ixf)=iuwk
         ibuf(19,ixf)=juwk
         ibuf(20,ixf)=itrak
         ibuf(21,ixf)=ipt1
         ibuf(22,ixf)=ipt2
         ibuf(23,ixf)=iming
         ibuf(24,ixf)=imaxg
         ibuf(25,ixf)=jmin
         ibuf(26,ixf)=jmax
         ibuf(27,ixf)=jsmt
         ibuf(28,ixf)=ifsym
         ibuf(29,ixf)=jfsym
         ibuf(30,ixf)=1
c
      if(ixf.eq.maxbuf) then
       ibuf(1,1)=ixf
       ibuf(2,1)=iblock
       ibuf(3,1)=0
       call iosys('write real "'//path//'" to ciform without rewinding',
     $ 30*maxbuf,xbuf,0,' ')
       ixf=1
       noutf=noutf+1
      end if
c..bhl
c         call loopin
c..bhl
         go to 320
  360    continue
         if(lev.eq.levl) go to 440
         levm=lev
         lev=levm+1
         imn=imain(lev)
         isb=isub(lev)
         kseg=jsegm(lev)
         go to 320
 380     nlwki = nlwks(kmn)
         iuwk = iuwksb(levm)
         juwk = iuwkmn(levm)
         itrak = itrack(levm)
         acf = acoef(levm)
         ladt = lad
c..bhl
         ixf=ixf+1
         ibuf(1,ixf)=ladt
         xbuf(2,ixf)=acf
         xbuf(3,ixf)=d
         xbuf(4,ixf)=ccf
         ibuf(5,ixf)=itr1
         ibuf(6,ixf)=itr2
         ibuf(7,ixf)=ia
         ibuf(8,ixf)=ja
         ibuf(9,ixf)=itype
         ibuf(10,ixf)=isegt
         ibuf(11,ixf)=lvfrm1
         ibuf(12,ixf)=nlwki
         ibuf(13,ixf)=nlwkj
         ibuf(14,ixf)=imax
         ibuf(15,ixf)=imin
         ibuf(16,ixf)=nuwk
         ibuf(17,ixf)=puwkt
         ibuf(18,ixf)=iuwk
         ibuf(19,ixf)=juwk
         ibuf(20,ixf)=itrak
         ibuf(21,ixf)=ipt1
         ibuf(22,ixf)=ipt2
         ibuf(23,ixf)=iming
         ibuf(24,ixf)=imaxg
         ibuf(25,ixf)=jmin
         ibuf(26,ixf)=jmax
         ibuf(27,ixf)=jsmt
         ibuf(28,ixf)=ifsym
         ibuf(29,ixf)=jfsym
         ibuf(30,ixf)=1
c
      if(ixf.eq.maxbuf) then
       ibuf(1,1)=ixf
       ibuf(2,1)=iblock
       ibuf(3,1)=0
       call iosys('write real "'//path//'" to ciform without rewinding',
     $ 30*maxbuf,xbuf,0,' ')
       ixf=1
       noutf=noutf+1
      end if
c..bhl
c         call loopin
c..bhl
         go to 200
  400    continue
         if(levm.le.lmin(levm))goto 200
         iseg=nxtseg(kseg)
         if(nlcsmn(iseg).gt.0.and.(nabca(kmn)+nabc(kmn)).eq.0)goto 200
         if(levm.le.levfrm)goto 460
  410    lev=levm
         levm=lev-1
         isegm(lev)=iseg
         kseg=jsegpt(iseg)
         imn=kmn
         imain(lev)=kmn
         isb=ksb
         isub(lev)=ksb
         ksegmx=jsegnr(iseg)
         go to 200
440      continue
         if(lev.eq.levi) goto 1010
         levm=lev
         lev=levm+1
         iseg=isegm(lev)
         imn=imain(lev)
         isb=isub(lev)
         kseg=jsegm(lev)
         ksegmx=jsegnr(iseg)
         goto 200
460      continue
c
c     ----- finish loops with 3 or fewer external indices -----
c
         iamn=nabca(kmn)
         iasb=nabca(ksb)
         ibmn=nabc(kmn)
         ibsb=nabc(ksb)
         ipt1=iamn+iamn+ibmn
         ipt2=iasb+iasb+ibsb
         if(ipt1.gt.2.or.ipt2.gt.2) go to 410
         ipt1=4-ipt1-iamn
         ipt2=4-ipt2-iasb
         lvfrm1=levm-1
         iuwk=iuwksb(levm)
         juwk=iuwkmn(levm)
         itrak=itrack(levm)
         jfsym=nabcs(kmn)
         ifsym=nabcs(ksb)
         acf=acoef(levm)
         ccf=bcoef(levm)
         isegt=iseg
         if(iseg.gt.14) go to 470
         if(iseg.gt.3)  go to 465
         jsmt=jsm
         ladt=iad
         goto 475
 465     continue
         jsmt=ksm
         ladt=jad
         go to 475
 470     continue
         ladt=kad
         jsmt=lsm
 475     continue
c..bhl
         ixf=ixf+1
         ibuf(1,ixf)=ladt
         xbuf(2,ixf)=acf
         xbuf(3,ixf)=d
         xbuf(4,ixf)=ccf
         ibuf(5,ixf)=itr1
         ibuf(6,ixf)=itr2
         ibuf(7,ixf)=ia
         ibuf(8,ixf)=ja
         ibuf(9,ixf)=itype
         ibuf(10,ixf)=isegt
         ibuf(11,ixf)=lvfrm1
         ibuf(12,ixf)=nlwki
         ibuf(13,ixf)=nlwkj
         ibuf(14,ixf)=imax
         ibuf(15,ixf)=imin
         ibuf(16,ixf)=nuwk
         ibuf(17,ixf)=puwkt
         ibuf(18,ixf)=iuwk
         ibuf(19,ixf)=juwk
         ibuf(20,ixf)=itrak
         ibuf(21,ixf)=ipt1
         ibuf(22,ixf)=ipt2
         ibuf(23,ixf)=iming
         ibuf(24,ixf)=imaxg
         ibuf(25,ixf)=jmin
         ibuf(26,ixf)=jmax
         ibuf(27,ixf)=jsmt
         ibuf(28,ixf)=ifsym
         ibuf(29,ixf)=jfsym
         ibuf(30,ixf)=2
c
      if(ixf.eq.maxbuf) then
       ibuf(1,1)=ixf
       ibuf(2,1)=iblock
       ibuf(3,1)=0
       call iosys('write real "'//path//'" to ciform without rewinding',
     $ 30*maxbuf,xbuf,0,' ')
       ixf=1
       noutf=noutf+1
      end if
c..bhl
c     call loopex
c..bhl
         goto 200
 480     continue
         if(pagein.and.nwks.gt.nwksmx) call pageout
         if(irstrt.ne.0) call save
         i=i-1
         if (i.ge.iming) go to 190
c
c..bhl
      if(ixf.gt.0) then
       ibuf(1,1)=ixf
       ibuf(2,1)=iblock
       ibuf(3,1)=-1
       call iosys('write real "'//path//'" to ciform without rewinding',
     $ 30*maxbuf,xbuf,0,' ')
       ixf=1
       noutf=noutf+1
      else
       ibuf(1,1)=-1
       ibuf(2,1)=iblock
       ibuf(3,1)=-1
       call iosys('write real "'//path//'" to ciform without rewinding',
     $ 30*maxbuf,xbuf,0,' ')
       ixf=1
       noutf=noutf+1
      end if
c..bhl
 500  continue
c
      if(ibhl.eq.1) then
       write(itape6,501) noutf
      else
       write(itape6,502) noutf
      end if
 501  format(/,3x,i8,' buffers of diagonal formulas were written')
 502  format(/,3x,i8,' buffers of ci       formulas were written')
c
      iblck1=1
      inxt=1
      return
      end
