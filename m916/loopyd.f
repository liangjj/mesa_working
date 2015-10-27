*deck @(#)loopyd.f	1.1  11/30/90
      subroutine loopyd(nabc,nlwks,           iwght,nabca,ijadd
     *,                 ijgrp,kadd,ladd,inext,iarc,jmnnxt,jmxnxt
     *,                 nabcs,isym,            isegm,jsegm,imain
     *,                 isub,iuwkmn,iuwksb,itrack,lmin,acoef,bcoef
     *,                 hdwgt,hdseg,hdpnt,imxgrp,imngrp,jmxgrp,jmngrp
     *,                 aint,file,unit,sigma,vector)
c
c
c
      implicit real*8 (a-h,o-z)
c
      integer xor
      character*16 file,unit
c..bhl
      character*16 path
c..bhl
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
c
      dimension sigma(*),vector(*)
c
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
      equivalence(ibuf(1,1),xbuf(1,1))
c
c..bhl
      common /nbhl/ ibhl
c..bhl
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
c
       ibhl=0
c
c--------------------------------------------------------------------loopy
c
      entry loopy
c
       ibhl=ibhl+1
       if(ibhl.eq.1) then
        path='diagonal form'
       else
        path='ci form'
       end if
c
       iout=itape6
c
      maxbuf=1000
c
      call iosys('rewind "'//file//'" on '//unit,0,0,0,' ')
      call iosys('rewind "'//path//'" on ciform',0,0,0,' ')
c
      do 500 iblock=1,ngroup
         call iosys('read real "'//file//'" from '//unit//
     #              ' without rewinding',nmax,aint,0,' ')
c..bhl
 498  continue
c
         call iosys('read real "'//path//'" from ciform '//
     #         ' without rewinding',30*maxbuf,xbuf,0,' ')
c
       ikx=ibuf(1,1)
       iblck=ibuf(2,1)
       last=ibuf(3,1)
c
        if(iblock.ne.iblck) call lnkerr(' loopy: iblck ne iblock ')
c
c..        write(itape6,*)'  ikx  last ',ikx,last
c
        if(ikx.gt.0) then
c
        do 499 i=2,ikx
         ladt=ibuf(1,i)
         acf=xbuf(2,i)
         d=xbuf(3,i)
         ccf=xbuf(4,i)
         itr1=ibuf(5,i)
         itr2=ibuf(6,i)
         ia=ibuf(7,i)
         ja=ibuf(8,i)
         itype=ibuf(9,i)
         isegt=ibuf(10,i)
         lvfrm1=ibuf(11,i)
         nlwki=ibuf(12,i)
         nlwkj=ibuf(13,i)
         imax=ibuf(14,i)
         imin=ibuf(15,i)
         nuwk=ibuf(16,i)
         puwkt=ibuf(17,i)
         iuwk=ibuf(18,i)
         juwk=ibuf(19,i)
         itrak=ibuf(20,i)
         ipt1=ibuf(21,i)
         ipt2=ibuf(22,i)
         iming=ibuf(23,i)
         imaxg=ibuf(24,i)
         jmin=ibuf(25,i)
         jmax=ibuf(26,i)
         jsmt=ibuf(27,i)
         ifsym=ibuf(28,i)
         jfsym=ibuf(29,i)
c
         index=ibuf(30,i)
c
90011    format(5(2x,f12.8))
c
         go to (1,2,3),index
   1     continue
c          write(iout,90012)nuwk,puwkt,iuwk,juwk,itrak,ipt1,ipt2
c          write(iout,90012)ladt,itr1,itr2,ia,ja,itype,isegt,lvfrm1
c          write(iout,90012)nlwki,nlwkj,imax,imin,iming,imaxg
c          write(iout,90012) jmin,jmax,jsmt,ifsym,jfsym
c          write(iout,90011) acf,d,ccf
90012     format(8(2x,i6))
c          write(iout,*)'  calling loopin  '
          call loopin
c       write(iout,90011)(sigma(ibl),ibl=1,136)
          go to 4
   2     continue
c          write(iout,*)'  calling loopex '
          call loopex
c       write(iout,90011)(sigma(ibl),ibl=1,136)
          go to 4
   3     continue
c          write(iout,*)'  calling allext '
          call allext
c       write(iout,90011)(sigma(ibl),ibl=1,136)
   4     continue
c         ibhl=ibhl+1
c         if(ibhl.gt.4) call lnkerr(' test stop ')
c
 499    continue
c
        end if
c
         if(last.eq.0) go to 498
c
 500  continue
      iblck1=1
      inxt=1
      return
      end
