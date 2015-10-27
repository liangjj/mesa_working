*deck @(#)pm911.f	5.1  11/6/94
      subroutine pm911(z,ia,maxi)
c
      implicit integer (a-z)
      integer wpadti
      integer ia(*)
      integer iii(3),ipciv(2)
      integer nr3,nr6,nw6
c
      integer nacpair(2,20)
c
      logical pagein
      logical debug
c
      real*8 z(*)
      real*8 x,rep,fzcore,eguess,eci,cnverg,sqcdif
      real*8 czero,drtver,civer,prtciv
      real*8 time,total,smaln,smalo
      real*8 cutof1,cutof2
c
      character*2 entry(20)
c
      common /statpg/ ir3(100),nr3(100),ir6(100),nr6(100),iw6(100)
     #,               nw6(100),ipgs(100),ipgd(100),ipgo(100),ipgij(100)
     #,               ipgout(100)
c     common /level/  ilev
c
      common /tapes/itap20,itape5,itape6,itape8,itap12,itap03,itap04
     *,             itape3,itap05,itap06
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,nrow4o,nwksoc,nlevoc,norboc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /optns/  iguess,irstrt,irooti,irootf,i34x
      common /all/  x(3),arr,ix(6),lvfrm1,nlwki,nlwkj,imax,imin
      common /symm/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
     #,             numsym(8)
      common /diag/ rep,fzcore,eguess,eci,cnverg,sqcdif,czero
     *,             refwlk,mxiter,icnvg,iter,nroot
      common /page/ iword3,iword4,pagein,noffi,noffj,npass,nwordl,iword6
      common /lbls/ lblint(26),lbldrt(26),civer,drtver
      common /rstrt/iblck1,inxt,lowtri,ndvdit,iblock
      common /count/  icount,ixx4,ixx5,ixx6,ixx8,iww4,iww5,iww6,iww8
     #,               iwx7,ixw9,ixy3,ixy16,ixy18,ixy22,iwy2,iwy15,iwy17
     #,               iwy21,iyx20,iyw19
      common /count2/ ientry(20),time(20)
      common/count3/miseg(40,20),mxm,mxseg,msoutr
c*
      common /xx4xbl/   nijvir
      common/dry2/prtciv,navail,nrused(2)
      common/dry1/idry,jdry,ksegsv,indxs(100),magic
c
      common /matrix/ cutof1,cutof2,nmultp,nmultd
      common /io/ inp,iout
      character*3 answer
      character*128 gden
c
      character*4096 ops
      logical logkey
c
      data itoln,itolo/8,4/, ipciv/1,1/,itrdm/0/,ioptci/1/
      data entry /'ww','xx','yy','xw','wx','wy','yx','yw','xy','zy'
     #,           'wz','xz','yz','4x','4x','4x','zz','4i','i1','20'/
c
      save itoln,itolo,ipciv,itrdm,ioptci,entry
c

      parameter (debug=.false.)
c
      write(iout,10)
   10 format(1x,'m911:density matrix construction')
c
c..bhl
      call iosys('read character options from rwf',-1,0,0,ops)
c
      if(logkey(ops,'nacme',.false.,' ')) then
         write(iout,1000)
 1000    format(5x,'nacme run ')
         itrdm=1
         nacpair(1,1)=1
         nacpair(2,1)=2
         call iosys('does "nacme counter" exist on rwf',0,0,0,answer) 
         if(answer.eq.'no') then
            write(iout,1010)
 1010       format(5x,'initializing the nacme counter on the rwf ')
            knacme=1
            numnac=intkey(ops,'nacme=couplings',1,' ')
            call intarr(ops,'nacme=states',nacpair,2*numnac,' ')
            call iosys('write integer "nacme pairs" to rwf',
     $                 1,numnac,0,' ')
            call iosys('write integer "nacme states" to rwf',
     $                 2*numnac,nacpair,0,' ')
            state1=nacpair(1,knacme)
            state2=nacpair(2,knacme)
         else
            call iosys('read integer "nacme pairs" from rwf',
     $                 1,numnac,0,' ')
            call iosys('read integer "nacme counter" from rwf',
     $                 1,knacme,0,' ')
            knacme=knacme+1
            write(iout,1050) knacme,numnac
 1050       format(5x,'nacme counter  :',i5,/5x,'total couplings:',i5)
            if(knacme.gt.numnac) then
               call chainx(0)
               stop
            end if
            call iosys('read integer "nacme states" from rwf',
     $                 2*numnac,nacpair,0,' ')
            state1=nacpair(1,knacme)
            state2=nacpair(2,knacme)
         end if 
         write(iout,1020) state1,state2
 1020    format(5x,'density being computed for states:',2i5)
      else
         itrdm=0
      end if
c
      maxr=intowp(maxi)
      cutof1=10.0d+00**(-10)
      cutof2=10.0d+00**(-11)
      nmultp=0
      nmultd=0
      itape5=inp
      itape6=iout
c
      do 1895 i=1,20
      do 1895 j=1,40
 1895    miseg(j,i)=0
      mxm=0
      mxseg=0
      msoutr=0
c
      do 2 i=1,20
         ientry(i)=0
         time(i)=0.0d+00
    2 continue
      icount=0
      icntd=0
c
      if(debug) then
         write(iout,*)' maxi ',maxi
      endif
c
      do 1313 i=1,maxi
         ia(i)=0
 1313 continue
c
c     write(iout,5)
c     write(iout,15) maxr
c  15 format(1h0,18x,'dimension of real common is',i8)
c
c     -----read in options for this run -----
c
c
      navail=2
      nrused(1)=1
      nrused(2)=2
      magic=0
      ipflg=0
c
      prtciv=float(ipciv(1)) * 10.d0 ** (-ipciv(2) )
      smaln=10.d0**(-itoln)
      smalo=10.d0**(-itolo)
      iguess=1
      mxiter=0
      irooti=0
      nroots=0
      ntol=0
      irstrt=0
      i34x=0
      ilvfrm=0
   16 format(10i5)
      if (iguess.le.0) iguess=1
      if(ntol.le.0) ntol=8
      cnverg=10.0d+00**(-ntol)
      if(nroots.eq.0) nroots=1
      if(mxiter.eq.0) mxiter=10
      mxiter=mxiter*nroots
      ndvdmx=mxiter+iguess+1
      if(irooti.eq.0) irooti=1
      if(irootf.eq.0) irootf=irooti+nroots-1
c     write(iout,17) navail,nrused
c  17 format(1h0,'options',11x,'navail nrused:',3i5)
c
c     write(iout,20)
c  20 format(1h0,'read in distinct row table')
c
c     ----- get the constants for dividing up core -----
c
      call iosys('read integer ngroup from rwf',1,ngroup,0,' ')
      call iosys('read integer "symmetries in ci" from rwf',
     $            1,nsym,0,' ')
      call iosys('read integer norbs from rwf',1,norbs,0,' ')
      call iosys('read integer "number of basis functions" from rwf',
     #            1,nbf,0,' ')
      call iosys('read integer nrows from rwf',1,nrows,0,' ')
      call iosys('read integer nlevs from rwf',1,nlevs,0,' ')
      call iosys('read integer nrefs from rwf',1,nrefs,0,' ')
      call iosys('read integer symorb from rwf',1,symorb,0,' ')
      call iosys('read integer nwks from rwf',1,nwks,0,' ')
      call iosys('read integer orbfrm from rwf',1,orbfrm,0,' ')
      call iosys('read integer numij from rwf',1,numij,0,' ')
      call iosys('read integer nmax from rwf',1,nmax,0,' ')
      call iosys('read integer nijvir from rwf',1,nijvir,0,' ')
      call iosys('read integer maxb from rwf',1,bmax,0,' ')
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,rep,0,' ')
      nrows4=4*nrows
c
      call iosys('read character "guga density filename" from rwf',
     $            0,0,0,gden)
      call iosys('open gden as unknown',0,0,0,gden)
      call iosys('does "guga density matrix" exist on gden',0,
     #            0,0,answer)
      if(answer.eq.'no') then
         write(iout,20)
   20 format(5x,'creating guga density matrix')
         call iosys('create real "guga density matrix" on gden',
     1               ngroup*nmax,0,0,' ')
      end if
c
      call iosys('rewind "guga density matrix" on gden',0,0,0,' ')
c
      call iosys('does "guga square 1pdm" exist on rwf',0,0,0,answer)
      if(answer.eq.'no') then
         call iosys('create real "guga square 1pdm" on rwf',nbf*nbf,
     1               0,0,' ')
      end if
c
      call iosys('rewind "guga square 1pdm" on rwf',0,0,0,' ')
c
c
      nmax2=nmax+nmax
      nwks2=nwks+nwks
      lvfrm1=orbfrm
      levfrm=lvfrm1+1
c
c
      if (ilvfrm.ne.0) levfrm=ilvfrm
c
c
      norboc=nlevoc-1
      if(irstrt.gt.0) go to 30
c
c
   30 write(iout,31)lbldrt
   31 format(5x,'label from drt...',26a3)
      write(iout,32)nwks
   32 format(5x,'number of configurations is',i11)
      write(iout,33)nbf
   33 format(5x,'number of basis functions is',i10)
      write(iout,34) norbs
   34 format(5x,'number of allowed orbitals is',i9)
      write(iout,35)lvfrm1
   35 format(5x,'fermi level is ',i23)
c
      iii(2)=norbs
c     calculate pointers into blank common
c
      arc=1
c*    weight=arc+nrows4
      ijxx=arc+nrows4
      klxx=ijxx+numij
      nklxx=klxx+nijvir
      ijww=nklxx+nsym
      klww=ijww+numij
      nklww=klww+nijvir
      ist4x=nklww+nsym
      len4x=ist4x+norbs
      weight=len4x+norbs
c*
      wtw=weight+nrows4
      wtx=wtw+nsym*lvfrm1
      wty=wtx+nsym*lvfrm1
      wab=wty+lvfrm1
      nlwks=wab+lvfrm1
      aval=nlwks+nrows
      bval=aval+nrows
      rowsym=bval+nrows
      ijadd=rowsym+nrows
      ijgrp=ijadd+numij
      kadd=ijgrp+numij
      ladd=kadd+symorb
      jmxnxt=ladd+symorb
      jmnnxt=jmxnxt+norbs
      inint=jmnnxt+norbs
      inext=inint+norbs
      ningrp=inext+norbs
      imngrp=ningrp+ngroup
      imxgrp=imngrp+ngroup
      jmngrp=imxgrp+ngroup
      jmxgrp=jmngrp+ngroup
      orbsym=jmxgrp+ngroup
      isegm=orbsym+norbs
      jsegm=isegm+nlevs
      imain=jsegm+nlevs
      isub=imain+nlevs
      iuwkmn=isub+nlevs
      iuwksb=iuwkmn+nlevs
      itrack=iuwksb+nlevs
      lmin=itrack+nlevs
      hdwgt=lmin+nlevs
      hdseg=hdwgt+nlevs
      hdpnt=hdseg+nlevs
      acoef=iadtwp(hdpnt+nlevs)
      bcoef=acoef+nlevs*2
      root=bcoef+nlevs*2
      dvdvec=root+nroots*2
      dvdmat=dvdvec+(nroots+1)*ndvdmx*2
      lowtri=ndvdmx*(ndvdmx+1)
      leni=wpadti(dvdmat+lowtri)
      lowtri=lowtri/2
      if(debug) then
         write(iout,41) leni
   41    format(' length of integer common is ',i8)
      endif
      if(leni.gt.maxi) then
         write(iout,935)leni
  935    format(' insufficient space in integer blank common need',i10,
     *          ' words.')
         call lnkerr(' m911: insufficient memory')
      end if
      if(debug) then
         write(iout,42)
   42    format(' read drt file')
      endif
      call getdrt(ia(isegm),ia(aval),ia(bval)
     #,        ia(rowsym),ia(arc),ia(nlwks)
     #,        ia(weight),ia(wab),ia(wtw),ia(wtx),ia(wty)
     #,                   ia(kadd),ia(ladd),ia(ijadd),ia(ijgrp)
     #,        ia(inint),ia(inext),ia(jmnnxt),ia(jmxnxt),ia(ningrp)
     #,        ia(orbsym),ia(imngrp),ia(imxgrp),ia(jmngrp),ia(jmxgrp)
     #,        ia(ijxx),ia(klxx),ia(nklxx),ia(ijww),ia(klww),ia(nklww))
      nwksmx=nwks
      ln=0
      do 43 i=1,nsym
         if (minsym(i).ge.1000) numsym(i)=0
         if (minsym(i).lt.1000) numsym(i)=maxsym(i)-minsym(i)+1
         if (numsym(i).gt.ln) ln=numsym(i)
   43 continue
c
      lnsq=ln**2
      norb2=norbs**2
      ints=iadtwp(leni)
      ints1=ints+nmax
      ci=ints1+norb2
      call rzero(z(ints1),norb2)
      si=ci+lnsq
      cj=si+lnsq
      sj=cj+lnsq
      h1=sj+lnsq
      h2=h1+lnsq
      sigma=h2+lnsq
      vector=sigma+nwks
      lenr=vector+nwks
c
      if(debug) then
         write (iout,90911) lenr,maxr
90911    format(/,' need',i8,' real words of core to avoid'
     #,           ' paging',/,' have',i8,' real words.')
      endif
c
      if(lenr.ge.maxr) then
         write(iout,45)
   45    format(1h ,'insufficient core to hold all vectors. paging will'
     *             ,' be tried,')
         nwksmx=(maxr-sigma)/2
         if(2*nlwkmx.gt.nwksmx) then
            write(iout,945)4*nlwkmx+sigma
  945       format(' insufficient space in blank common need',
     $               i10,' words.')
            call lnkerr(' m911: insufficient memory')
         end if
         write(iout,46) nwksmx
  46     format(1x,i8,' elements of the vector can be held at once.')
         vector=sigma+nwksmx
         lenr=vector+nwksmx
         write (iout,49) lenr
  49     format (1x,' amount of real core used:',i8)
      end if
c pws
      npass=(nwks-1)/nwksmx+1
      nwkmx2=intowp(nwksmx)
      nwordl=intowp(nwks)-(npass-1)*nwkmx2
c      write (iout,80802) nwordl
c80802 format(' nwordl:',i8)
c
      call initnx(indxs,norbs)
      call initpg(z(vector),z(sigma))
      call initlp(z(ints),z(ints1),z(vector),z(sigma),ia(ijadd),
     $            ia(kadd),ia(ladd),ia(wtw),ia(wtx),ia(wty),
     $            ia(wab),ia(orbsym),
     $            z(ci),z(si),z(cj),z(sj),z(h1),z(h2),
     $            ia(ijxx),ia(nklxx),ia(ijww),ia(nklww),
     $            ia(klxx),ia(klww))
      call initdg(z(ints),z(ints1),z(vector),z(sigma),ia(ijadd),
     $            ia(kadd),ia(ladd),ia(wtw),ia(wtx),ia(wty),
     $            ia(wab),ia(orbsym))
      call initex(z(ints),z(ints1),z(vector),z(sigma),ia(ijadd),
     $            ia(kadd),ia(ladd),ia(wtw),ia(wtx),ia(wty),
     $            ia(wab),ia(orbsym))
      call initx2(z(ints),z(ints1),z(vector),z(sigma),ia(ijadd),
     $            ia(kadd),ia(ladd),ia(wtw),ia(wtx),ia(wty),
     $            ia(wab),ia(orbsym),
     $            z(ci),z(si),z(cj),z(sj),z(h1),z(h2),
     $            ia(ijxx),ia(nklxx),ia(ijww),ia(nklww))
      call initst(z(ints),z(vector),z(sigma),ia(ijadd),ia(kadd),
     $            ia(ladd),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym),
     $            z(ci),z(si),z(cj),z(sj),z(h1),z(h2),
     $            ia(ijxx),ia(nklxx),ia(ijww),ia(nklww))
      call inimat(z(ints),ia(ijadd),ia(kadd),ia(ladd))
      call init4x(z(ints),z(ints1),z(vector),z(sigma),ia(ijadd),
     $            ia(kadd),ia(ladd),ia(wtw),ia(wtx),ia(wty),
     $            ia(wab),ia(orbsym))
c     call start(ia(dvdmat),z(ints),z(sigma),z(vector))
c     if(irstrt.le.0) write(iout,70)
c  70 format(1h0,'calculate diagonal elements')
      call loopyd(ia(bval),ia(nlwks)                      ,ia(weight)
     *,           ia(aval),ia(ijadd),ia(ijgrp),ia(kadd),ia(ladd)
     *,           ia(inext),ia(arc),ia(jmnnxt),ia(jmxnxt),ia(rowsym)
     *,           ia(orbsym)                    ,ia(isegm)
     *,           ia(jsegm),ia(imain),ia(isub),ia(iuwkmn),ia(iuwksb)
     *,           ia(itrack),ia(lmin),z(acoef),z(bcoef)
     *,           ia(hdwgt),ia(hdseg),ia(hdpnt),ia(imxgrp),ia(imngrp)
     #,           ia(jmxgrp),ia(jmngrp)
     #,           z(ints),z(ints1),file,unit,z(vector),z(sigma),
     #            ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym),
     #            z(ci),z(si),z(cj),z(sj),z(h1),z(h2),
     #            ia(ijxx),ia(nklxx),ia(ijww),ia(nklww),
     #            ia(klxx),ia(klww))
c     write(iout,5)
c     write(iout,1)
      call fdiag(z(vector),z(sigma),smaln,smalo,itrdm,ioptci,
     #           state1,state2)
c
c     ----- form two-particle density matrix -----
c
      call loopy(ia(bval),ia(nlwks)                      ,ia(weight)
     *,           ia(aval),ia(ijadd),ia(ijgrp),ia(kadd),ia(ladd)
     *,           ia(inext),ia(arc),ia(jmnnxt),ia(jmxnxt),ia(rowsym)
     *,           ia(orbsym)                    ,ia(isegm)
     *,           ia(jsegm),ia(imain),ia(isub),ia(iuwkmn),ia(iuwksb)
     *,           ia(itrack),ia(lmin),z(acoef),z(bcoef)
     *,           ia(hdwgt),ia(hdseg),ia(hdpnt),ia(imxgrp),ia(imngrp)
     #,           ia(jmxgrp),ia(jmngrp)
     #,           z(ints),z(ints1),file,unit,z(vector),z(sigma),
     #            ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym),
     #            z(ci),z(si),z(cj),z(sj),z(h1),z(h2),
     #            ia(ijxx),ia(nklxx),ia(ijww),ia(nklww),
     #            ia(klxx),ia(klww))
c..rlm     call loopy
c
c     call trace(z(ints1),norbs)
c      if(iand(ipflg,1).ne.0)then
c      write(iout,70701)
c70701 format(/5x,'from dm: square 1-particle dm',/)
c     call prtsm(z(ints1),norbs,iout)
c      end if
c
c
c     write(iout,90)
c  90 format(1h0,5x,'timing information')
cc final timing values are returned here
c     write (iout,891) icount,ixx4,ixx5,ixx6,ixx8,iww4,iww5,iww6,iww8
c    #,                  iwx7,ixw9,ixy3,ixy16,ixy18,ixy22,iwy2,iwy15
c    #,                  iwy17,iwy21,iyx20,iyw19
  891 format (//,' icount is:',i8,/
     #,                        t15,'xx entries, iseg 4, 5, 6, 8:',4i8,/
     #,                        t15,'ww entries, iseg 4, 5, 6, 8:',4i8,/
     #,                        t15,'wx entries, iseg 7         :',i8,/
     #,                        t15,'xw entries, iseg 9         :',i8,/
     #,                        t15,'xy entries, iseg 3,16,18,22:',4i8,/
     #,                        t15,'wy entries, iseg 2,15,17,21:',4i8,/
     #,                        t15,'yx entries, iseg 20        :',i8,/
     #,                        t15,'yw entries, iseg 19        :',i8)
c     write(iout,1891)mxm,mxseg,msoutr
c     do 1892 i=1,mxm
c1892    write(iout,1893)i,(miseg(j,i),j=1,mxseg)
c1891 format(2x,'from initex:max m=',i5,' max iseg=',i5,' out of ran='
c    1         ,i3)
c1893 format(5x,'m=',i3,/(25i4))
c     total=0.0d+00
c     do 95 i=1,20
c        total=total+time(i)
c  95 continue
c     write (iout,892) total,(entry(i),ientry(i),time(i)
c    #,                   time(i)/total*100,i=1,20)
c 892 format(1h1,'  entries and times in seconds',f10.3,//,(t5,a2,i8
c    #,       f10.3,f7.2))
cc
c     write (iout,893) (i-1,ipgs(i),ipgd(i),ipgo(i),ipgij(i),ipgout(i)
c    #,                   i=nlevs,2,-1)
c 893 format(1h1,"  paging summary: calls to routines by levels",//
c    #,     "     level   pages   paged   pageo  pageij  pageout"
c    #,     //,(6i8))
c     write (iout,894) (i-1,ir3(i),nr3(i)/2,nr3(i)/(ir3(i)*2)
c    #,                      ir6(i),nr6(i)/2,nr6(i)/(ir6(i)*2)
c    #,                      iw6(i),nw6(i)/2,nw6(i)/(iw6(i)*2)
c    #,                   i=nlevs,2,-1)
c 894 format(1h1,"  paging summary: real numbers transferred per level"
c    #,//,"                 vector reads                    sigma reads"
c    #,         "------------------------writes"
c    #,   /,"  level  ",3("     number     #transferred      average  ")
c    #,       /,(i5,9i14))
c     do 895 i=2,nlevs
c     ir3(1)=ir3(1)+ir3(i)
c     nr3(1)=nr3(1)+nr3(i)/2
c     ir6(1)=ir6(1)+ir6(i)
c     nr6(1)=nr6(1)+nr6(i)/2
c     iw6(1)=iw6(1)+iw6(i)
c     nw6(1)=nw6(1)+nw6(i)/2
c 895 continue
c     write (iout,896) ir3(1),nr3(1),ir6(1),nr6(1),iw6(1),nw6(1)
c 896 format (/,5x,3(2i14,14x))
c     call etime
c     call wtime
c     write(iout,5)
c
c     ----- exit gracefully -----
      call iosys('close gden',0,0,0,' ')
c
      if(itrdm.ne.0) then
         write(iout,1040)
 1040    format(5x,'writing the nacme counter to the rwf ')
         call iosys('write integer "nacme counter" to rwf',
     #               1,knacme,0,' ')
         call chainx(1)
         stop
      end if
c
      call chainx(0)
c
c
      return
      end
