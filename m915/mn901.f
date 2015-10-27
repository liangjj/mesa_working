*deck @(#)mn901.f	1.3  8/3/91
      subroutine mn901(ia,z,maxcor,operat,cvec,svec,tfile,tunit)
c
c***begin prologue     mn901
c***date written       yymmdd   (yymmdd)
c***revision date      880730   (yymmdd)
c
c     30 july   1988   bhl at llnl
c     setting nattim = nroots as the default
c
c     14 august 1987   pws at lanl
c     adding nroots-at-a-time as an option, and allowing
c     bliu to use rsp by expanding root to maxvec
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)mn901.f	1.3   8/3/91
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       mn901
c
      implicit integer (a-z)
c
      character*4096 ops
      character*(*) operat,tfile,tunit
      character*16 file,unit
      character*3 ians
      integer ia(*)
      logical pagein
      real*8 z(*),x,rep,fzcore,eguess,eci,cnverg,thresh,sqcdif
      real*8 cvec(*),svec(*)
      real*8 czero,drtver,civer
      real*8 time
      character*2 entry(20)
      character*8 prtflg
      integer nr3,nr6,nw6
c
      common /statpg/ ir3(100),nr3(100),ir6(100),nr6(100),iw6(100)
     #,               nw6(100),ipgs(100),ipgd(100),ipgo(100),ipgij(100)
     #,               ipgout(100)
      common /io/     itape5,itape6
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /optns/  iguess,irstrt,irooti,irootf,i34x
      common /all/  x(3),arr,ix(6),lvfrm1,nlwki,nlwkj,imax,imin
      common /symq/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
     #,             numsym(8)
      common /d901/ rep,fzcore,eguess,eci,refwlk,mxiter,cnverg,icnvg
     *,             iter,sqcdif,czero,nroot
      common /page/ iword3,iword4,pagein,noffi,noffj,npass,nwordl,iword6
      common /lbls/ lblint(26),lbldrt(26),civer,drtver
      common /rstrt/iblck1,inxt,lowtri,ndvdit,iblock
      common /count/  icount,ixx4,ixx5,ixx6,ixx8,iww4,iww5,iww6,iww8
     #,               iwx7,ixw9,ixy3,ixy16,ixy18,ixy22,iwy2,iwy15,iwy17
     #,               iwy21,iyx20,iyw19
      common /count2/ ientry(20),time(20)
      common /x4x901/ nijvir,nrefs
      common /nvectr/ nvc,mxvc
c
c
      equivalence (drtver,iver)
c
      data entry /'ww','xx','yy','xw','wx','wy','yx','yw','xy','zy'
     #,           'wz','xz','yz','4x','4x','4x','zz','4i','i1','20'/
c
      maxi=maxcor
      maxr=iadtwp(maxi)
      file=tfile
      unit=tunit
c
      do 2 i=1,20
         ientry(i)=0
         time(i)=0.0d+00
    2 continue
      icount=0
      icntd=0
c
c     write (3,'(' please enter maxr')')
c     read (3,*) maxr
      civer=3.290485
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
c
c     ----- start the timing routines -----
c
c
      nroots=intkey(ops,'ci=nroots',1,' ')
      mxiter=intkey(ops,'ci=cycles',15,' ')
      ntol=intkey(ops,'ci=tolerance',5,' ')
      ithrsh=intkey(ops,'ci=threshhold',5,' ')
      nguess=intkey(ops,'ci=nguess',nroots,' ')
      nattim=intkey(ops,'ci=nroots-at-a-time',nroots,' ')
c
   16 format(10i5)
      if(ntol.le.0) ntol=8
      cnverg=10.0d+00**(-ntol)
      thresh=10.0d+00**(-ithrsh)
      if(nroots.eq.0) nroots=1
      mxiter=mxiter*nattim+nguess
      ndvdmx=mxiter+iguess+1
      if(irooti.eq.0) irooti=1
      if(irootf.eq.0) irootf=irooti+nroots-1
c..bhl
          write(itape6,17)
  17      format(/,1x,'m915: ci formula link   ')
c..bhl
c      if(prtflg.ne.'minimum') then
c         write(itape6,17) nroots,nguess,mxiter,cnverg,thresh
c   17    format(1x,'m901:',
c     $         /5x,'number roots sought:          ',i8,
c     $         /5x,'number guess vectors:         ',i8,
c     $         /5x,'maximum number iterations:    ',i8,
c     $         /5x,'convergence criterion:        ',e8.1,
c     $         /5x,'davidson vector threshhold:   ',e8.1)
c..bhl
c        write(itape6,17) iguess,mxiter,irooti,nroots,ntol,irstrt,i34x,
c    $                    ilvfrm
c  17 format(1x,'m901:',
c    $          /5x,'guess:',i6,
c    $          /5x,'mxiter:',i6,
c    $          /5x,'rooti:',i6,
c    $          /5x,'nroot:',i6,
c    $          /5x,'ntol:',i6,
c    $          /5x,'rstrt:',i6,
c    $          /5x,'34x:',i6,
c    $          /5x,'fermi level:',i6)
c      end if
c
c     ----- get all the dimensions etc. needed -----
c
      call iosys('read integer "number of drt functions" from rwf',
     $     1,nbf,0,' ')
      call iosys('read integer "symmetries in ci" from rwf',
     $     1,nsym,0,' ')
      call iosys('read integer norbs from rwf',1,norbs,0,' ')
      call iosys('read integer nrows from rwf',1,nrows,0,' ')
      nrows4=nrows*4
      call iosys('read integer nlevs from rwf',1,nlevs,0,' ')
      call iosys('read integer nrefs from rwf',1,nrefs,0,' ')
      call iosys('read integer nwks from rwf',1,nwks,0,' ')
      call iosys('read integer orbfrm from rwf',1,orbfrm,0,' ')
      call iosys('read integer symorb from rwf',1,symorb,0,' ')
      call iosys('read integer numij from rwf',1,numij,0,' ')
      call iosys('read integer ngroup from rwf',1,ngroup,0,' ')
      call iosys('read integer nmax from rwf',1,nmax,0,' ')
      call iosys('read integer maxb from rwf',1,bmax,0,' ')
      call iosys('read integer nijvir from rwf',1,nijvir,0,' ')
c
      call iosys('read real "nuclear repulsion energy" from rwf',
     $     1,rep,0,' ')
c
      iout=itape6
c..bhl      write(iout,*)' mn901: nuclear repulsion energy ',rep
c
      call iosys('does "frozen core energy" exist on rwf',0,0,0,ians)
      if (ians.eq.'no') then
         fzcore=0.0d+00
      else
         call iosys('read real "frozen core energy" from rwf',
     #                    1,fzcore,0,' ')
      end if
c
      nroots=min(nroots,nwks)
c
c
      nmax2=intowp(nmax)
      nwks2=intowp(nwks)
      lvfrm1=orbfrm
      levfrm=lvfrm1+1
c
c
c      if (prtflg.ne.'minimum') then
c      write(itape6,32)nwks
c   32 format(5x,'number of configurations:     ',i8)
c      write(itape6,33)nbf
c   33 format(5x,'number of basis functions:    ',i8)
c      write(itape6,34) norbs
c   34 format(5x,'number of allowed orbitals:   ',i8)
c      end if
c
c     calculate pointers into blank common
c
      arc=1
c*    weight=arc+nrows4
      ijxx=arc+nrows4
      klxx=ijxx+numij
      nklxx=klxx+nijvir
      ijww=nklxx+nsym*orbfrm
      klww=ijww+numij
      nklww=klww+nijvir
      ist4x=nklww+nsym*orbfrm
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
      bcoef=acoef+nlevs
      root=bcoef+nlevs
      dvdvec=root+ndvdmx
      dvdmat=dvdvec+ndvdmx**2
      lowtri=ndvdmx*(ndvdmx+1)/2
      leni=wpadti(dvdmat+lowtri)
      need=iadtwp(leni)
      if (need.gt.maxi) then
         call getscm(need,z,maxr,'m901',0)
      end if
c
      call rzero(z,need)
c
      save=need
c
c
      call getdrt(ia(isegm),ia(aval),ia(bval)
     #,        ia(rowsym),ia(arc),ia(nlwks)
     #,        ia(weight),ia(wab),ia(wtw),ia(wtx),ia(wty)
     #,                   ia(kadd),ia(ladd),ia(ijadd),ia(ijgrp)
     #,        ia(inint),ia(inext),ia(jmnnxt),ia(jmxnxt),ia(ningrp)
     #,        ia(orbsym),ia(imngrp),ia(imxgrp),ia(jmngrp),ia(jmxgrp)
     #,        ia(ijxx),ia(klxx),ia(nklxx),ia(ijww),ia(klww),ia(nklww))
c
c
      nwksmx=nwks
      ln=0
      do 43 i=1,nsym
         if (minsym(i).ge.1000) numsym(i)=0
         if (minsym(i).lt.1000) numsym(i)=maxsym(i)-minsym(i)+1
         if (numsym(i).gt.ln) ln=numsym(i)
   43 continue
      lnsq=ln**2
      ints=need
      ci=ints+nmax
      si=ci+lnsq
      cj=si+lnsq
      sj=cj+lnsq
      h1=sj+lnsq
      h2=h1+lnsq
      sigma=h2+lnsq
c
      if (operat.eq.'h times c') then
         sigma=loc(svec)-loc(z)+1
         vector=loc(cvec)-loc(z)+1
         mxvc=1
         nvc=1
         lenr=h2+lnsq
      else
         if (sigma+2*nwks*min(64,nroots).gt.maxr) then
            call getscm(0,z,itop,'m901: size of core',0)
            mxvc=min(64,nroots,(itop-10-sigma)/(2*nwks))
         else
            mxvc=min(64,nroots)
         end if
         if (mxvc.lt.1) call lnkerr('not enough room to hold vectors')
         if (prtflg.ne.'minimum') then
            write (itape6,700) mxvc,nroots
  700       format(5x,'davidson vectors in use:      ',i8,
     $            /5x,'davidson vectors possible:    ',i8)
         end if
c
         vector=sigma+nwks*mxvc
         lenr=vector+nwks*mxvc
      end if
c
      if (lenr.gt.maxr) then
         call getscm(lenr+10,z,maxr,'m901',0)
      end if
      call rzero(z(save+1),lenr-save)
      if (operat.eq.'h times c') call rzero(svec,nwks)
c
      call initpg(z(vector),z(sigma))
c
      call initlp
c
      call initdg(z(ints),z(sigma),ia(ijadd),ia(kadd),ia(ladd)
     *,           ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym))
c
      call initex(z(ints),z(vector),z(sigma),ia(ijadd),ia(kadd)
     *,           ia(ladd),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym))
c
      call initx2(z(ints),z(vector),z(sigma),ia(ijadd),ia(kadd)
     *,           ia(ladd),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym))
c
      call initst(z(ints),z(vector),z(sigma),ia(ijadd),ia(kadd)
     *,           ia(ladd),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym)
     *,           z(ci),z(si),z(cj),z(sj),z(h1),z(h2)
     #,           ia(ijxx),ia(nklxx),ia(ijww),ia(nklww)
     #,           ia(klxx),ia(klww))
c
      call inimat(z(ints),ia(ijadd),ia(kadd),ia(ladd))
c
      call init4x(z(ints),z(vector),z(sigma),ia(ijadd),ia(kadd)
     *,           ia(ladd),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym))
c
      call start(z(dvdmat),z(ints),z(sigma),z(vector))
c
      call loopyd(ia(bval),ia(nlwks)                      ,ia(weight)
     *,           ia(aval),ia(ijadd),ia(ijgrp),ia(kadd),ia(ladd)
     *,           ia(inext),ia(arc),ia(jmnnxt),ia(jmxnxt),ia(rowsym)
     *,           ia(orbsym)                    ,ia(isegm)
     *,           ia(jsegm),ia(imain),ia(isub),ia(iuwkmn),ia(iuwksb)
     *,           ia(itrack),ia(lmin), z(acoef), z(bcoef)
     *,           ia(hdwgt),ia(hdseg),ia(hdpnt),ia(imxgrp),ia(imngrp)
     #,           ia(jmxgrp),ia(jmngrp),z(ints),file,unit)
c
       ndvdit=3
       nvc=1
       call loopy
c
       return
c
c
c     ----- error handling -----
c
  920 write(itape6,925)norbs,norb2
  925 format(' number of orbitals does not match ',2i10)
      stop
  930 write(itape6,935)leni
  935 format(' insufficient space in integer blank common need',i10,
     *' words.')
      stop
  940 write(itape6,945)4*nlwkmx+sigma
  945 format(' insufficient space in blank common need',i10,' words.')
      stop
      end
