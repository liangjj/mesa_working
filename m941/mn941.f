*deck @(#)mn941.f	5.1  11/6/94
      subroutine mn941(ops,operat,cvec,svec,tfile,tunit)
c
c***begin prologue     mn941
c***date written       yymmdd   (yymmdd)
c***revision date      880730   (yymmdd)
c
c     30 july   1988   bhl at llnl
c       setting nattim=nroots as the default.
c
c     14 august 1987   pws at lanl
c       adding nroots-at-a-time as an option, and allowing
c       bliu to use rsp by expanding root to maxvec
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)mn941.f	5.1  11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       mn941
c
      implicit integer (a-z)
c
      character*(*) ops
      character*(*) operat,tfile,tunit
      character*16 file,unit
      character*16 chrkey,gestyp
      character*8 prtflg
      character*3 ians
      logical pagein
      dimension ngot(5)
      real*8 za, zb, dvd
      integer ia, ib, idvd
      real*8 x,rep,fzcore,eguess,eci,cnverg,thresh,sqcdif
      real*8 cvec(*),svec(*)
      real*8 czero
      pointer(pa,za(1)),(pa,ia(1))
      pointer(pb,zb(1)),(pb,ib(1))
      pointer(pdvd,dvd(1)),(pdvd,idvd(1))
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
      common /rstrt/iblck1,inxt,lowtri,ndvdit,iblock
      common /x4x901/ nijvir,nrefs
      common /nvectr/ nvc,mxvc
c
c
c
c     find maximum number of integer words available.  this
c     is a parameter which is set in stone in the code.  it
c     can be changed if the physical memory of the machine is
c     larger.
c
      call iosys('read integer mxcore from rwf',1,maxcor,0,' ')
      maxi=maxcor
      maxr=iadtwp(maxi)
      file=tfile
      unit=tunit
c
      icount=0
      icntd=0
c
c
      call iosys('read character "print flag" from rwf',
#                 -1,0,0,prtflg)
c
c
c     recover ci control parameters
      call iosys('read integer nwks from rwf',1,nwks,0,' ')
      call iosys('read integer norbs from rwf',1,norbs,0,' ')
      call iosys('read integer "number of drt functions" from rwf',
     $     1,nbf,0,' ')
      nroots=intkey(ops,'ci=number-of-roots',1,' ')
      nroots=min(nroots,nwks)
      mxiter=intkey(ops,'ci=iterations',15,' ')
      thresh=intkey(ops,'ci=tolerance',1.d-10,' ')
      cnverg=intkey(ops,'ci=convergence',1.d-08,' ')
      nwksg=intkey(ops,'ci=guess-size',0,' ')
      nattim=intkey(ops,'ci=nroots-at-a-time',nroots,' ')
      gestyp=chrkey(ops,'ci=guess=','diagonals',' ')
      if(nroots.eq.0) then
         nroots=1
         nattim=1
      endif
      if(mxiter.lt.2*nattim) then
         mxiter=2*nattim
      endif
      ndvdmx=mxiter+iguess+1
c
   16 format(10i5)
      if(irooti.eq.0) irooti=1
      if(irootf.eq.0) irootf=irooti+nroots-1
c
      if(prtflg.ne.'minimum') then
         write(itape6,17) nwks,nbf,norbs,nroots,mxiter,cnverg,thresh
   17    format(1x,'m901:',
     $         /5x,'number configurations      = ',i8,
     $         /5x,'number orbitals            = ',i8,
     $         /5x,'number orbitals allowed    = ',i8,
     $         /5x,'number roots sought        = ',i8,
     $         /5x,'maximum number iterations  = ',i8,
     $         /5x,'convergence criterion      = ',e8.1,
     $         /5x,'davidson vector threshhold = ',e8.1)
      end if
c
c     ----- get all the dimensions etc. needed -----
c

      call iosys('read integer "symmetries in ci" from rwf',
     $     1,nsym,0,' ')
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
c
      call iosys('does "frozen core energy" exist on rwf',0,0,0,ians)
      if (ians.eq.'no') then
         fzcore=0.0d+00
      else
         call iosys('read real "frozen core energy" from rwf',
     #                    1,fzcore,0,' ')
      end if
c
c
c
      nmax2=intowp(nmax)
      nwks2=intowp(nwks)
      lvfrm1=orbfrm
      levfrm=lvfrm1+1
c
c
c
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
      need=wpadti(bcoef+nlevs)
      maxi=maxi-need
c
c     get the memory for the pointers
c
      call getmem(need,pa,ngot(1),'m941-pointers',0)
c
c
c
      call rzero(za,ngot(1))
c
      save=need
c
c
      call getdrt(ia(isegm),ia(aval),ia(bval),
#                 ia(rowsym),ia(arc),ia(nlwks),
#                 ia(weight),ia(wab),ia(wtw),
#                 ia(wtx),ia(wty),ia(kadd),
#                 ia(ladd),ia(ijadd),ia(ijgrp),
#                 ia(inint),ia(inext),ia(jmnnxt),
#                 ia(jmxnxt),ia(ningrp),ia(orbsym),
#                 ia(imngrp),ia(imxgrp),ia(jmngrp),
#                 ia(jmxgrp),ia(ijxx),ia(klxx),
#                 ia(nklxx),ia(ijww),ia(klww),
#                 ia(nklww))
c
c     now get the memory needed for the integrals and the 
c     init routines.
c
      nwksmx=nwks
      ln=0
      do 43 i=1,nsym
         if (minsym(i).ge.1000) numsym(i)=0
         if (minsym(i).lt.1000) numsym(i)=maxsym(i)-minsym(i)+1
         if (numsym(i).gt.ln) ln=numsym(i)
   43 continue
      lnsq=ln**2
      ints=1
      ci=ints+nmax
      si=ci+lnsq
      cj=si+lnsq
      sj=cj+lnsq
      h1=sj+lnsq
      h2=h1+lnsq
      need=wpadti(h2+lnsq)
      call getmem(need,pb,ngot(2),'m941-ints',0)
      maxi=maxi-need
c
c     get the memory needed for the davidson matrices
c
      dvdeig=1
      dvdvec=dvdeig + ndvdmx
      dvdmat=dvdvec+ndvdmx**2
      work=dvdmat+ndvdmx**2
      wdvdeig=work+5*ndvdmx
      wdvdvec=wdvdeig + ndvdmx
      wdvdmat=wdvdvec+ndvdmx*ndvdmx
      wvec=wdvdmat+ndvdmx*ndvdmx
      need=wpadti(wvec+ndvdmx*ndvdmx)
      call getmem(need,pdvd,ngot(3),'m941-dvd',0)
      maxi=maxi-need
c
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
            itop=iadtwp(itop)
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
      call initpg
c
c..rlm     call initlp
      call initlp(zb(ints),z(vector),z(sigma),
     $            ia(ijadd),ia(kadd),ia(ladd),
     $            ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym),
     $            zb(ci),zb(si),zb(cj),zb(sj),zb(h1),zb(h2),
     $            ia(ijxx),ia(nklxx),ia(ijww),ia(nklww),
     $            ia(klxx),ia(klww))
c
c..rlm     call initdg
      call initdg(zb(ints),z(sigma),ia(ijadd),ia(kadd),ia(ladd),
     $            ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym))
c
c..rlm      call initex
      call initex(zb(ints),z(vector),z(sigma),ia(ijadd),ia(kadd),
     $           ia(ladd),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym))
c
c..rlm     call initx2
      call initx2(zb(ints),z(vector),z(sigma),ia(ijadd),ia(kadd),
     *            ia(ladd),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym))
c
c..rlm     call initst
      call initst(zb(ints),z(vector),z(sigma),ia(kadd),
     $           ia(ladd),ia(wtw),ia(wtx),ia(wab),ia(orbsym),
     $           zb(ci),zb(si),zb(cj),zb(sj),zb(h1),zb(h2),
     $           ia(ijxx),ia(nklxx),ia(ijww),ia(nklww),
     $           ia(klxx),ia(klww))
c
cpws      call inimat(zb(ints),ia(ijadd),ia(kadd),ia(ladd))
c
c..rlm     call init4x
      call init4x(zb(ints),z(vector),z(sigma),ia(ijadd),ia(kadd),
     $           ia(ladd),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym))
c
      call start(z(dvdmat),zb(ints),z(sigma),z(vector))
c
      call loopyd(ia(bval),ia(nlwks)                      ,ia(weight),
     $            ia(aval),ia(ijadd),ia(ijgrp),ia(kadd),ia(ladd),
     $            ia(inext),ia(arc),ia(jmnnxt),ia(jmxnxt),ia(rowsym),
     $            ia(orbsym)                    ,ia(isegm),
     $            ia(jsegm),ia(imain),ia(isub),ia(iuwkmn),ia(iuwksb),
     $            ia(itrack),ia(lmin), z(acoef), z(bcoef),
     $            ia(hdwgt),ia(hdseg),ia(hdpnt),ia(imxgrp),ia(imngrp),
     $            ia(jmxgrp),ia(jmngrp),zb(ints),file,unit,z(vector),
     $            z(sigma),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym),
     $            zb(ci),zb(si),zb(cj),zb(sj),zb(h1),zb(h2),
     $            ia(ijxx),ia(nklxx),ia(ijww),ia(nklww),
     $            ia(klxx),ia(klww))
c
c
      if (operat.eq.'h times c') then
c
c        ----- multiply vector in cvec by diagonals in svec -----
c
         call iosys('write real h_diagonals to rwf',nwks,svec,0,' ')
c
         do 345 i=1,nwks
            svec(i)=svec(i)*cvec(i)
  345    continue
c
         ndvdit=3
         nvc=1
         call loopy(ia(bval),ia(nlwks),ia(weight)
     *,           ia(aval),ia(ijadd),ia(ijgrp),ia(kadd),ia(ladd)
     *,           ia(inext),ia(arc),ia(jmnnxt),ia(jmxnxt),ia(rowsym)
     *,           ia(orbsym),ia(isegm)
     *,           ia(jsegm),ia(imain),ia(isub),ia(iuwkmn),ia(iuwksb)
     *,           ia(itrack),ia(lmin), z(acoef), z(bcoef)
     *,           ia(hdwgt),ia(hdseg),ia(hdpnt),ia(imxgrp),ia(imngrp)
     #,           ia(jmxgrp),ia(jmngrp),zb(ints),file,unit,
     #            z(vector),z(sigma),ia(wtw),ia(wtx),ia(wty),
     #            ia(wab),ia(orbsym)
     *,           zb(ci),zb(si),zb(cj),zb(sj),zb(h1),zb(h2)
     #,           ia(ijxx),ia(nklxx),ia(ijww),ia(nklww)
     #,           ia(klxx),ia(klww))
      else
c
c        initialize the diagonalization routine and generate the guess vectors.
         call fdiag(zb(ints),z(vector),z(sigma),z(dvdmat),z(root),
     #           z(dvdvec),ndvdmx,thresh,nguess,prtflg,nattim,
     #           ia(bval),ia(nlwks),ia(weight)
     *,          ia(aval),ia(ijadd),ia(ijgrp),ia(kadd),ia(ladd)
     *,          ia(inext),ia(arc),ia(jmnnxt),ia(jmxnxt),ia(rowsym)
     *,          ia(orbsym),ia(isegm)
     *,          ia(jsegm),ia(imain),ia(isub),ia(iuwkmn),ia(iuwksb)
     *,          ia(itrack),ia(lmin), z(acoef), z(bcoef)
     *,          ia(hdwgt),ia(hdseg),ia(hdpnt),ia(imxgrp),ia(imngrp)
     #,          ia(jmxgrp),ia(jmngrp),file,unit,ia(wtw),ia(wtx),
     #           ia(wty),ia(wab),ia(orbsym)
     *,          zb(ci),zb(si),zb(cj),zb(sj),zb(h1),zb(h2)
     #,          ia(ijxx),ia(nklxx),ia(ijww),ia(nklww)
     #,          ia(klxx),ia(klww),gestyp)
c
c        iterate the diagonalization procedure.
         call roots(zb(ints),z(vector),z(sigma),z(root),z(dvdvec),
     #           z(dvdmat),ndvdmx,prtflg,
     #           ia(bval),ia(nlwks),ia(weight)
     *,          ia(aval),ia(ijadd),ia(ijgrp),ia(kadd),ia(ladd)
     *,          ia(inext),ia(arc),ia(jmnnxt),ia(jmxnxt),ia(rowsym)
     *,          ia(orbsym),ia(isegm)
     *,          ia(jsegm),ia(imain),ia(isub),ia(iuwkmn),ia(iuwksb)
     *,          ia(itrack),ia(lmin), z(acoef), z(bcoef)
     *,          ia(hdwgt),ia(hdseg),ia(hdpnt),ia(imxgrp),ia(imngrp)
     #,          ia(jmxgrp),ia(jmngrp),file,unit,ia(wtw),ia(wtx),
     #           ia(wty),ia(wab),ia(orbsym)
     *,          zb(ci),zb(si),zb(cj),zb(sj),zb(h1),zb(h2)
     #,          ia(ijxx),ia(nklxx),ia(ijww),ia(nklww)
     #,          ia(klxx),ia(klww))
c
         call iosys('destroy bliu',0,0,0,' ')
      end if
c
c
      return
      end
