*deck @(#)mn901.f	5.1  11/6/94
      subroutine mn901(tfile,tunit)
c
c***begin prologue     mn901
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
c***source             @(#)mn901.f	5.1  11/6/94
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
      character*(*) tfile,tunit
      character*16 file,unit, key
      character*16 chrkey,gestyp
      character*8 prtflg, dsk
      character*3 ians
      character*128 nmcnfg
      integer ia
      integer ngot(3)
      logical pagein
      real*8 x,rep,fzcore,eguess,eci,cnverg,thresh,sqcdif
      real*8 czero
      real*8 za, z, y
      pointer (p,ia(1)), (p,za(1))
      pointer (pz,z(1))
      pointer (py,y(1))
c
      common /statpg/ ir3(100),nr3(100),ir6(100),nr6(100),iw6(100)
     #,               nw6(100),ipgs(100),ipgd(100),ipgo(100),ipgij(100)
     #,               ipgout(100)
      common /io/     inp,iout
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
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      key=chrkey(ops,'int=drt=key','drt',' ')
      call pakstr(key,lenkey)
      call iosys('read character "drt file name '//key(1:lenkey)
     1            //'" from rwf',0,0,0,dsk)
      write(iout,*) '            reading drt information from '//dsk
      call iosys('read character "hamiltonian manipulation '//
     $           'filename" from rwf',0,0,0,nmcnfg)
      call iosys('open hconfig as old',0,0,0,nmcnfg)
      call getmem(0,p,ngot(1),'first',0)
      call iosys('read integer mncore from rwf',1,bottom,0,' ')
      call iosys('read integer mxcore from rwf',1,top,0,' ')
      write(iout,*) 'minimum core available = ',bottom
      write(iout,*) 'maximum core available = ',top
      if(chrkey(ops,'core-use=m901','minimum',' ').eq.'minimum') then
         top=bottom
      endif
      maxi=top
      maxr=iadtwp(maxi)
      file=tfile
      unit=tunit
c
      icount=0
      icntd=0
c
c
c     recover ci control parameters
      nroots=intkey(ops,'ci=nroots',1,' ')
      mxiter=intkey(ops,'ci=cycles',15,' ')
      ntol=intkey(ops,'ci=tolerance',5,' ')
      ithrsh=intkey(ops,'ci=threshhold',5,' ')
      nguess=intkey(ops,'ci=nguess',nroots,' ')
      nattim=intkey(ops,'ci=nroots-at-a-time',nroots,' ')
      gestyp=chrkey(ops,'ci=guess=','diagonals',' ')
c
      if(ntol.le.0) ntol=8
      cnverg=10.0d+00**(-ntol)
      thresh=10.0d+00**(-ithrsh)
      if(nroots.eq.0) nroots=1
      mxiter=mxiter*nattim+nguess
      ndvdmx=mxiter+iguess+1
      if(irooti.eq.0) irooti=1
      if(irootf.eq.0) irootf=irooti+nroots-1
c
      if(prtflg.ne.'minimum') then
         write(iout,17) nroots,nguess,mxiter,cnverg,thresh
      end if
c
c     ----- get all the dimensions etc. needed -----
c
      call iosys('read integer "number of drt functions" from '//dsk,
     $     1,nbf,0,' ')
      call iosys('read integer "symmetries in ci" from '//dsk,
     $     1,nsym,0,' ')
      call iosys('read integer norbs from '//dsk,1,norbs,0,' ')
      call iosys('read integer nrows from '//dsk,1,nrows,0,' ')
      nrows4=nrows*4
      call iosys('read integer nlevs from '//dsk,1,nlevs,0,' ')
      call iosys('read integer nrefs from '//dsk,1,nrefs,0,' ')
      call iosys('read integer nwks from '//dsk,1,nwks,0,' ')
      call iosys('read integer orbfrm from '//dsk,1,orbfrm,0,' ')
      call iosys('read integer symorb from '//dsk,1,symorb,0,' ')
      call iosys('read integer numij from '//dsk,1,numij,0,' ')
      call iosys('read integer ngroup from '//dsk,1,ngroup,0,' ')
      call iosys('read integer nmax from '//dsk,1,nmax,0,' ')
      call iosys('read integer maxb from '//dsk,1,bmax,0,' ')
      call iosys('read integer nijvir from '//dsk,1,nijvir,0,' ')
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
      nroots=min(nroots,nwks)
c
c
      nmax2=intowp(nmax)
      nwks2=intowp(nwks)
      lvfrm1=orbfrm
      levfrm=lvfrm1+1
c
c
c     if (ilvfrm.ne.0) levfrm=ilvfrm
c
c
      if (prtflg.ne.'minimum') then
         write(iout,32)nwks
         write(iout,33)nbf
         write(iout,34) norbs
      end if
c
c     calculate pointers into blank common
c
c     integer arrays located in ia 
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
c     real arrays located in za 
      acoef=iadtwp(hdpnt+nlevs)
      bcoef=acoef+nlevs
      root=bcoef+nlevs
      dvdvec=root+ndvdmx
      dvdmat=dvdvec+ndvdmx**2
      lowtri=ndvdmx*(ndvdmx+1)/2
      need=dvdmat+lowtri
      leni=wpadti(need)
      call getmem(leni,p,ngot(1),'getdrt',1)
      maxi=maxi-ngot(1)
      maxr=iadtwp(maxi)
c
      call izero(ia,ngot(1))
c
c
c
      call getdrt(ia(isegm),ia(aval),ia(bval)
     #,        ia(rowsym),ia(arc),ia(nlwks)
     #,        ia(weight),ia(wab),ia(wtw),ia(wtx),ia(wty)
     #,                   ia(kadd),ia(ladd),ia(ijadd),ia(ijgrp)
     #,        ia(inint),ia(inext),ia(jmnnxt),ia(jmxnxt),ia(ningrp)
     #,        ia(orbsym),ia(imngrp),ia(imxgrp),ia(jmngrp),ia(jmxgrp)
     #,        ia(ijxx),ia(klxx),ia(nklxx),ia(ijww),ia(klww),
     #         ia(nklww),dsk)
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
      ints=1
      ci=ints+nmax
      si=ci+lnsq
      cj=si+lnsq
      sj=cj+lnsq
      h1=sj+lnsq
      h2=h1+lnsq
      need=h2+lnsq
      leni=wpadti(need)
      if(need.gt.maxr) then
         write(iout,35)
         call lnkerr('quit')
      else
         call getmem(leni,pz,ngot(2),'getint',1)
         call rzero(z,iadtwp(ngot(2)))
         maxi=maxi-ngot(2)
         maxr=iadtwp(maxi)
      endif
c
c     figure out how many vectors can be held in remaing core.
c
      mxvc=(maxr-5000)/(2*nwks)
      mxvc=min(mxvc,64,nroots)
      sigma=1
      vector = sigma + mxvc*nwks
      need = vector + mxvc*nwks
      leni = wpadti(need)
      call getmem(leni,py,ngot(3),'vectors',1)
c
      if (mxvc.lt.1) then
          call lnkerr('not enough room to hold vectors')
      endif
      if (prtflg.ne.'minimum') then
          write (iout,36) mxvc,nroots
      end if
c
c
      call rzero(y,2*mxvc*nwks)
c
      call initpg(y(vector),y(sigma))
c
c..rlm     call initlp
      call initlp(z(ints),y(vector),y(sigma),
     $            ia(ijadd),ia(kadd),ia(ladd),
     $            ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym),
     $            z(ci),z(si),z(cj),z(sj),z(h1),z(h2),
     $            ia(ijxx),ia(nklxx),ia(ijww),ia(nklww),
     $            ia(klxx),ia(klww))
c
c..rlm     call initdg
      call initdg(z(ints),y(sigma),ia(ijadd),ia(kadd),ia(ladd),
     $            ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym))
c
c..rlm      call initex
      call initex(z(ints),y(vector),y(sigma),ia(ijadd),ia(kadd),
     $           ia(ladd),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym))
c
c..rlm     call initx2
      call initx2(z(ints),y(vector),y(sigma),ia(ijadd),ia(kadd),
     *            ia(ladd),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym))
c
c..rlm     call initst
      call initst(z(ints),y(vector),y(sigma),ia(kadd),
     $           ia(ladd),ia(wtw),ia(wtx),ia(wab),ia(orbsym),
     $           z(ci),z(si),z(cj),z(sj),z(h1),z(h2),
     $           ia(ijxx),ia(nklxx),ia(ijww),ia(nklww),
     $           ia(klxx),ia(klww))
c
cpws      call inimat(z(ints),ia(ijadd),ia(kadd),ia(ladd))
c
c..rlm     call init4x
      call init4x(z(ints),y(vector),y(sigma),ia(ijadd),ia(kadd),
     $           ia(ladd),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym))
c
      call start(za(dvdmat),z(ints),y(sigma),y(vector))
c
      call loopyd(ia(bval),ia(nlwks)                      ,ia(weight),
     $            ia(aval),ia(ijadd),ia(ijgrp),ia(kadd),ia(ladd),
     $            ia(inext),ia(arc),ia(jmnnxt),ia(jmxnxt),ia(rowsym),
     $            ia(orbsym)                    ,ia(isegm),
     $            ia(jsegm),ia(imain),ia(isub),ia(iuwkmn),ia(iuwksb),
     $            ia(itrack),ia(lmin), za(acoef), za(bcoef),
     $            ia(hdwgt),ia(hdseg),ia(hdpnt),ia(imxgrp),ia(imngrp),
     $            ia(jmxgrp),ia(jmngrp),z(ints),file,unit,y(vector),
     $            y(sigma),ia(wtw),ia(wtx),ia(wty),ia(wab),ia(orbsym),
     $            z(ci),z(si),z(cj),z(sj),z(h1),z(h2),
     $            ia(ijxx),ia(nklxx),ia(ijww),ia(nklww),
     $            ia(klxx),ia(klww))
c
c
c
c     initialize the diagonalization routine and generate the 
c     guess vectors.
      call fdiag(z(ints),y(vector),y(sigma),za(dvdmat),za(root),
     #           za(dvdvec),ndvdmx,thresh,nguess,prtflg,nattim,
     #           ia(bval),ia(nlwks),ia(weight)
     *,          ia(aval),ia(ijadd),ia(ijgrp),ia(kadd),ia(ladd)
     *,          ia(inext),ia(arc),ia(jmnnxt),ia(jmxnxt),ia(rowsym)
     *,          ia(orbsym),ia(isegm)
     *,          ia(jsegm),ia(imain),ia(isub),ia(iuwkmn),ia(iuwksb)
     *,          ia(itrack),ia(lmin), za(acoef), za(bcoef)
     *,          ia(hdwgt),ia(hdseg),ia(hdpnt),ia(imxgrp),ia(imngrp)
     #,          ia(jmxgrp),ia(jmngrp),file,unit,ia(wtw),ia(wtx),
     #           ia(wty),ia(wab),ia(orbsym)
     *,          z(ci),z(si),z(cj),z(sj),z(h1),z(h2)
     #,          ia(ijxx),ia(nklxx),ia(ijww),ia(nklww)
     #,          ia(klxx),ia(klww),gestyp)
c
c        iterate the diagonalization procedure.
      call roots(z(ints),y(vector),y(sigma),za(root),za(dvdvec),
     #           za(dvdmat),ndvdmx,prtflg,
     #           ia(bval),ia(nlwks),ia(weight)
     *,          ia(aval),ia(ijadd),ia(ijgrp),ia(kadd),ia(ladd)
     *,          ia(inext),ia(arc),ia(jmnnxt),ia(jmxnxt),ia(rowsym)
     *,          ia(orbsym),ia(isegm)
     *,          ia(jsegm),ia(imain),ia(isub),ia(iuwkmn),ia(iuwksb)
     *,          ia(itrack),ia(lmin), za(acoef), za(bcoef)
     *,          ia(hdwgt),ia(hdseg),ia(hdpnt),ia(imxgrp),ia(imngrp)
     #,          ia(jmxgrp),ia(jmngrp),file,unit,ia(wtw),ia(wtx),
     #           ia(wty),ia(wab),ia(orbsym)
     *,          z(ci),z(si),z(cj),z(sj),z(h1),z(h2)
     #,          ia(ijxx),ia(nklxx),ia(ijww),ia(nklww)
     #,          ia(klxx),ia(klww))
c
      call iosys('destroy bliu',0,0,0,' ')
c
c
      call getmem(-ngot(3),py,idum,'vectors',idum)
      call getmem(-ngot(2),pz,idum,'getint',idum)
      call getmem(-ngot(1),p,idum,'getdrt',idum)
      return
 16   format(10i5)
 17   format(1x,'m901:',
     $         /5x,'number roots sought:          ',i8,
     $         /5x,'number guess vectors:         ',i8,
     $         /5x,'maximum number iterations:    ',i8,
     $         /5x,'convergence criterion:        ',e8.1,
     $         /5x,'davidson vector threshhold:   ',e8.1)
 32   format(5x,'number of configurations:     ',i8)
 33   format(5x,'number of basis functions:    ',i8)
 34   format(5x,'number of allowed orbitals:   ',i8)
 35   format(/,5x,'not enough core for integral arrays')
 36   format(5x,'davidson vectors in use:      ',i8,
     $          /5x,'davidson vectors possible:    ',i8)
      end




