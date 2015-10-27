*deck @(#)fdiag.f	5.1  11/6/94
      subroutine fdiag(ints,v,d,dvdmat,root,dvdvec,ndvdmx,thresh,nguess,
     #                 prtflg,nattim,
     #                 nabc,nlwks,iwght,nabca,ijadd
     #,                ijgrp,kadd,ladd,inext,iarc,jmnnxt,jmxnxt
     #,                nabcs,isym,isegm,jsegm,imain
     #,                isub,iuwkmn,iuwksb,itrack,lmin,acoef,bcoef
     #,                hdwgt,hdseg,hdpnt,imxgrp,imngrp,jmxgrp,jmngrp
     #,                file,unit,wtw,wtx,wty,wab,ss
     #,                ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                klxx,klww,gestyp)
c
c***begin prologue     fdiag.f
c***date written       850601
c***revision date      11/6/94      
c
c      7 march   1991  rlm at lanl
c        fixing calls to bliu - in particular the mixed mode
c        problem with initialize.
c     15 January 1991  rlm at lanl
c        adding option to read guess ci vectors from the chk.
c     14 august 1987   pws at lanl
c        adding nroots-at-a-time option (nattim)
c
c***keywords           ci, guess
c***author             saxe,paul (lanl) 
c***source             @(#)fdiag.f	5.1   11/6/94
c***purpose            prepares a guess for the ci wavefunction. 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       fdiag.f
      implicit real*8 (a-h,o-z)
c
      real*8 ci(1),cj(1),si(1),sj(1),h1(1),h2(1)
      real*8 d(nwksmx,mxvc),v(nwksmx,mxvc),ints(nmax),dvdmat(lowtri)
      real*8 root(*),dvdvec(*),thresh
      real*8 acoef(nlevs),bcoef(nlevs)
      real*8 sdot,norm,ovrlap
      integer arr,refwlk,symorb,bmax,orbfrm
      integer ijxx(numij ),nklxx(nsym,orbfrm),ijww(numij )
      integer nklww(nsym,orbfrm),klxx(1),klww(1)
      integer wtw(orbfrm,nsym),wtx(orbfrm,nsym),wty(orbfrm),
     #        wab(orbfrm),ss(norbs)
      integer old,guess
      character*8 prtflg,status
      character*16 file,unit
      character*16 cjunk
      character*128 namchk
      character*4 itoc
      character*(*) gestyp
c
c     level characteristics
      integer isegm(nlevs),jsegm(nlevs),imain(nlevs),isub(nlevs)
      integer iuwkmn(nlevs),iuwksb(nlevs),itrack(nlevs),isym(norbs)
      integer lmin(nlevs)
      integer hdseg(nlevs),hdpnt(nlevs),hdwgt(nlevs)
c     integer levpt(nlevs),levnr(nlevs)
c
c     graph description arrays
      integer nabc(nrows ),nlwks(nrows)              ,nabcs(nrows )
      integer              iwght(nrows4),iarc(nrows4),nabca(nrows )
c
c     integral addressing arrays
      integer inext(norbs),jmnnxt(norbs),jmxnxt(norbs)
      integer ijadd(numij),ijgrp(numij),kadd(symorb),ladd(symorb)
      integer imngrp(ngroup),imxgrp(ngroup),jmngrp(ngroup)
      integer jmxgrp(ngroup)
c
c
      logical pagein
c
c     common blocks
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /optns/  iguess,irstrt,irooti,irootf,ikl34
      common /symq/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8),ibl(8)
      common /d901/ rep,fzcore,eguess,eci,refwlk,mxiter,cnverg,icnvg
     *,             iter,sqcdif,czero,nroot
      common /io/   inp,iout
      common /all/  val1,val2,val3,arr,itr1,itr2,ia,ja,itype,isegt
     *,             lvfrm1,nlwki,nlwkj,imax,imin
      common /page/ iword3,iword4,pagein,noffi,noffj,npass,nwordl,iword6
      common /rstrt/  iblck1,inxt,lowtri,ndvdit,iblock
      common /nvectr/ nvc,mxvc
c
 1001 format(5x,'nuclear repulsion:    ',g20.12)
 1002 format(5x,'reference energies:   ',g20.12,'(',i8,')')
 1003 format(27x,g20.12,'(',i8,')')
 1004 format(5x,'# guess vectors read from chk',i9)
c
      if (prtflg.eq.'minimum') then
         status='noprint'
      else
         status='print'
      end if
c
c     initialize the diagonalization routine.
      sqcdif=1.0d+00
      icnvg=-1
      call bliu('initialize',status,d,thresh,nwks,mxiter,nroots,iout,
     #         nattim,rep+fzcore,cnverg,dvdvec)
c
c     generate initial guesses.
      ndvdit=3
      nvc=1
      if (gestyp.eq.'diagonals') then
c        prepare for unit guess on lowest diagonal matrix element
         call iosys('write real tempd on bliu',nwks,d,0,' ')
         do 60 guess=1,nguess
            call iosys('read real tempd from bliu',-1,d,0,' ')
            eguess=0.0d+00
            do 50 j=1,nwks
               if(d(j,1).gt.eguess) go to 40
               refwlk=j
               eguess=d(j,1)
   40          continue
   50       continue
            d(refwlk,1)=0.0
            call iosys('write real tempd to bliu',nwks,d,0,' ')
c
c
            call rzero(v,nwks)
            v(refwlk,1)=1.0d+00
            if (prtflg.ne.'minimum') then
               if(guess.eq.1) then
                  write(iout,1001) rep
                  write(iout,1002) eguess+rep+fzcore,refwlk
               else
                  write(iout,1003) eguess+rep+fzcore,refwlk
               endif
            end if
c
            call rzero(d,nwks)
            call loopy(nabc,nlwks,iwght,nabca,ijadd
     *,                 ijgrp,kadd,ladd,inext,iarc,jmnnxt,jmxnxt
     *,                 nabcs,isym,isegm,jsegm,imain
     *,                 isub,iuwkmn,iuwksb,itrack,lmin,acoef,bcoef
     *,                 hdwgt,hdseg,hdpnt,imxgrp,imngrp,jmxgrp,jmngrp
     *,                 ints,file,unit,v,d,wtw,wtx,wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
            call bliu('with vectors',0,v,d,nwks,mxiter,0,0,0,0,0,0)
   60       continue
      else if (gestyp.eq.'references') then
      else if (gestyp.eq.'chk') then
c
c        read vectors from a previous calculation.
c        open the checkpoint file
         call iosys('read character "checkpoint filename" from rwf',
     $              0,0,0,namchk)
         call iosys('open chk as old',0,0,0,namchk)
c
         do 80 guess=1,nroots
            call rzero(v(1,guess),nwks)
            call iosys('read real "ci root '//itoc(guess)
     $                 //'" from chk',-1,v(1,guess),0,' ')
c
c           orthonormalize this root to previous ones.
            do 70 old=1,guess
               ovrlap=sdot(n,v(1,guess),1,v(1,old),1)
               call saxpy(n,-ovrlap,v(1,old),1,v(1,guess),1)
               norm=sqrt(sdot(n,v(1,guess),1,v(1,guess),1))
               call sscal(n,1.0d+00/norm,v(1,guess),1)
   70       continue
c
            if (prtflg.ne.'minimum'.and.guess.eq.1) then
               write(iout,1001) rep
               write(iout,1004) nroots
            end if
c
            call rzero(d,nwks)
            call loopy(nabc,nlwks,iwght,nabca,ijadd
     *,                ijgrp,kadd,ladd,inext,iarc,jmnnxt,jmxnxt
     *,                nabcs,isym,isegm,jsegm,imain
     *,                isub,iuwkmn,iuwksb,itrack,lmin,acoef,bcoef
     *,                hdwgt,hdseg,hdpnt,imxgrp,imngrp,jmxgrp,jmngrp
     *,                ints,file,unit,v(1,guess),d,wtw,wtx,wty,wab,ss
     #,                ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                klxx,klww)
            call bliu('with vectors',0,v(1,guess),d,nwks,mxiter,
     $                 0,0,0,0,0,0)
   80    continue
         call iosys('close chk',namchk,0,0,' ')
      else
         cjunk=gestyp
         call lnkerr('unrecognized guess type for ci wavefunction'
     $               //cjunk)
      end if
c
c
      return
      end
