*deck @(#)roots.f	5.1  11/6/94
      subroutine roots(ints,v,s,root,dvdvec,dvdmat,ndvdmx,prtflg,
     #                  nabc,nlwks,iwght,nabca,ijadd
     *,                 ijgrp,kadd,ladd,inext,iarc,jmnnxt,jmxnxt
     *,                 nabcs,isym,isegm,jsegm,imain
     *,                 isub,iuwkmn,iuwksb,itrack,lmin,acoef,bcoef
     *,                 hdwgt,hdseg,hdpnt,imxgrp,imngrp,jmxgrp,jmngrp
     *,                 file,unit,wtw,wtx,wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
c
      implicit integer (a-z)
c
      integer symorb
      real*8 ints(*),v(nwks,mxvc),s(nwks,mxvc),root(nroots)
      real*8 dvdvec(*),dvdmat(*)
      real*8 acoef(nlevs),bcoef(nlevs)
      real*8 rep,fzcore,eguess,eci,cnverg,sqcdif,czero,edav
      character*20 status,itoc*4
      character*8 prtflg
      character*16 file,unit
      real*8 ci(1),cj(1),si(1),sj(1),h1(1),h2(1)
      integer ijxx(numij ),nklxx(nsym,orbfrm),ijww(numij )
      integer nklww(nsym,orbfrm),klxx(1),klww(1)
c
c     level characteristics
      dimension isegm(nlevs),jsegm(nlevs),imain(nlevs),isub(nlevs)
      dimension iuwkmn(nlevs),iuwksb(nlevs),itrack(nlevs),isym(norbs)
      dimension lmin(nlevs)
      dimension hdseg(nlevs),hdpnt(nlevs),hdwgt(nlevs)
c
c     graph description arrays
      dimension nabc(nrows ),nlwks(nrows)              ,nabcs(nrows )
      dimension              iwght(nrows4),iarc(nrows4),nabca(nrows )
c
c     integral addressing arrays
      dimension inext(norbs),jmnnxt(norbs),jmxnxt(norbs)
      dimension ijadd(numij),ijgrp(numij),kadd(symorb),ladd(symorb)
      dimension imngrp(ngroup),imxgrp(ngroup),jmngrp(ngroup)
      dimension jmxgrp(ngroup)
      integer wtw(orbfrm,nsym),wtx(orbfrm,nsym),wty(orbfrm),
     #        wab(orbfrm),ss(norbs)
      logical debug
c
      common /io/     inp,iout
      common /dims/   nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs,
     #                nrowoc,levfrm,nwksmx,nlwkmx,nuwkmx,bmax,nroots,
     #                orbfrm
      common /d901/   rep,fzcore,eguess,eci,refwlk,mxiter,cnverg,icnvg,
     #                iter,sqcdif,czero,nroot
      common /nvectr/ nvc,mxvc
c
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
c
c
      parameter (debug=.false.)
c
c     iteration loop for ci solution
   10 continue
         call bliu('solve',status,v,s,nwks,mxiter,0,0,0,dvdmat,root,
     #             dvdvec)
         if (status.eq.'converged') go to 100
c
         nvc=0
   20    continue
            nvc=nvc+1
            call bliu('new trial',status,v(1,nvc),s(1,nvc),nwks,
     #                mxiter,0,0,0,0,root,dvdvec)
            if (status.ne.'done'.and.nvc.lt.mxvc) go to 20
            if (status.eq.'done'.and.nvc.eq.1) go to 10
            if (status.eq.'done') nvc=nvc-1
c
            call rzero(s,nwks*nvc)
            call loopy(nabc,nlwks,iwght,nabca,ijadd
     *,                 ijgrp,kadd,ladd,inext,iarc,jmnnxt,jmxnxt
     *,                 nabcs,isym,isegm,jsegm,imain
     *,                 isub,iuwkmn,iuwksb,itrack,lmin,acoef,bcoef
     *,                 hdwgt,hdseg,hdpnt,imxgrp,imngrp,jmxgrp,jmngrp
     *,                 ints,file,unit,v,s,wtw,wtx,
     *                  wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
            do 30 ivc=1,nvc
               if (debug) then
                  write (iout,83) ivc,(v(iq,ivc),iq=1,nwks)
   83             format ('1'//,' vector ivc=',i5,/,(1x,10f12.6))
                  write (iout,84) ivc,(s(iq,ivc),iq=1,nwks)
   84             format (//,' sigma, ivc=',i5,/,(1x,10f12.6))
               end if
               call bliu('with vectors',0,v(1,ivc),s(1,ivc),nwks,
     #                    mxiter,0,0,0,0,0,0)
   30       continue
            if (status.eq.'done') go to 10
            nvc=0
         go to 20
c
  100 continue
c
c     ----- recover each root's ci vector, find the most important
c           configuration, and write vector to rwf
c
      call iosys('read real diagonals from bliu',-1,s,0,' ')
c
      if (prtflg.ne.'minimum') then
c        write(iout,*)'  roots:  nuclear repulsion energy ',rep
         write (iout,90)
   90    format (///,' root   reference     energy         ci energy ',
     #           '  davidson energy        c(0)')
      end if
c
      do 200 iroot=1,nroots
         call bliu('get vector',status,v,0,nwks,0,iroot,0,0,0,0,0)
         if (status.ne.'ok') go to 200
         call iosys('write real "ci root '//itoc(iroot)//'" to rwf',
     #               nwks,v,0,' ')
c
         czero=0.0d+00
         refwlk=0
         do 110 i=1,nwks
            if (abs(v(i,1)).gt.abs(czero)) then
               czero=v(i,1)
               refwlk=i
            end if
  110    continue
c
         eguess=s(refwlk,1)+rep+fzcore
         eci=root(iroot)+rep+fzcore
         edav=eci+(eci-eguess)*(1.0d+00-czero**2)
         if (iroot.eq.1) then
            call iosys('write real energy to rwf',1,eci,0,' ')
         end if
c
         call iosys('write real "ci energy '//itoc(iroot)//'" to rwf',
     #              1,eci,0,' ')
c
c
         if (prtflg.ne.'minimum') then
            write (iout,120) iroot,refwlk,eguess,eci,edav,czero
  120       format (1x,i3,i10,3g18.9,f8.4)
         end if
c
  200 continue
c
      call bliu('finish',0,0,0,0,0, 0,0,0,0,0,0)
c
c
      return
      end
