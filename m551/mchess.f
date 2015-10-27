*deck @(#)mchess.f	1.2  7/30/91
      subroutine mchess(xjbuf,nbufj,nbufk,xj,txk,xk,xkt,den,hess,hmo,
     $     dab,ldab,f,lfm,fab,lfab,grad,lgrad, noba,lnoba, c,isymc,
     $     aii,aij,bii,bij,bijt, nco,nao,nvo,nocc,nob,nbf, nsym, deg,
     $     temp,tv,  iscr, nfiv,nfav,nfhm,ndab,
     $     noci, itran, thrsh, mixinv,isymm,
     $     mix,leni,mixo,lok,locsym,iocore,itflag,nmix,ldahab,ldafab,
     $     linear, lhbuf,lhmo,icphf, fstrt,jstrt,kstrt,
     $     jnij,jmij,jnkl,jmkl, knij,kmij,knkl,kmkl, jkcore,
     $     bufix,bufax,lbufso,ilast,nbufcc,nbufac,incorh,
     $     hbuf,shmo,lshmo,lbhmo,asort,ncors,mblkd,lblkd)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      900417   (yymmdd)
c
c   17 april    1990   rlm at lanl
c     passing mblkd, lblkd to remove mclden entry point.
c
c   16 november 1987   bhl      lanl
c     modified to skip subroutines if nao=0
c     mjac kstsav jstsav  defined for closed shell scf runs
c     changes denoted by c..bhl..close shell scf
c   01 december 1987   bhl      brl
c     if (ilast.ne.0.and.iocore.ne.0) write mc_ao_ymatrix to rwf
c     for use by m1001 to solve cp-mcscf equations
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mchess.f	1.2   7/30/91
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c
      implicit real*8(a-h,o-z)
c
      dimension xjbuf(2), xj(2),xk(2),xkt(2),txk(2)
      dimension den(2),hess(2), dab(2),ldab(2), fab(2),lfab(2)
      dimension grad(2)
      dimension f(2), aii(2),aij(2),bii(2),bij(2),bijt(2)
      dimension nbf(2), hmo(2)
      dimension lgrad(2), noba(2),lnoba(2), nco(2),nao(2)
      dimension nvo(2),nocc(2),nob(2)
      dimension deg(2), temp(2), tv(2), isymc(2),c(2),lfm(2)
      dimension mix(2),leni(2),mixo(2),locsym(2),lok(2),mixinv(2)
      dimension isymm(2)
      integer mblkd(51,2), lblkd(2)
      real*8 bufix(lbufso),bufax(lbufso)
c
c
      dimension hbuf(*),shmo(*),lshmo(*),lbhmo(*),asort(*)
c
      common / number / zero,pt5,one,two,four,eight
c
      integer fstrt
c
      common /io/ inp,iout
c
c
c---------------------------------------------------------------------c
c     read the core fock matrix ( one-electron ints. if no core orbs)
c---------------------------------------------------------------------cc
c
      mmtot=lfm(nsym+1)-1
c
      call iosys('read real mc_core_fock from mcscr',mmtot,fab,0,' ')
c

      nbft=0
      noct=0
      do 100 ii=1,nsym
         nbft=nbft+nbf(ii)
         noct=noct+nocc(ii)
         if(nbf(ii).eq.0) go to 100
         call mcsqf(fab(lfm(ii)),f(lfab(ii)),nbf(ii))
 100  continue
c
      do 105 ii=1,nsym
         if(nbf(ii).eq.0) go to 105
         call mczero(fab(lfab(ii)),nbf(ii),nbf(ii))
         call mczero(grad(lgrad(ii)),nbf(ii),nocc(ii))
 105  continue
c
c-----------------------------------------------------
c     process the coulomb contributions to the hessian
c-----------------------------------------------------
c
      do 1000 mr=1,nsym
         do 1000 ms=1,mr
c
            naor=nao(mr)
            naos=nao(ms)
            ncor=nco(mr)
            ncos=nco(ms)
            noccr=nocc(mr)
            noccs=nocc(ms)
            nbfs=nbf(ms)
            nbfr=nbf(mr)
            mrs=nbfr*nbfs
c
c..bhl..closed shell scf
c
            mhac=1
            jstsav=jstrt
            kstsav=kstrt
            if(naor.eq.0.or.naos.eq.0)go to 40
c
            if(mr.eq.ms)go to 10
            nzero=naor*naos*mrs
            go to 15
 10         continue
            nzero=((naor*(naor+1))/2)*mrs
 15         continue
            do 20 i=1,nzero
               hess(i)=zero
 20         continue
c
            if(incorh.eq.0) then
               mhac=1
            else
               mhac=nzero+1
            endif
c
c
c
            call mcbigj(xj,hess,den,aii,aij,dab,ldab,deg,fab,lfab,
     $           noba,lnoba,nco,nao,nvo,nob,c,isymc,
     $           nsym,noci,temp,tv,itran,thrsh,mr,ms,ndab,nfiv,nfav,
     $           lhaa,lhac,lhca,lhcc, xjbuf,noct,nbft,nbufj,
     $           jstrt,jnij,jmij,jnkl,jmkl, jkcore,bufix,bufax,lbufso,
     $           lxjsav,ireadj,jstsav,irdj,nbufjj,mblkd,lblkd)
c
c
c
 40         continue
c
c
c----------------------------------------------------
c  add the fock operator contributions to the hessian
c----------------------------------------------------
c
            if(mr.ne.ms)go to 50
c
            lhaa=1
            call mcfock(dab(ldab(mr)),f(lfab(mr)),
     $           ncor,naor,mrs,hess,lhaa,lhcc,deg(mr))
c
c
 50         continue
c
c---------------------------------------------------
c    the coulomb contributions to hessian are used to
c    construct the gradient
c----------------------------------------------------
c
            call mcgrad(grad(lgrad(mr)),grad(lgrad(ms)),
     $           nco(mr),nco(ms),nao(mr),nao(ms),nob(mr),nob(ms),
     $           c(isymc(mr)),c(isymc(ms)), nbfr,nbfs, hess, mr,ms, tv,
     $           temp,xj)
c
c
c-----------------------------------------------------
c     add the exchange contributions to the hessian
c-----------------------------------------------------
c
c   array xk and txk are flip-flopped in the in-core and out-of-core
c   version of mcdnk.... xk(nbufjk)  txk(mrs)
c
c..bhl..closed shell scf
c
          if(naor.ne.0.and.naos.ne.0) then
c
               call mcbigk(xk,txk,hess,den,bii,bij,dab,ldab,deg,fab,
     $              lfab,noba,lnoba,nco,nao,nvo,nob,c,isymc,
     $              nsym,noci,temp,tv,itran,thrsh,xkt,bijt,mr,ms,
     $              ndab,nfiv,lhaa,lhac,lhca,lhcc, noct,nbft,nbufk,
     $              kstrt, knij,kmij,knkl,kmkl, jkcore,bufix,bufax,
     $              lbufso,lxksav,ireadk,kstsav,irdk,nbufkk,mblkd,lblkd)
c
          endif
c
c
c
 65         continue
c
            if(mr.ne.ms)go to 60
c
c-----------------------------------------------c
c      add fab operator to the hessian
c-----------------------------------------------c
c
c      call mcfab(fab(lfab(mr)),ncor,mrs,hess,lhcc)
c
c
 60         continue
c
            if(incorh.eq.0) then
c
               call iosys('write real mc_ao_ymatrix to rwf',
     $              nzero,hess,0,' ')
c
            endif
c
            call mccor(xj,xk,txk,hess(mhac),den,aii,aij,bii,bij,bijt,
     $           dab,ldab,deg,fab,lfab,noba,lnoba,nco,nao,nvo,nob,c,
     $           isymc,nsym,noci,temp,tv,itran,thrsh,mr,ms,ndab,nfiv,
     $           nfav,lhaa,lhac,lhca,lhcc,xjbuf,noct,nbf,
     $           nbufj,jstsav,jnij,jmij,jnkl,jmkl,jkcore,
     $           nbufk,kstsav,knij,kmij,knkl,kmkl,
     $           lxjsav,lxksav,ireadj,ireadk,nbufcc,nbufac,
     $           f,incorh,irdj,irdk,nbufjj,nbufkk)
c
c
 1000    continue
c
         ncor1=ncor+1
         ncos1=ncos+1
c
c------------------------------------------------
c  transform hessian to the mo basis
c------------------------------------------------
c
c
c=====
c     save ao y-matrix for cp-mcscf
c=====
      if(ilast.ne.0) then
      call iosys('write integer mc_incorh to rwf',1,incorh,0,' ')
       if(iocore.ne.0) then
        nthes=(nocc(1)*(nocc(1)+1)/2)*nbf(1)*nbf(1)
        call iosys('write real mc_ao_ymatrix to rwf',
     #  nthes,hess,0,' ')
       endif
      endif
c=====
c
         ms=1
         mr=1
c
c-3      write(iout,*)'  iocore ',iocore
c
         if(iocore.eq.0)go to 400
c
c
c
         call rzero(hmo,lhmo)
c
         ix=1
c
         ls1=isymm(ms)+1
         lr1=isymm(mr)+1
c
c
         if(naor.eq.0)go to 330
         ipr=isymc(mr)
         ips=isymc(ms)
         lr=locsym(mr)
         ls=locsym(ms)
         iv=lgrad(mr)
c
         do 320 is=ncos1,noccs
            locss=lok(ls+is)
            lenis=leni(ls+is)
            if(is.eq.ncos1)go to 315
            is1=is-1
            do 310 ir=ncos1,is1
               locrr=lok(lr+ir)
               lenir=leni(lr+ir)
c
c
               call mcthij(mixinv(iv),hess(ix),c(ips),c(ipr),nbfs,nbfr,
     $              lenis,lenir,
     $              mixo(locss+ls1),mixo(locrr+lr1),
     $              mix(locss+ls1),mix(locrr+lr1),
     $              locss,locrr,ms,is,mr,ir,
     $              tv,temp,hmo,itflag,nmix)
c
c
               ix=ix+mrs
 310        continue
 315        continue
c
            call mcthii(hess(ix),c(ips),nbfs,lenis,
     $           mixo(locss+ls1),mix(locss+ls1),
     $           locss,ms,is,
     $           tv,temp,hmo,itflag,nmix)
c
            ix=ix+mrs
 320     continue
c
 330     continue
c
         if(naor.eq.0.or.ncor.eq.0)go to 360
c
c
         ipr=isymc(mr)
         ips=isymc(ms)
         lr=locsym(mr)
         ls=locsym(ms)
         iv=lgrad(mr)
         do 350 is=1,ncos
            lenis=leni(ls+is)
            locss=lok(ls+is)
            do 340 ir=ncor1,noccr
               lenir=leni(lr+ir)
               locrr=lok(lr+ir)
c
c     mcthij --> mxthij
c
               call mxthij(mixinv(iv),hess(ix),c(ipr),c(ips),nbfr,nbfs,
     $              lenir,lenis,
     $              mixo(locrr+lr1),mixo(locss+ls1),
     $              mix(locrr+lr1),mix(locss+ls1),
     $              locrr,locss,mr,ir,ms,is,
     $              tv,temp,hmo,itflag,nmix)
c
c
               ix=ix+mrs
 340        continue
 350     continue
 360     continue
c
c
         if(ncor.eq.0.or.ncos.eq.0)go to 390
c
c
         ipr=isymc(mr)
         ips=isymc(ms)
         lr=locsym(mr)
         ls=locsym(ms)
         iv=lgrad(mr)
         do 380 is=1,ncos
            lenis=leni(ls+is)
            locss=lok(ls+is)
            if(is.eq.1)go to 375
            is1=is-1
            do 370 ir=1,is1
               lenir=leni(lr+ir)
               locrr=lok(lr+ir)
c
               call mcthij(mixinv(iv),hess(ix),c(ips),c(ipr),nbfs,nbfr,
     $              lenis,lenir,
     $              mixo(locss+ls1),mixo(locrr+lr1),
     $              mix(locss+ls1),mix(locrr+lr1),
     $              locss,locrr,ms,is,mr,ir,
     $              tv,temp,hmo,itflag,nmix)
c
               ix=ix+mrs
 370        continue
 375        continue
c
            call mcthii(hess(ix),c(ips),nbfs,lenis,
     $           mixo(locss+ls1),mix(locss+ls1),
     $           locss,ms,is,
     $           tv,temp,hmo,itflag,nmix)
c
            ix=ix+mrs
 380     continue
c
 390     continue
c
c
         go to 2000
 400     continue
c
c
c-3      write(iout,*)'  out-of-core construction of hessian '
c
c
c accume=-1   sort with accumulation
c
ct      call mysort('start',asort,asort,ncors,lhmo,
ct     1 0,0,0,0,'hessian','rwf',0)
c
         call sorter('start',asort,asort,ncors,lhmo,0,-1,0,0,
     $        'mcscf_hessian','rwf',.false.)
c
c
         ix=1
c
         naor=nao(mr)
         naos=nao(ms)
         ncor=nco(mr)
         ncos=nco(ms)
         noccr=nocc(mr)
         noccs=nocc(ms)
         nbfs=nbf(ms)
         nbfr=nbf(mr)
         mrs=nbfr*nbfs
         ls1=isymm(ms)+1
         lr1=isymm(mr)+1
         ncos1=ncos+1
         ncor1=ncor+1
c
c-2       write(iout,*)'  naor ncor ',naor,ncor
c
 1300    continue
c
         if(naor.eq.0)go to 1330
c=====
         ipr=isymc(mr)
         ips=isymc(ms)
         lr=locsym(mr)
         ls=locsym(ms)
         iv=lgrad(mr)
c=====
         call iosys('rewind mc_ao_ymatrix on rwf',0,0,0,' ')
c
         ijrd=(lhbuf/mrs)
         ijmx=(naor*(naor+1))/2
         ijrd=min(ijrd,ijmx)
         lread=ijrd*mrs
         call iosys('read real mc_ao_ymatrix from rwf '//
     $        'without rewinding',lread,hbuf,0,' ')
         ijleft=ijmx-ijrd
         ic=0
         ix=1
c=====
         do 1320 is=ncos1,noccs
            locss=lok(ls+is)
            lenis=leni(ls+is)
            if(is.eq.ncos1)go to 1315
            is1=is-1
            do 1310 ir=ncos1,is1
               locrr=lok(lr+ir)
               lenir=leni(lr+ir)
               lenir=leni(lr+ir)
c
               if(ic.lt.ijrd)go to 1311
               ix=1
               ic=0
               ijrd=min(ijleft,ijrd)
               lread=ijrd*mrs
               ijleft=ijleft-ijrd
               call iosys('read real mc_ao_ymatrix from rwf '//
     $              'without rewinding',lread,hbuf,0,' ')
 1311          continue
c
               ic=ic+1
c
c
               call scthij(mixinv(iv),hbuf(ix),c(ips),c(ipr),nbfs,nbfr,
     $              lenis,lenir,
     $              mixo(locss+ls1),mixo(locrr+lr1),
     $              mix(locss+ls1),mix(locrr+lr1),
     $              locss,locrr,ms,is,mr,ir,
     $              tv,temp,
     $              shmo,lshmo,lbhmo,asort,ncors,nmix)
               ix=ix+mrs
 1310       continue
 1315       continue
c
            if(ic.lt.ijrd)go to 1312
            ix=1
            ic=0
            ijrd=min(ijrd,ijleft)
            ijleft=ijleft-ijrd
            lread=ijrd*mrs
            call iosys('read real mc_ao_ymatrix from rwf '//
     $           'without rewinding',lread,hbuf,0,' ')
 1312       continue
c
            ic=ic+1
c
c
            call scthii(hbuf(ix),c(ips),nbfs,lenis,
     $           mixo(locss+ls1),mix(locss+ls1),
     $           locss,ms,is,
     $           tv,temp,
     $           shmo,lshmo,lbhmo,asort,ncors,nmix)
            ix=ix+mrs
 1320    continue
c
 1330    continue
c
         if(naor.eq.0.or.ncor.eq.0)go to 1360
c
         ix=1
         ic=0
         ijrd=(lhbuf/mrs)
         ijmx=naor*ncor
         ijrd=min(ijrd,ijmx)
         ijleft=ijmx-ijrd
         lread=ijrd*mrs
         call iosys('rewind hess_ac on rwf',0,0,0,' ')
         call iosys('read real hess_ac from rwf without rewinding',
     $        lread,hbuf,0,' ')
c
c
         ipr=isymc(mr)
         ips=isymc(ms)
         lr=locsym(mr)
         ls=locsym(ms)
         iv=lgrad(mr)
         do 1350 is=1,ncos
            lenis=leni(ls+is)
            locss=lok(ls+is)
            do 1340 ir=ncor1,noccr
               lenir=leni(lr+ir)
               locrr=lok(lr+ir)
c
               if(ic.lt.ijrd)go to 1341
               ix=1
               ic=0
               ijrd=min(ijrd,ijleft)
               ijleft=ijleft-ijrd
               lread=ijrd*mrs
               call iosys('read real hess_ac from rwf '//
     $              'without rewinding',lread,hbuf,0,' ')
 1341          continue
c
c
               ic=ic+1
               call sxthij(mixinv(iv),hbuf(ix),c(ipr),c(ips),nbfr,nbfs,
     $              lenir,lenis,
     $              mixo(locrr+lr1),mixo(locss+ls1),
     $              mix(locrr+lr1),mix(locss+ls1),
     $              locrr,locss,mr,ir,ms,is,
     $              tv,temp,
     $              shmo,lshmo,lbhmo,asort,ncors,nmix)
               ix=ix+mrs
 1340       continue
 1350    continue
 1360    continue
c
c
         if(ncor.eq.0.or.ncos.eq.0)go to 1390
c
         ix=1
         ic=0
         ijrd=(lhbuf/mrs)
         ijmx=ncor*(ncor+1)/2
         ijrd=min(ijrd,ijmx)
         ijleft=ijmx-ijrd
         lread=ijrd*mrs
         call iosys('rewind hess_cc on rwf',0,0,0,' ')
         call iosys('read real hess_cc from rwf without rewinding',
     $        lread,hbuf,0,' ')
c
c
         ipr=isymc(mr)
         ips=isymc(ms)
         lr=locsym(mr)
         ls=locsym(ms)
         iv=lgrad(mr)
         do 1380 is=1,ncos
            lenis=leni(ls+is)
            locss=lok(ls+is)
            if(is.eq.1)go to 1375
            is1=is-1
            do 1370 ir=1,is1
               lenir=leni(lr+ir)
               locrr=lok(lr+ir)
c
               if(ic.lt.ijrd)go to 1371
               ix=1
               ic=0
               ijrd=min(ijrd,ijleft)
               ijleft=ijleft-ijrd
               lread=ijrd*mrs
               call iosys('read real hess_cc from rwf '//
     $              'without rewinding',lread,hbuf,0,' ')
 1371          continue
c
c
               ic=ic+1
               call scthij(mixinv(iv),hbuf(ix),c(ips),c(ipr),nbfs,nbfr,
     $              lenis,lenir,
     $              mixo(locss+ls1),mixo(locrr+lr1),
     $              mix(locss+ls1),mix(locrr+lr1),
     $              locss,locrr,ms,is,mr,ir,
     $              tv,temp,
     $              shmo,lshmo,lbhmo,asort,ncors,nmix)
               ix=ix+mrs
 1370       continue
 1375       continue
c
            if(ic.lt.ijrd)go to 1376
            ix=1
            ic=0
            ijrd=min(ijrd,ijleft)
            ijleft=ijleft-ijrd
            lread=ijrd*mrs
            call iosys('read real hess_cc from rwf without rewinding',
     $           lread,hbuf,0,' ')
 1376       continue
c
c
            ic=ic+1
            call scthii(hbuf(ix),c(ips),nbfs,lenis,
     $           mixo(locss+ls1),mix(locss+ls1),
     $           locss,ms,is,
     $           tv,temp,
     $           shmo,lshmo,lbhmo,asort,ncors,nmix)
            ix=ix+mrs
 1380    continue
c
 1390    continue
c
 2000    continue
c
c
         return
         end
