*deck @(#)mciter.f	1.3  8/3/91
      subroutine mciter(nsym,nbf,nfo,nco,nao,nvo, nmix,
     $     mix,mixh,mixi,mixo,len,lok,locp,icas,
     $     nstate,istate,
     $     iocore,itflag,thrsh,itran,
     $     icicup,nleqps,nleqit,thrleq,
     $     iaugh,nheig,thre,nviter,innorm,
     $     nf14, nf16, nf30, nf36, nf39,
     $     nf46, nf49, nf81, nhd, nf82, nf83, nf84, lbuf84,
     $     nf91, nf92, nf93, lbufh, nf94,
     $     sqcdf, otout, ipfdm, igeom, opnscf, iprtg,smlld,shftd,
     $     cr,icr, nwint, nphes, ldar, icutah,nmaxah,thrcah,
     $     osqcdf,icphf,niter,jktrn,lsect,lbufso,ilast,itcorh,
     $     mcroot,ops)
c
c***begin prologue     mciter
c***date written       871022   (yymmdd)
c***revision date      910605   (yymmdd)
c
c    5 june    1991    rlm at lanl
c      passing mcroot.
c   17 april   1990    rlm at lanl
c      getting rid of the mclden entry point.
c   23 october 1987    bhl at brl
c      fixing  second-order scf by setting noci=1
c      before call to mcledr.
c   18 november 1987   bhl      lanl
c      skip over routines if nao=0
c      changed denoted by c..bhl..closed shell scf
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mciter.f	1.3   8/3/91
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       mciter
c
      implicit real*8(a-h,o-z)
c
c
      parameter (maxsym=20,mxsym1=21)
c
      character*(*) ops
      character*8 mcscf
      character*8 citype
      logical otout,debug
c
      integer nbf(nsym),nfo(nsym),nco(nsym),nao(nsym),nvo(nsym)
      integer mix(2),mixh(2),mixi(2),mixo(2),len(2),lok(2),locp(2)
      integer icas(nsym)
c
      integer nobt(maxsym)
      integer nocc(maxsym)
      integer isymc(maxsym+1),isymm(maxsym+1),locsym(maxsym+1)
      integer ldab(maxsym+1),lfab(maxsym+1)
      integer lfmm(maxsym+1),lgrad(maxsym+1)
      integer nobtt(maxsym),nob(maxsym)
      integer fstrt,opnscf,wpadti
      integer ihd(20)
c
      real*8 aii(maxsym),bii(maxsym),aij(200),bij(200),bijt(200)
      real*8 deg(10)
c
c
      common /io/ inp,iout
      common /number/ zero,pt5,one,two,four,eight
      integer icr(2)
      real*8 cr(1)
c
      data locsym/mxsym1*0/,isymc/mxsym1*0 /
      data isymm /mxsym1*0/,lfab/mxsym1*0/
      data ldab/mxsym1*0/,lfmm/mxsym1*0/
      data debug/.false./
c
c
      call iosys('read character mcscf_type from mcscr',-1,0,0,mcscf)
      if(opnscf.ne.0) mcscf='opnscf'
      call iosys('read character "mcscf: ci used" from rwf',
     $           -1,0,0,citype)
c
c------------------------------
c  initialize symmetry pointers
c------------------------------
c
c
      ls=0
      ic=1
      im=0
      ld=1
      lf=1
      lfm=1
      nocct=0
      naot = 0
c
      do 10 i=1,nsym
         naot = naot + nao(i)
c
         nocc(i)=nco(i)+nao(i)
         nob(i)=nocc(i)+nvo(i)
         nobtt(i)=nob(i)+nfo(i)
c
         nocct=nocct+nocc(i)
c
         ic=ic+nfo(i)*nbf(i)
         isymc(i)=ic
         ic=ic+nob(i)*nbf(i)
c
         isymm(i)=im
         lgrad(i)=im+1
         im=im+nbf(i)*nocc(i)
c
         locsym(i)=ls
         ls=ls+nocc(i)
c
         ldab(i)=ld
         ld=ld+nao(i)*nao(i)
c
         lfab(i)=lf
         lf=lf+nbf(i)*nbf(i)
c
         lfmm(i)=lfm
         lfm=lfm+(nbf(i)*(nbf(i)+1))/2
c
 10   continue
c
      isymc(nsym+1)=ic
      isymm(nsym+1)=im
      lgrad(nsym+1)=im+1
      locsym(nsym+1)=ls
      ldab(nsym+1)=ld
      lfab(nsym+1)=lf
      lfmm(nsym+1)=lfm
c
c---------------------------------------------------
c  call partial-transform program
c---------------------------------------------------
c
      ncco=nco(1)
      naao=nao(1)
      nobb=nob(1)
      nbbf=nbf(1)
      nbufjk=20000
c
      nnoc=ncco+naao
      nnb=(nnoc*(nnoc+1))/2
      lcore=nnb*(nbbf*nbbf+(nbbf*(nbbf+1))/2)+8*nbbf*nbbf+2*nbufjk+
     $      naao*(naao+1)/2
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(0,icr,maxint,'transformations',0)
      nwwp=iadtwp(maxint)
c
      if (lcore.lt.nwwp.and.jktrn.eq.0) then
         need=wpadti(lcore)
         call getscm(need,icr,ncore,'in-core transformation',0)
         nwwp=iadtwp(ncore)
         call intrn(cr,icr,nwwp,ncco,naao,nobb,nbbf,nbufjk,
     $              niter,ecore,enuc,ncsf,lsect,opnscf)
      else
         need=wpadti(maxscm)
         call getscm(need,icr,ncore,'out-of-core transformation',0)
         nwwp=iadtwp(ncore)
         call outtrn(cr,icr,nwwp,ncco,naao,nobb,nbbf,nbufjk,
     $               niter,ecore,enuc,ncsf,lsect,opnscf)
      end if
c
      if(opnscf.eq.0) then
         call iosys('read integer nwks from rwf',1,ncsf,0,' ')
      else
         ncsf=1
      endif
c
      call iosys('write integer mc_ncsf to rwf',1,ncsf,0,' ')
      call iosys('write integer mc_nmix to rwf',1,nmix,0,' ')
c
c
c
      lpqrs=(nsym*(nsym+1))/2
      lpqrs=(lpqrs*(lpqrs+1))/2
c
c---------------------------------
c     call the  direct-ci  package
c---------------------------------
c
      if(mcscf.ne.'opnscf') then
         ncore=ncore-10
         if (citype.eq.'m902') then
            call mn820(cr,cr,ncore,'mcscr','mcscr','mcscf','ci')
            call mn902(cr,cr,ncore,'ci',0,0,'mcscr','mcscr','mcscf',
     #                 mcroot,'mc')
         else if (citype.eq.'m903') then
            call mn903(cr,cr,ncore,'ci',0,0,'mcscr','mcscr','ci',
     #                 mcroot,'mc')
         else
            call lnkerr('unknown ci type: '//citype//' specified')
         end if
      else
         call iosys('write real "mc root 1" on rwf',
     $              3,cr,0,' ')
      endif
c
      if(ilast.ne.0.and.opnscf.eq.0) then
         nnao=naao*(naao+1)/2
         ntnao=nnao*nnao
         call iosys('read real "mo one-electron integrals" from mcscr',
     $               nnao,cr,0,' ')
         call iosys('write real "mo one-electron integrals" to rwf',
     $               nnao,cr,0,' ')
         call iosys('read real "mo two-electron integrals" from mcscr',
     $               ntnao,cr,0,' ')
         call iosys('write real "mo two-electron integrals" to rwf',
     $               ntnao,cr,0,' ')
      end if
c
c-------------------------------------
c     call the  density matrix package
c-------------------------------------
c
      if(mcscf.ne.'opnscf') then
         if (citype.eq.'m902') then
            call mn902(cr,cr,ncore,'density',cr,cr,'mcscr','mcscr',
     $                 'mcscf',mcroot,'mc')
            call mn830(cr,cr,ncore,'mcscr','mcscr','mcscf',
     $                 'guga density matrix','mo 1pdm','mo 2pdm')
         else if (citype.eq.'m903') then
            call mn903(cr,cr,ncore,'density',cr,cr,'mcscr','mcscr',' ',
     $                 mcroot,'mc')
         end if
      else if (naao.ne.0) then
c        closed shell scf
         call scfden(cr,naao,naao*(naao+1)/2)
      endif
c
c--------------------------------------------
c     square-up the 2-electron density matrix
c--------------------------------------------
c
      nnsym=(nsym*(nsym+1))/2
      nnsym=(nnsym*(nnsym+1))/2
cc
c      temporary storage
cc
c  ----- nblock = 1  for nsym = 1  on cray   bhl 8/14/85
cc
      nblock = 1
      numint=naot*(naot+1)/2
      numint=numint*(numint+1)/2
cc
      idab = 1
      ipare=wpadti(idab    + naot*naot   + 1)
      ipqrs= ipare   + naot+1
      imblk= ipqrs   + nnsym + 1
      ibuf = iadtwp(imblk   + 51*nblock + 1)
      iden = ibuf    + numint  + 1
      need = wpadti(iden)
c
      call getscm(need,icr,ncore,'mcsqdm call',0)
c
      nwwp =iadtwp(ncore - wpadti(iden))
      lastd = 1
c
c
      call iosys('write integer mc_nbasis to rwf',1,nbbf,0,' ')
      call iosys('write integer mc_ncore to rwf',1,ncco,0,' ')
      call iosys('write integer mc_nactive to rwf',1,naao,0,' ')
      call iosys('write integer mc_norbs to rwf',1,nobb,0,' ')
      call iosys('write integer mc_lbufso to rwf',1,lbufso,0,' ')
      call iosys('write integer mc_lbufh to rwf',1,lbufh,0,' ')
c
c..bhl..closed shell scf
c
      if (naao.gt.0) then
         call mcsqdm(nsym,naot,nao,ncco,icr(ipare), nf36,cr(idab),
     $        cr(iden),cr(ibuf),nwwp,icr(imblk),icr(ipqrs), nblock,
     $        lastd,linear,ipfdm,niter,ecore,enuc,ops)
      else
         call iosys('write real energy to rwf',1,ecore+enuc,0,' ')
         call iosys('write real "mcscf energy" to rwf',
     $              1,ecore+enuc,0,' ')
      end if

c
c
c-----------------------------------------------
c     construct the orbital hessian and gradient
c-----------------------------------------------
c
c     core allocation for hessian construction
c     header record from transformed integral file
      call iosys('read integer mc_jk_hdr from mcscr',20,ihd,0,' ')
c
      jstrt=ihd(16)
      kstrt=ihd(17)
      fstrt=ihd(18)
c
c
      maxmrs = 0
      maxhes = 0
      maxnbf = 0
c
      nocci=nocc(1)
      ncoi=nco(1)
      naoi=nao(1)
      nbfi=nbf(1)
      mrs = nbfi*nbfi
      nocct=nocci*(nocci+1)/2
c
      maxmrs=mrs
      maxnbf=nbfi
      maxhes=nocct*mrs
c
      lahes=(naoi *(naoi+1)/2)*mrs
      lches=(naoi+1)*mrs
      minhes=max(lahes,lches)
c
      nwwp  =  nwwp - lastd
c
      lhmo= nmix * nmix
c
      icmo  =  iden  +  lastd + 1
      inoba = wpadti(icmo  +  ic)
      lnoba =  inoba +  nocct + 1
cc
cc
      ifc   =  iadtwp(lnoba + nsym + 1)
      ifab  =  ifc   +  lf
      ibufax=  ifab  +  lf
      ibufix= ibufax + lbufso
      igrad = ibufix+lbufso
c
      noc2=nocci*(nocci+1)/2
      nbf2=nbfi*(nbfi+1)/2
      nbk=nbufjk/maxmrs
      if(nbk.gt.noc2) nbk=noc2
      nbufkj=nbk*maxmrs
      if(nbufjk.lt.maxmrs) nbufjk=maxmrs
      nbufk=nbufjk
      nbj=nbufjk/nbf2
      if(nbj.gt.noc2) nbj=noc2
      nbufj=nbj*nbf2
c
      ixk   =  igrad +  isymm(nsym+1)
      ixkt = ixk + maxmrs
      ixj =  ixkt + maxmrs
      itemp =  ixj  +  maxmrs
      itv   =  itemp +  maxmrs
      iczmo =  itv   +  maxnbf
      ixjbuf =  iczmo +  ic
      itxk = ixjbuf + nbufjk
      ihess = itxk + nbufjk
c

c
      call getscm(0,cr,maxscm,'siz of core?',0)
c
      maxscm=iadtwp(maxscm)
      iwork = maxscm- ihess
      ineed = ihess + minhes
c
      if(debug) then
         write(iout,*)' mciter:iwork  maxscm ',iwork,maxscm
         if(ineed.gt.maxscm) then
            write(iout,*)' ineed maxscm ',ineed,maxscm
            call lnkerr(' increase memory for the hessian construction')
         endif
      endif
c
      if(iwork.gt.maxhes+lhmo.and.itcorh.eq.1) then
         incorh=1
         iocore=1
         icc=ncoi*(ncoi+1)/2
         iac=ncoi*naoi
         ineed=ihess+maxhes+lhmo
         ihmo=ihess+maxhes
         ineedh=ihmo+lhmo
      else
         incorh=0
         iocore=0
         call fixbuf(iwork,maxmrs,naoi,ncoi,icc,iac)
         ineed=maxscm-100
         ihmo=ixjbuf
         ihbuf=ihmo+nmix
         ishmo=ihbuf+lbufh
         ilshmo=ishmo+2*maxmrs
         ilbhmo=ilshmo+2*maxmrs
         isort=ilbhmo+2*maxmrs
         ncors=ineed-isort
         if(ncors.lt.1) then
            call lnkerr(' increase the core for hessian sort')
         end if
      endif
c
      if(debug) then
         write(iout,*)' mciter:incorh icc iac ',incorh,icc,iac
      endif
c
      if (ineed.gt.maxscm) then
         write(iout,*)' ineed maxscm ',ineed,maxscm
         call lnkerr(' increase core for hessian construction')
      endif
c
c
      call getscm(wpadti(ineed),cr,ncor,'iterations',0)
c
c
c---------------------------------------c
c   read orbitals from rwf
c---------------------------------------c
c
      nobbt=nob(1)
      nbfft=nbf(1)
      ntot=nobbt*nbfft
c
      call iosys('read real mo_old from mcscr',ntot,cr(icmo),0,' ')
c
c
      noci = 0
      if (icicup .eq. 0) noci = 1
c
      if(debug) then
         write(iout,*)' mciter:before mchsdr   nphes ',nphes
      endif
c
      call mchsdr( cr(ixjbuf), nbufj,nbufk, cr(ixj),cr(itxk),cr(ixk),
     $     cr(ixkt),cr(iden),cr(ihess),cr(ihmo),cr(idab),ldab,cr(ifc),
     $     lfmm,cr(ifab),lfab,cr(igrad),lgrad,
     $     icr(inoba),icr(lnoba),cr(icmo),isymc,
     $     aii,aij,bii,bij,bijt,
     $     nco,nao,nvo,nocc,nob,nbf, nsym,deg, cr(itemp),cr(itv),
     $     iscr,nf91,nf92,nf93,nf16, noci, itran, thrsh,
     $     mix,len,mixh,mixi,mixo,lok,locsym,locp,
     $     iocore,itflag,icas,nmix,lbufh,isymm,linear,nphes,
     $     cr(iczmo), nf46, igeom, iprtg, lhbuf, lhmo,
     $     icphf, nblock, nhd,ldar, fstrt,jstrt,kstrt,
     $     jnij,jmij,jnkl,jmkl, knij,kmij,knkl,kmkl, jkcore,
     $     cr(ibufax),cr(ibufix),lbufso,ilast,icc,iac,incorh,
     $     cr(ihbuf),cr(ishmo),cr(ilshmo),cr(ilbhmo),cr(isort),ncors,
     $     icr(imblk),icr(ipqrs))
c
      if(ilast.ne.0) then
c
         nocct=ncco+naao
         nocnob=nocct*nobb
         call iosys('write integer mc_locsym to rwf',
     $        nsym+1,locsym,0,' ')
         call iosys('write integer mc_mixo to rwf',
     $        nocnob,mixo,0,' ')
         call iosys('write integer mc_locp to rwf',
     $        nocct,locp,0,' ')
         call iosys('write integer mc_len to rwf',
     $        nocct,len,0,' ')
c
         return
c
      end if
c
      if (iaugh .eq. 0) then
c------------------------------
c        solve linear equations
c------------------------------
         nao2=(nao(1)*(nao(1)+1))/2
         nao22=nao2*nao2
         ldel = iadtwp(inoba)
         lg1 = ldel + nmix + ncsf
         lg2 = lg1 + nao2
         lsg = lg2+nao22*2
         lbufix = lsg+nao22*2
         lstor=lbufix+lbufso
         ncor = iadtwp(ncore) - lstor
c
         write(iout,90901)
90901    format('  second order option invoked .. ')
c
         if (opnscf.ne.0) noci=1
cbhl
c
         call mcledr(nsym,nbf,nob,nfo,nco,nao,cr(icmo),cr(lg1),
     $               cr(lg2),locsym,len,locp,mixo,nmix,
     $               nf14,nf49,nf81,nf82,nhd,nf16,ldafab,nf93,
     $               lbufh,nf92,nf91,nf94,nf83,nf84,lbuf84,mdstv,
     $               noci,nleqps,nleqit,thrleq,
     $               cr(ldel),sqcdf,cr(lstor),icr(wpadti(lstor)),
     $               ncor,ncsf,cr(lbufix),lbufso,mcroot)
c
cbhl
         osqcdf=sqcdf
      else 
c-------------------------------------
c        diagonalize augmented hessian
c-------------------------------------
         ldel = iadtwp(inoba)
         lstor = ldel + nmix
         ncor = iadtwp(ncore) - lstor
         call mcaugh(nmix,nf93,lbufh,nheig,thre,nviter,innorm,
     $               cr(ldel),sqcdf,smlld,shftd, cr(lstor),
     $               icr(wpadti(lstor)), ncor,osqcdf,icutah,
     $               nmaxah,thrcah)
         osqcdf=sqcdf
      endif
c
c----------------------------------
c        form new mcscf orbitals
c----------------------------------
c
      id2h = 0
cc
      do 901 iv=1,nsym
         nobt(iv)= nob(iv)+nfo(iv)
c        write(iout,902) nob(iv),nobt(iv),nbf(iv),nocc(iv)
c 902   format(' nob nobt nbf nocc ',5i8)
 901  continue
c
      call mcvec(nsym,nbf,nob,nocc,cr(icmo),locsym,locp,len,
     $           mixo,cr(ldel),id2h,cr(lstor),ncor,nfo)
c
c----------------------------------
c     output mcscf orbitals
c----------------------------------
c
      ntot=nbfft*nobbt
      call iosys('write real mo_new to mcscr',ntot,cr(icmo),0,' ')
c
      return
      end
