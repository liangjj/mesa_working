      subroutine mcledr(nsym,nbf,nob,nfob,ncob,naob,cv,g1,g2,
     $     locsym,len,locp,mix,nmix,
     $     nf16,nfsm,ndshd1,ndshd2,nhd,nda1,lda1,nf41,lenb,
     $     nf35,nf36,nf37,nf45,nf46,lb46,ndstv,
     $     noci,npassj,istopj,stopj,
     $     b, sqcdf,
     $     cr,icr,ncore,ncsf,bufix,lbufso,mcroot)
c
c***begin prologue     mcledr
c***date written       871022   (yymmdd)
c***revision date      900417   (yymmdd)
c
c   17 april    1990   rlm at lanl
c    removing need for the enmxmc entry point.
c   17 november 1987   bhl at lanl
c   skip over ci-gradient routine if noci > 0
c   changes denoted by c..bhl..scf
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcledr.f	1.3  7/30/91
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       mcledr
c
      implicit real*8 (a-h,o-z)
      integer wpadti
      character*4 itoc
      character*8 citype
      character*8 mcscf,prtflg
      real*8 bufix(lbufso)
c
      dimension nbf(2),nob(2),ncob(2),naob(2),cv(2),g1(2),g2(2)
      dimension locsym(2),len(2),locp(2),mix(2),nfob(2)
      dimension cr(2),icr(2)
      dimension b(2)
c
      common /io/ inp,iout
c
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
c
      call iosys('read character mcscf_type from mcscr',-1,0,0,mcscf)
c..bhl..scf
      if(noci.eq.0) then
         call iosys('read character "mcscf: ci used" from rwf',
     $              -1,0,0,citype)
      else
         ncsf=0
      endif
c..bhl..scf
c
c----------------------------------------
c     storage allocation for jacobi
c----------------------------------------
c
      call getscm(0,cr,maxcor,'size',0)
      maxcor=iadtwp(maxcor-100000)
c
      mdim = nmix + ncsf
      lgrad = 1
      ldiag = lgrad + mdim
      lt = ldiag + mdim
      lc = lt + mdim*(istopj+1)
      lc0 = lc + mdim
      lhess = lc0 + ncsf
c
      nhess = istopj * (istopj + 1) / 2
c
      lthess = lhess + nhess
      lgl = lthess + nhess
      lolap = lgl + istopj
      lsg = lolap + istopj
c
      naa=naob(1)*(naob(1)+1)/2
      nai=naob(1)*nob(1)
      nci=ncob(1)*nob(1)
      naa2=naa*naob(1)*naob(1)
c
      ltg=lsg+naa2
      labcx=ltg+naa2
      labix=labcx+naa*nai
      laibx=labix+naa*nci
      lstor=laibx+1
c
      incor=1
      ihss=lstor
      if(lstor.gt.maxcor) then
         lstor=labcx
         labix=lstor
         laibx=lstor
         incor=0
      else
         lhss=nmix*nmix
         if(lstor+lhss+lenb.lt.maxcor) then
            ihss=lstor
            ibf=ihss+lhss
            lstor=ibf+lenb
            incor=2
         end if
      end if
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(wpadti(lstor),cr,ncore,'mc linear equations',0)
c
      ncor = iadtwp(ncore) - lstor
c
      if(incor.eq.2) then
         call gethes(cr(ihss),lhss,cr(ibf),lenb,nf41)
      end if
c
c..bhl..scf..start
      if (naob(1).ne.0) then
         if(incor.ne.0) then
            call getrab(cr(labcx),bufix,lbufso,'abcx',
     $           naob(1),nob(1),naob(1))
cbhl
c fix for ncob = 0
cbhl
            if(ncob(1).ne.0) then
               call getrab(cr(labix),bufix,lbufso,'abix',
     $              ncob(1),nob(1),naob(1))
            endif
c
c      call addrab(cr(labix),bufix,lbufso,'aibx',
c    #   ncob(1),nob(1),naob(1))
c
         endif
      end if
c..bhl..scf.start
      if(noci.eq.0) then
cc
c--------------------------------
c     read ci vector
c--------------------------------
cc
      call iosys('read real "mc root '//itoc(mcroot)//'" from rwf',
     $            ncsf,cr(lc0),0,' ')
c
c---------------------------------
c     calculate ci gradient
c---------------------------------
c
      junk=0
      if (citype.eq.'m902') then
         call mn820(cr(lstor),cr(lstor),ncor,'mcscr','mcscr',
     $              'mcscf','ci')
         call mn902(cr(lstor),cr(lstor),ncor,'h times c',cr(lc0),
     $        cr(lgrad+nmix),'mcscr','mcscr','ci',junk,'mc')
      else if (citype.eq.'m903') then
         call mn903(cr(lstor),cr(lstor),ncor,'h times c',
     $        cr(lc0),cr(lgrad+nmix),'mcscr','mcscr','ci',junk,'mc')
      end if
c
c-------------------------
c     calculate ci energy
c-------------------------
      energy = 0.d0
      lc00 = lc0 - 1
      lg0 = lgrad + nmix - 1
      do 310 i=1, ncsf
 310     energy=energy + cr(lc00+i) * cr(lg0+i)
c
      if(prtflg.ne.'minimum') then
         write (iout,3001) energy
 3001    format(/,'  energy(mcledr)',1x,f24.16)
      endif
c
c-----------------------------------
c     read ci diagonal elements
c-----------------------------------
c
      ld0 = ldiag + nmix - 1
      call iosys('read real "h diagonals" from rwf',ncsf,
     $     cr(ld0+1),0,' ')
c
c----------------------------------------------
c     modify ci gradient and diagonal elements
c----------------------------------------------
c
      do 510 i = 1, ncsf
         cr(lg0+i) = cr(lg0+i) - energy * cr(lc00+i)
         cr(ld0+i) = cr(ld0+i) - energy
 510  continue
c
      endif
c
c..bhl..scf.end
c
c---------------------------------
c     initialize mcmxvc
c---------------------------------
c
      lbuf = lstor + 1
      lstor2 = lbuf + lenb
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(wpadti(lstor2),cr,ncore,'mc linear equations 2',0)
c
      ncor2 = iadtwp(ncore) - lstor2
c
      naobt = 0
      do 610 i = 1, nsym
 610  naobt = naobt + naob(i)
      lij = naobt * (naobt + 1) / 2
      lijkl=(lij*(lij+1))/2
c
c     call enmxvc(nsym,nbf,nob,nfob,ncob,naob,cv,cr(lc0),
c    $     locsym,len,locp,mix,
c    $     lij,lijkl,
c    $     nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
c    $     nfsm,ndshd1,ndshd2,nhd,
c    $     cr(lbuf),cr(lstor2),icr(wpadti(lstor2)),
c    $     ncor2)
c
c-------------------------------------------------------------
c     read diagonal elements of hessian and orbital gradient
c-------------------------------------------------------------
c
      call rdiagh(cr(ldiag),cr(lbuf),lenb,nmix)
c
      call iosys('rewind mcscf_gradient on rwf',0,0,0,' ')
      call iosys('read real mcscf_gradient from rwf without rewinding',
     $     nmix,cr(lgrad),0,' ')
c
c---------------------------------
c     solve linear equations
c---------------------------------
c
      call mcjcbi(cr(lgrad),cr(ldiag),cr(lt),b,cr(lc0),
     $     cr(lgl),cr(lolap),cr(lhess),cr(lthess),cr(lc),
     $     nmix,ncsf,mdim,energy,
     $     ndstv,nf45,nf46,lb46,noci,npassj,istopj,stopj,sqcdf,
     $     bufix,lbufso,cr(labcx),cr(labix),cr(laibx),
     $     cr(ihss),incor,cr(lsg),cr(ltg),
     $     nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,locp,mix,
     $     lij,lijkl,
     $     nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
     $     nfsm,ndshd1,ndshd2,nhd,
     $     cr(lbuf),cr(lstor2),icr(wpadti(lstor2)),
     $     ncor2)
c
c
      return
      end
