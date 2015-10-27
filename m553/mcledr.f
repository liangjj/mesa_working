*deck @(#)mcledr.f	5.1  11/6/94
      subroutine mcledr(nsym,nbf,nob,nfob,ncob,naob,cv,g1,g2,
     $     locsym,len,locp,mix,nmix,
     $     nf16,nfsm,ndshd1,ndshd2,nhd,nda1,lda1,nf41,lenb,
     $     nf35,nf36,nf37,nf45,nf46,lb46,ndstv,
     $     noci,npassj,istopj,stopj,
     $     b, sqcdf,
     $     cr,icr,ncore,ncsf,bufix,lbufso,mcroot,
     $     nwtm,wt,ncsfm,nspinm,navg)
c***begin prologue     mcledr
c***date written       871022   (yymmdd)
c***revision date      920202   (yymmdd)
c
c    2 february 1992   rlm at lanl
c      adding bhl state-averaged mcscf
c   17 april    1990   rlm at lanl
c      removing need for the enmxmc entry point.
c   17 november 1987   bhl at lanl
c     skip over ci-gradient routine if noci > 0
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcledr.f	5.1   11/6/94
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
      character*3 answer
      character*4 itoc
      character*8 citype
      character*8 mcscf,prtflg
c
      integer wpadti
      integer nbf(2),nob(2),ncob(2),naob(2)
      integer locsym(2),len(2),locp(2),mix(2),nfob(2)
      integer icr(*)
      integer lmx(10,10,10), kmx(10,10,10), iptci(2,10,10)
      integer lmixm(10),nwtm(*),ncsfm(*)
      real*8 wt(10,*),wener(10,10)
      real*8 bufix(lbufso)
      real*8 cv(2),g1(2),g2(2)
      real*8 cr(*),b(*)
      logical debug,keyval
c
      common /io/ inp,iout
c
      parameter (debug=.false.)
c
      if(debug) then
         write(iout,*) 'in mcledr'
      endif
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      call iosys('read character mcscf_type from mcscr',-1,0,0,mcscf)
      if(noci.eq.0) then
         call iosys('read character "mcscf: ci used" from rwf',
     $              -1,0,0,citype)
      else
         ncsf=0
      endif
c
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,enuc,0,' ')
      call iosys('does "frozen core energy" exist on rwf',0,0,0,
     $           answer)
      if (answer.eq.'no') then
         fzcore=0.0d+00
      else
         call iosys('read real "frozen core energy" from rwf',1,
     #               fzcore,0,' ')
      end if

c
c     ----- setup ci mixings pointers -----
      if(navg.gt.1) then
         if(debug) then
            write(iout,*)' setting-up state-averaged ci pointers '
         endif
         kmixt=0
         lmixt=0
         ntcsf=0
         nlcsf=0
         ntwt=0
         do 19003 k=1,nspinm
            ntcsf=ntcsf+ncsfm(k)
            nwt=nwtm(k)
            ntwt=ntwt+nwt
            lmixm(k)=0
            nlcsf=nlcsf+nwt*ncsfm(k)
            kmix=0
            lmix=0
            if(nwt.ne.1) then
               nwt1=nwt-1
               do 19001 i=1,nwt1
                  i1=i+1
                  do 19002 j=i1,nwt
                     if(wt(i,k).eq.wt(j,k)) then
                        kmix=kmix+1
                        kmx(i,j,k)=kmix
                        kmx(j,i,k)=kmix
                     else
                        lmix=lmix+1
                        lmx(i,j,k)=lmix
                        lmx(j,i,k)=lmix
                        iptci(1,lmix,k)=i
                        iptci(2,lmix,k)=j
                     endif
19002             continue
19001          continue
               lmixm(k)=lmix
               lmixt=lmixt+lmix
               kmixt=kmixt+kmix
            endif
19003    continue
      end if
c
c     ----- storage allocation for jacobi -----
c
cbhl
      if (.not.keyval(ops,'nroots',nroots))     nroots=1
      if (mcroot.ne.0) nroots=mcroot
cbhl
c
c
      call getscm(0,cr,maxcor,'size',0)
      maxcor=iadtwp(maxcor-100000)
c..samc
      if(navg.le.1) then
         nspinm=1
         nwtm(1)=1
         ncsfm(1)=ncsf
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
      else
c..samc
         mdim = nmix + nlcsf + lmixt
         itgrad=1
         ithc=itgrad+lmixt*nmix
         lgrad = ithc+nlcsf
         ldiag = lgrad + mdim
         lt = ldiag + mdim
         lc = lt + mdim*(istopj+1)
         lc0 = lc + mdim
         lhess = lc0 + nlcsf
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
      end if
c
      naobt = 0
      do 610 i = 1, nsym
 610     naobt = naobt + naob(i)
      lij = naobt * (naobt + 1) / 2
      lijkl=(lij*(lij+1))/2
c
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
     $                  naob(1),nob(1),naob(1))
cbhl
c fix for ncob = 0
cbhl
            if(ncob(1).ne.0) then
               call getrab(cr(labix),bufix,lbufso,'abix',
     $                     ncob(1),nob(1),naob(1))
            endif
c
c      call addrab(cr(labix),bufix,lbufso,'aibx',
c    #   ncob(1),nob(1),naob(1))
c
         endif
      end if
c
c..bhl..scf.start
c
      if(noci.eq.0) then
         if(navg.le.1) then
cc
c--------------------------------
c           read ci vector
c--------------------------------
cc
            call iosys('read real "mc root '//itoc(nroots)
     $               //'" from rwf',ncsf,cr(lc0),0,' ')
c
c---------------------------------
c           calculate ci gradient
c---------------------------------
c
            zero=0.d0
            junk=0
            if (citype.eq.'m902') then
               call mn820(cr(lstor),cr(lstor),ncor,'mcscr',
     $                    'mcscr','mcscf','ci')
               call mn902(cr(lstor),cr(lstor),ncor,'h times c',
     $                    cr(lc0),cr(lgrad+nmix),'mcscr','mcscr',
     $                    'ci',junk,'mc')
            else if (citype.eq.'m903') then
               call mn903(cr(lstor),cr(lstor),ncor,'h times c',
     $                    cr(lc0),cr(lgrad+nmix),'mcscr','mcscr',
     $                    'ci',junk,'mc')
            end if
c
c-------------------------
c           calculate ci energy
c-------------------------
            energy = 0.d0
            lc00 = lc0 - 1
            lg0 = lgrad + nmix - 1
            do 310 i = 1, ncsf
 310           energy = energy + cr(lc00+i) * cr(lg0+i)
c
            if(prtflg.ne.'minimum') then
               write (iout,3001) energy
 3001          format(/,'  energy(mcledr)',1x,f24.16)
            endif
c
c-----------------------------------
c           read ci diagonal elements
c-----------------------------------
c
            ld0 = ldiag + nmix - 1
            call iosys('read real "h diagonals" from rwf',ncsf,
     $                  cr(ld0+1),0,' ')
c
c----------------------------------------------
c        modify ci gradient and diagonal elements
c----------------------------------------------
c
            do 510 i = 1, ncsf
               cr(lg0+i) = cr(lg0+i) - energy * cr(lc00+i)
               cr(ld0+i) = cr(ld0+i) - energy
 510        continue
c..samc
         else
            if(debug) then
               write(iout,*)' state-averaged set-up '
            endif
c--------------------------------
c           read ci vector
c--------------------------------
cc
            ewt=0.d0
            nroots=0
c
            ithcw=ithc-1
            lcw=lc0-1
            lgw=lgrad-1
            itgw=itgrad-1
            ldw=ldiag-1
c
            do 19110 k=1,nspinm
               ncsf=ncsfm(k)
               nwt=nwtm(k)
c 
               do 19111 j=1,nwt
                  nroots=nroots+1
                  call iosys('read real "mc root '//itoc(nroots)
     $                     //'" from rwf',ncsf,cr(lcw+1),0,' ')
c
c---------------------------------
c                 calculate ci gradient
c---------------------------------
c
                  zero=0.d0
                  junk=0
                  if (citype.eq.'m902') then
                     call mn820(cr(lstor),cr(lstor),ncor,'mcscr',
     $                          'mcscr','mcscf','ci')
                     call mn902(cr(lstor),cr(lstor),ncor,'h times c',
     $                          cr(lcw+1),cr(lgw+nmix+1),'mcscr',
     $                          'mcscr','ci',junk,'mc')
                  else if (citype.eq.'m903') then
                     call mn903(cr(lstor),cr(lstor),ncor,'h times c',
     $                          cr(lcw+1),cr(lgw+nmix+1),'mcscr',
     $                          'mcscr','ci',junk,'mc')
                  end if
c
c-------------------------
c                 calculate ci energy
c-------------------------
                  energy = 0.d0
                  lg0 = lgw + nmix 
                  do 311 i = 1, ncsf
 311                 energy = energy + cr(lcw+i) * cr(lg0+i)
                  wener(j,k)=energy
                  ewt=ewt+wt(j,k)*energy
c
c                 write (iout,3011) j,energy,energy+enuc+fzcore
c3011             format(/,' mcledr: root  e(elec)  e(total) ',
c    $                   i4,2(3x,f18.12))
                  write (iout,3011) j,energy,energy+enuc+fzcore
 3011             format(26x,' root:',
     $                   i4,2(3x,f18.12))
c
c-----------------------------------
c                 read ci diagonal elements
c-----------------------------------
c
                  ld0 = ldw + nmix 
                  call iosys('read real "h diagonals" from rwf',ncsf,
     $                        cr(ld0+1),0,' ')
c
c----------------------------------------------
c                 modify ci gradient and diagonal elements
c----------------------------------------------
c
                  do 5510 i = 1, ncsf
                     cr(ithcw+i) = cr(lg0+i)
                     cr(ld0+i) = wt(j,k)*(cr(ld0+i) - energy)
 5510             continue
                  lgw=lgw+ncsf
                  ithcw=ithcw+ncsf
                  ldw=ldw+ncsf
                  lcw=lcw+ncsf
19111          continue
c
19110       continue
c
            call projct(cr(lgrad+nmix),cr(lc0),ncsfm,nwtm,nspinm)
c
            lg0=lgrad+nmix-1
            ix=1
            do 59130 k=1,nspinm
               nwt=nwtm(k)
               ncsf=ncsfm(k)
               do 59030 j=1,nwt
                  call sscal(ncsf,wt(j,k),cr(lg0+ix),1)
                  ix=ix+ncsf
59030          continue
59130       continue
c
            if(lmixt.eq.0) go to 19032
c
c  note.. this code must be changed if trial ci vectors
c         are not contrained to be eigenfunctions of the 
c         hamiltonian
c
            ix=nmix+nlcsf+lgrad-1
            do 19031 i=1,lmixt
               cr(ix+i)=0.d0
19031       continue
19032       continue
c
            if(debug) then
               write (iout,33001) ewt+enuc+fzcore
33001          format(/,' mcledr: state-averaged energy ',20x,f18.12)
            endif
c
            if(lmixt.eq.0)go to 88006
            write(iout,91001)
91001       format(/,'  making tgrad ')
c
            lstor2=lstor
            ncor2=ncor
c
            lg1=lstor2
            lg2=lg1+lij
            lgg2=lg2+lijkl
            lstor3=lgg2+lij*lij
            ncor3=ncore-lstor3
c
            jx=0
            ix=0
            do 19140 k=1,nspinm
               ncsf=ncsfm(k)
               lmix=lmixm(k)
               nwt=nwtm(k)
               if(lmix.eq.0)go to 19041
               do 19040 m=1,lmix
                  ici=(iptci(1,m,k)-1)*ncsf+lc0+ix
                  jci=(iptci(2,m,k)-1)*ncsf+lc0+ix
                  if (citype.eq.'m902') then
                     call mn902(cr(lstor2),cr(lstor2),ncor2,
     $                          'transition density',cr(ici),cr(jci),
     $                          'mcscr','mcscr','mcscf',junk,'mc')
                     call mn830(cr(lstor2),cr(lstor2),ncor2,'mcscr',
     $                          'mcscr','mcscf',
     $                          'guga transition density matrix',
     $                          'mo 1ptdm','mo 2ptdm')
                  else if (citype.eq.'m903') then
                     call mn903(cr(lstor2),cr(lstor2),ncor2,
     $                          'transition density',cr(ici),cr(jci),
     $                          'mcscr','mcscr',' ',junk,'mc')
                  end if
c
c                 ----- first .. read transition density
c
                  call iosys('read real "mo 1ptdm" from mcscr',
     $                       lij,cr(lg1),0,' ')
                  call iosys('read real "mo 2ptdm" from mcscr',
     $                        lij**2,cr(lgg2),0,' ')
c
c---------------------------------------------
c                 call cg with transition density matrices
c---------------------------------------------
c
                  call dcan(junk,cr(lgg2),cr(lg2),naob(1),lij,
     $                     (lij+1)*lij/2)
c
c                 1pdm from guga needs to be scaled
c
                  call scdm1(cr(lg1),naob(1),lij)
c
c---------------------------------------------
c                 call cg with transition density matrices
c---------------------------------------------
c
                  call rzero(cr(itgrad+jx),nmix)
                  call mccga(nsym,nbf,nob,nfob,ncob,naob,cv,
     $                       locsym,len,locp,mix,cr(lg1),cr(lg2),
     $                       nda1,lda1,nf35,nf36,
     $                       cr(itgrad+jx),cr(lstor3),icr(lstor3),ncor3,
     $                       bufix,lbufso,cr(labcx),cr(labix),cr(laibx),
     $                       incor,cr(lsg),cr(ltg))
c
                  jx=jx+nmix
c
19040          continue
19041          continue
               ix=ix+nwt*ncsf
19140       continue
88006       continue
cbhl
c..samc
         end if
c
      end if
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
c
c      call enmxvc(nsym,nbf,nob,nfob,ncob,naob,cv,cr(lc0),
c     $     locsym,len,locp,mix,
c     $     lij,lijkl,
c     $     nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
c     $     nfsm,ndshd1,ndshd2,nhd,
c     $     cr(lbuf),cr(lstor2),icr(lstor2),ncor2,
c     $     wener,wt,iptci,ntwt,lmixt,cr(ithc),cr(itgrad),
c     $     ntcsf,nwtm,lmixm,nspinm,nlcsf)
c
c-------------------------------------------------------------
c     read diagonal elements of hessian and orbital gradient
c-------------------------------------------------------------
c
      call rdiagh(cr(ldiag),cr(lbuf),lenb,nmix)
c
      call iosys('rewind mcscf_gradient on rwf',0,0,0,' ')
      call iosys('read real mcscf_gradient from rwf without rewinding',
     $           nmix,cr(lgrad),0,' ')
c
c..for wtdiag in state-average mc
c
      if(noci.eq.0) then
         call iosys('read real "h diagonals" from rwf',ncsf,
     $               cr(lt),0,' ')
      end if
c
c---------------------------------
c     solve linear equations
c---------------------------------
c
      call mcjcbi(cr(lgrad),cr(ldiag),cr(lt),b,cr(lc0),
     $            cr(lgl),cr(lolap),cr(lhess),cr(lthess),cr(lc),
     $            nmix,ncsf,mdim,energy,
     $            ndstv,nf45,nf46,lb46,noci,npassj,istopj,stopj,sqcdf,
     $            bufix,lbufso,cr(labcx),cr(labix),cr(laibx),
     $            cr(ihss),incor,cr(lsg),cr(ltg),
     $            cr(ithc),wener,wt,iptci,nwtm,lmixm,ncsfm,nspinm,
     $            lmixt,nlcsf,navg,
     $            nsym,nbf,nob,nfob,ncob,naob,cv,
     $            locsym,len,locp,mix,lij,lijkl,
     $            nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
     $            nfsm,ndshd1,ndshd2,nhd,
     $            cr(lbuf),cr(lstor2),icr(wpadti(lstor2)),ncor2,
     $            ntwt,cr(itgrad),ntcsf)
c
c
      return
      end
