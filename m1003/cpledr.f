*deck @(#)cpledr.f	5.1  11/6/94
      subroutine cpledr(nsym,nbf,nob,nfob,ncob,naob,cv,
     $      nmix,ncsf,mdim,noci,istopj,stopj,b,sqcdf,
     $      cr,icr,ncore,lbufso,nograd,ndf,locsym,locp,len,mix,
     $      prtena,mcroot,
     $      nwtm,wt,ncsfm,nspinm,navg,iptci,lmixm,lmixt,nlcsf)
c***begin prologue     cpledr
c***date written       880322   (yymmdd)  
c***revision date      920210   (yymmdd)
c   10 february 1992  rlm at lanl
c      adding bhl state-averaged cpmc.
c   22 march 1988     bhl at llnl
c      read ci vector from "mc root" instead of "ci root"
c
c***keywords           
c***author             lengsfield, byron 
c***source             @(#)cpledr.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       cpledr
c
      implicit real*8(a-h,o-z)
      integer wpadti
      character*4 itoc
      character*8 mcscf,citype,prtflg
      character*8 answer
      logical prtena
      logical debug
c
      integer nbf(2),nob(2),ncob(2),naob(2),nfob(2)
      real*8 cv(2)
      real*8 cr(2)
      integer icr(2)
      real*8 b(2)
      integer locp(2),len(2),mix(2),locsym(2)
c..samc
      integer iptci(2,10,10),lmixm(10),
     $        nwtm(*),ncsfm(*)
      real*8 wt(10,*),wener(10,10)
c..samc
      data kount/0/
      save kount
      parameter (debug=.false.)
c
      common /io/ inp,iout
      common /number/zero,pt5,one,two,four,eight
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      call iosys('read character mcscf_type from mcscr',-1,0,0,mcscf)
c
c..bhl..scf
      if(noci.eq.0) then
         call iosys('read character "mcscf: ci used" from rwf',
     $              -1,0,0,citype)
         if(kount.eq.0) then
            if(debug) then
               write(iout,10) citype
  10           format(' mcscf type = ',a8)
            endif
            call iosys('read real "nuclear repulsion energy" from rwf',
     $                 1,enuc,0,' ')
            call iosys('does "frozen core energy" exist on rwf',0,0,0,
     $                 answer)
c.
            if (answer.eq.'no') then
               fzcore=0.0d+00
            else
               call iosys('read real "frozen core energy" from rwf',1,
     #                    fzcore,0,' ')
            end if
         end if
      endif
c
c     ----- storage allocation for jacobi -----
      if(navg.ne.1) then
         nroots=navg
      else
         nroots=mcroot
      end if
c
      nao=naob(1)
      nnao=nao*(nao+1)/2
      ntnao=nnao*nnao
c
      lij=nnao
      lijkl=nnao*(nnao+1)/2
c
      npassj=1
      lenb=lbufso
c
c     allocate core
      call getscm(0,cr,maxcor,'size',0)
      maxcor=iadtwp(maxcor-100000)
c
c..samc
      if(navg.le.1) then
         ibufix=1
         lgrad=ibufix+lbufso
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
c..samc
      else
c..samc
         itgrad=1
         ithc=itgrad+lmixt*nmix
         ibufix=ithc+nlcsf
         lgrad=ibufix+lbufso
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
c..samc
c
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(wpadti(lstor),cr,ncore,'mc linear equations',0)
      ncor = iadtwp(ncore)- lstor
c
      if(incor.eq.2) then
         call gethes(cr(ihss),lhss,cr(ibf),lenb)
      end if
c
      if(incor.ne.0) then
          call getrab(cr(labcx),cr(ibufix),lbufso,'abcx',
     #                naob(1),nob(1),naob(1))
          if(ncob(1).ne.0) then
             call getrab(cr(labix),cr(ibufix),lbufso,'abix',
     #                   ncob(1),nob(1),naob(1))
          endif
      endif
c
c
c
c     ----- read ci vector -----
      if(noci.eq.0) then
c..samc
         if(navg.le.1) then
            call iosys('read real "mc root '//itoc(nroots)
     #                  //'" from rwf',ncsf,cr(lc0),0,' ')
c
c           ----- calculate ci gradient -----
            zero=0.d0
            junk=0
            if (citype.eq.'m902') then
               call mn820(cr(lstor),cr(lstor),ncor,
     $                    'mcscr','mcscr','mcscf','ci')
               call mn902(cr(lstor),cr(lstor),ncor,'h times c',
     $                    cr(lc0),cr(lgrad+nmix),
     $                    'mcscr','mcscr','ci',junk,'mc')
            else if (citype.eq.'m903') then
               call mn903(cr(lstor),cr(lstor),ncor,'h times c',
     $                    cr(lc0),cr(lgrad+nmix),
     $                    'mcscr','mcscr','ci',junk,'mc')
            end if
c
c           ----- calculate ci energy -----
            energy = 0.d0
            lc00 = lc0 - 1
            lg0 = lgrad + nmix - 1
            do 310 i = 1, ncsf
  310       energy = energy + cr(lc00+i) * cr(lg0+i)
c
            if(prtena.or.debug) then
               write (iout,3001) energy
 3001          format(5x,'cpledr: energy = ',1x,f24.16)
            end if
c
c          ----- read ci diagonal elements -----
            ld0 = ldiag + nmix - 1
            call iosys('read real "h diagonals" from rwf',ncsf,
     #                 cr(ld0+1),0,' ')
c
c..samc(noci)     endif
c
c           ----- modify ci gradient and diagonal elements -----
            do 510 i = 1, ncsf
               cr(ld0+i) = cr(ld0+i) - energy
  510       continue
c
            if(nograd.eq.0) then
               do 511 i = 1, ncsf
                  cr(lg0+i) = cr(lg0+i) - energy * cr(lc00+i)
  511          continue
            end if
c..samc
         else
c..samc
            write(iout,*)' state-averaged set-up '
c           ----- read ci vector -----
            ewt=0.d0
            mroots=0
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
                  mroots=mroots+1
                  call iosys('read real "mc root '//itoc(mroots)
     $                       //'" from rwf',ncsf,cr(lcw+1),0,' ')
c
c                 ----- calculate ci gradient -----
                  zero=0.d0
                  junk=0
                  if (citype.eq.'m902') then
                     call mn820(cr(lstor),cr(lstor),ncor,
     $                          'mcscr','mcscr','mcscf','ci')
                     call mn902(cr(lstor),cr(lstor),ncor,'h times c',
     $                          cr(lcw+1),cr(lgw+nmix+1),'mcscr',
     $                          'mcscr','ci',junk,'mc')
                  else if (citype.eq.'m903') then
                     call mn903(cr(lstor),cr(lstor),ncor,'h times c',
     $                          cr(lcw+1),cr(lgw+nmix+1),'mcscr',
     $                          'mcscr','ci',junk,'mc')
                  end if
c
c                 ----- calculate ci energy -----
                  energy = 0.d0
                  lg0 = lgw + nmix 
                  do 311 i = 1, ncsf
 311                 energy = energy + cr(lcw+i) * cr(lg0+i)
                  wener(j,k)=energy
                  ewt=ewt+wt(j,k)*energy
c
                  if(kount.eq.0) then
                     write (iout,3011) j,energy,energy+enuc+fzcore
 3011                format(/,' mcledr: root  e(elec)  e(total)',
     $                      i4,2(3x,f18.12))
                  end if
c
c                 ----- read ci diagonal elements -----
                  ld0 = ldw + nmix 
                  call iosys('read real "h diagonals" from rwf',ncsf,
     $                        cr(ld0+1),0,' ')
c
c                 ----- modify ci gradient and diagonal elements -----
c                 write(iout,*)' diagonals ',j
c                 write(iout,55144)(cr(ld0+i),i=1,ncsf)
c55144            format(5(2x,f12.8))
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
c
            if(kount.eq.0) then
               kount=1
               write (iout,33001) ewt+enuc+fzcore
33001          format(/,' mcledr: state-averaged energy ',38x,f18.12,/)
            end if
c
            if(lmixt.eq.0)go to 88006
c           write(iout,91001)
c91001      format(/,'  making tgrad ')
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
               if(lmix.ne.0) then 
                  do 19040 m=1,lmix
                     ici=(iptci(1,m,k)-1)*ncsf+lc0+ix
                     jci=(iptci(2,m,k)-1)*ncsf+lc0+ix
                     if (citype.eq.'m902') then
                        call mn902(cr(lstor2),cr(lstor2),ncor2,
     $                            'transition density',cr(ici),cr(jci),
     $                            'mcscr','mcscr','mcscf',junk,'mc')
                        call mn830(cr(lstor2),cr(lstor2),ncor2,'mcscr',
     $                             'mcscr','mcscf',
     $                             'guga transition density matrix',
     $                            'mo 1ptdm','mo 2ptdm')
                     else if (citype.eq.'m903') then
                        call mn903(cr(lstor2),cr(lstor2),ncor2,
     $                             'transition density',cr(ici),cr(jci),
     $                             'mcscr','mcscr',' ',junk,'mc')
                     end if
c
c                    ----- first .. read transition density
c
                     call iosys('read real "mo 1ptdm" from mcscr',
     $                           lij,cr(lg1),0,' ')
                     call iosys('read real "mo 2ptdm" from mcscr',
     $                           lij**2,cr(lgg2),0,' ')
c
c                     ----- call cg with transition density matrices
                     call dcan(junk,cr(lgg2),cr(lg2),naob(1),lij,
     $                         (lij+1)*lij/2)
c
c                    ----- 1pdm from guga needs to be scaled
                     call scdm1(cr(lg1),naob(1),lij)
c
c                    ----- call cg with transition density matrices
                     call rzero(cr(itgrad+jx),nmix)
                     call mccga(nsym,nbf,nob,nfob,ncob,naob,cv,
     $                          locsym,len,locp,mix,cr(lg1),cr(lg2),
     $                          nda1,lda1,nf35,nf36,
     $                          cr(itgrad+jx),cr(lstor3),icr(lstor3),
     $                          ncor3,bufix,lbufso,cr(labcx),cr(labix),
     $                          cr(laibx),incor,cr(lsg),cr(ltg))
c
                     jx=jx+nmix
c
19040             continue
               endif
               ix=ix+nwt*ncsf
19140       continue
88006       continue
cbhl
c..samc
         end if
      end if
c
      if(debug) then
         write(iout,91929)(cr(lg0+i),i=1,ncsf)
91929    format(/,' ci gradient ',6(1x,f15.10))
      endif
c
c     ----- initialize mcmxvc -----
      lbuf = lstor + 1
      lstor2 = lbuf + lenb
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(wpadti(lstor2),cr,ncore,'mc linear equations 2',0)
      ncor2 = iadtwp(ncore) - lstor2
c
      naobt = 0
      do 610 i = 1, nsym
  610    naobt = naobt + naob(i)
      lij = naobt * (naobt + 1) / 2
      lijkl=(lij*(lij+1))/2
c
c     call enmxvc(nsym,nbf,nob,nfob,ncob,naob,cv,cr(lc0),
c    1      locsym,len,locp,mix,
c    2      lij,lijkl,
c    3      nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
c    4      nfsm,ndshd1,ndshd2,nhd,
c    5      cr(lbuf),cr(lstor2),icr(lstor2),
c    6      ncor2)
c
c     ----- read diagonal elements of hessian and orbital gradient
      call rdiagh(cr(ldiag),cr(lbuf),lenb,nmix)
c
      if(nograd.eq.0) then
         write(iout,*) ' reading test gradient mdim ',mdim
         call iosys('rewind test_gradient on rwf',0,0,0,' ')
         call iosys('read real test_gradient from rwf'
     $        //' without rewinding',mdim,cr(lgrad),0,' ')
         write(iout,22133)(cr(lgrad-1+li),li=1,mdim)
22133    format(5(2x,f12.8))
      else
         lg0=lgrad-1
         do 821 k=1,mdim
            cr(lg0+k)=b(k)
  821    continue
      end if
c
c
      if(noci.eq.0) then
      call iosys('read real "h diagonals" from rwf',ncsf,
     $     cr(lt),0,' ')
      end if
c
c     ----- solve linear equations -----
c
      call mcjcbi(cr(lgrad),cr(ldiag),cr(lt),b,cr(lc0),
     $     cr(lgl),cr(lolap),cr(lhess),cr(lthess),cr(lc),
     $     nmix,ncsf,mdim,energy,
     $     ndstv,nf45,nf46,lb46,noci,npassj,istopj,stopj,sqcdf,
     $     cr(ibufix),lbufso,cr(labcx),cr(labix),cr(laibx),
     $     cr(ihss),incor,cr(lsg),cr(ltg),
     $     cr(ithc),wener,wt,iptci,nwtm,lmixm,ncsfm,nspinm,
     $     lmixt,nlcsf,navg,
     $     nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,locp,mix,
     $     lij,lijkl,
     $     nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
     $     nfsm,ndshd1,ndshd2,nhd,
     $     cr(lbuf),cr(lstor2),icr(wpadti(lstor2)),ncor2,
     $     ntwt,cr(itgrad),ntcsf)

c
      if(debug) then
         write (iout,603) (b(i),i=1,mdim)
  603    format('0*mcledr solution '//4(1x,f16.8))
      endif
c
c
      return
      end
