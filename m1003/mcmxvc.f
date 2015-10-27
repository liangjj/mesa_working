*deck @(#)mcmxvc.f	5.1  11/6/94
      subroutine mcmxvc(b,nmix,mdim,energy,noci,iter,t,bufix,lbufso,
     $     rabcx,rabix,raibx,hss,incor,sg,tg,
     $     nsym,nbf,nob,nfob,ncob,naob,cv,c0,
     $     locsym,len,lok,mix,
     $     lij,lijkl,
     $     nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
     $     nfsm,ndshd1,ndshd2,nhd,
     $     buf,cr,icr,numdp,
     $     wener,wt,iptci,ntwt,lmixt,thc,tgrad,
     $     ntcsf,nwtm,lmixm,nspinm,ncsfm,nlcsf,navg)
c***begin prologue     mcmxvc
c***date written       871022   (yymmdd)
c***revision date      921004   (yymmdd)
c
c   4 october  1992    rlm at lanl
c     explicitly declaring aint to be external.  the stardent confused
c     this subroutine with a intrinsic function.
c  17 april    1990    rlm at lanl
c     removing entry points
c  16 november 1987    bhl at lanl
c     core allocation bug fixed.
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcmxvc.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       mcmxvc
c
      implicit real*8 (a-h,o-z)
c
      external aint
c
      character*8 citype
      character*8 mcscf
c
      dimension nbf(2),nob(2),ncob(2),naob(2),cv(2),c0(2)
      dimension locsym(2),len(2),lok(2),mix(2),nfob(2)
      dimension buf(2),cr(2),icr(2)
      dimension b(2),t(2)
      logical debug
      real*8 bufix(lbufso),sg(*),tg(*)
      real*8 rabcx(*),rabix(*),raibx(*),hss(*)
c
      real*8 wt(*),wener(*)
      real*8 thc(*),tgrad(*)
      integer iptci(*),lmixm(*),ncsfm(*),nspinm(*),nwtm(*)
c
      common /io/ inp,iout
c
      parameter (debug=.false.)
c
c
      call iosys('read character mcscf_type from mcscr',-1,0,0,mcscf)
      call iosys('read character "mcscf: ci used" from rwf',
     $     -1,0,0,citype)
c
      if(navg.le.1) then
         ncsf = mdim - nmix
      else
         ncsf=ncsfm(1)
      end if
c-----------------------------------------------------------------------
c
c     b * delta
c
c-----------------------------------------------------------------------
c
      do 210 i=1,mdim
         t(i)=0.d0
 210  continue
      if(incor.lt.2) then
         call hesvec(t,b,hss,lenb,nmix)
      else
         call ebc(t,hss,b,nmix,nmix,1)
      end if
c
      if(debug) then
         write(iout,1000) (t(iii),iii=1,nmix)
 1000    format(/,'   mcmxvc:t=b x delta',/4(2x,f14.8))
      endif
c
      if(noci.eq.1) return
c
c-----------------------------------------------------------------------
c
c     a * gamma
c
c-----------------------------------------------------------------------
      intape=2
      iflag=1
      junk=0
c..avg
      if(navg.le.1) then
         if (citype.eq.'m902') then
            call mn820(cr,cr,numdp,'mcscr','mcscr','mcscf','ci')
            call mn902(cr,cr,numdp,'h times c',b(nmix+1),t(nmix+1),
     $                 'mcscr','mcscr','ci',junk,'mc')
         else if (citype.eq.'m903') then
            call mn903(cr,cr,numdp,'h times c',b(nmix+1),t(nmix+1),
     $                 'mcscr','mcscr','ci',junk,'mc')
         end if
c
         do 310 i = 1, ncsf
 310     t(nmix+i) = t(nmix+i) - energy * b(nmix+i)
c
      else
c
         nm=nmix+1
         do 312 iavg=1,navg
c
            if (citype.eq.'m902') then
               call mn820(cr,cr,numdp,'mcscr','mcscr','mcscf','ci')
               call mn902(cr,cr,numdp,'h times c',b(nm),t(nm),
     $                   'mcscr','mcscr','ci',junk,'mc')
            else if (citype.eq.'m903') then
               call mn903(cr,cr,numdp,'h times c',b(nm),t(nm),
     $                    'mcscr','mcscr','ci',junk,'mc')
            end if
c
            do 311 i = 1, ncsf
 311        t(nm-1+i)=t(nm-1+i)-wener(iavg)*b(nm-1+i)
c
            nm=nm+ncsf
 312     continue
      end if
      if(debug) then
         write(iout,1010) (t(iii),iii=1,nmix+ncsf)
 1010    format(/,'   mcmxvc:t=a x gamma',/4(2x,f14.8))
      endif
c..avg
c
c-----------------------------------------------------------------------
c
c     c * delta
c
c-----------------------------------------------------------------------
 400  lstor1 = 1
      ncor1 = numdp - lstor1
c-------------------------------
c     generate updated integrals
c-------------------------------
      mmda1=lda1
      intape=50
      call mcgvca(nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,lok,mix,b,lijkl,
     $     nda1,lda1,nf35,nf36,nf37,ltrb,
     $     cr(lstor1),cr(lstor1),ncor1,intape,bufix,lbufso,
     $     rabcx,rabix,raibx,incor,sg)
      lda1=mmda1
c
c------------------------------------
c     call h*c with updated integrals
c------------------------------------
      ncor2 = numdp-ncsf
      nc1=ncsf+1
      if(navg.le.1) then 
         if (citype.eq.'m902') then
            call mn820(cr(nc1),cr(nc1),ncor2,'mcscr','mcscr',
     $                 'mcscf','mcscf')
            call mn902(cr(nc1),cr(nc1),ncor2,'h times c',c0,
     $                 cr,'mcscr','mcscr','mcscf',junk,'mc')
         else if (citype.eq.'m903') then
            call mn903(cr(nc1),cr(nc1),ncor2,'h times c',c0,cr,
     $                 'mcscr','mcscr','mcscf',junk,'mc')
         end if
c
         call vadd(t(nmix+1),t(nmix+1),cr,ncsf)
c
      else
c..avg
c        state averaged code
c..avg
         ia=1
         nm=nmix+1
         do 411 iavg=1,navg
            if (citype.eq.'m902') then
               call mn820(cr(nc1),cr(nc1),ncor2,'mcscr','mcscr',
     $                    'mcscf','mcscf')
               call mn902(cr(nc1),cr(nc1),ncor2,'h times c',c0(ia),
     $                    cr,'mcscr','mcscr','mcscf',junk,'mc')
            else if (citype.eq.'m903') then
               call mn903(cr(nc1),cr(nc1),ncor2,'h times c',c0(ia),
     $                    cr,'mcscr','mcscr','mcscf',junk,'mc')
            end if
            call vadd(t(nm),t(nm),cr,ncsf)
            nm=nm+ncsf
            ia=ia+ncsf 
  411    continue
      end if
      if(debug) then
         write(iout,1020) (t(iii),iii=1,nmix+ncsf)
 1020    format(/,'   mcmxvc:t=c x delta',/4(2x,f14.8))
      endif
c-----------------------------------------------------------------------
c
c     c(t) * gamma
c
c-----------------------------------------------------------------------
c
      if(navg.le.1) then
         lg1 = lstor1
         lg2 = lg1 + lij
         lag2 = lg2 + lijkl
         lstor2 = lag2 + lij*lij
         ncor2 = numdp - lstor2
      else
         lg1 = lstor1
         lag1 = lg1 + lij      
         lg2 = lag1 + lij
         lag2 = lg2 + lij*lij
         lstor2 = lag2 + lij*lij
         ncor2 = numdp - lstor2
      end if
c-----------------------------------------
c     generate transition density matrices
c-----------------------------------------
c..avg
      need=lstor2
      call getscm(need,icr,junk,'dcan',0)
c
      if(navg.le.1) then
         if (citype.eq.'m902') then
            call mn902(cr(lstor2),cr(lstor2),ncor2,
     $                 'transition density',c0,b(nmix+1),'mcscr',
     $                 'mcscr','mcscf',junk,'mc')
            call mn830(cr(lstor2),cr(lstor2),ncor2,'mcscr','mcscr',
     $                'mcscf','guga transition density matrix',
     $                'mo 1ptdm','mo 2ptdm')
         else if (citype.eq.'m903') then
            call mn903(cr(lstor2),cr(lstor2),ncor2,'transition density',
     $                 c0,b(nmix+1),'mcscr','mcscr',' ',junk,'mc')
         end if
c
c  ----- first .. read transition density
c
         call iosys('read real "mo 1ptdm" from mcscr',lij,cr(lg1),
     $              0,' ')
         call iosys('read real "mo 2ptdm" from mcscr',
     $               lij**2,cr(lag2),0,' ')
      else
         ia=1
         nm=nmix+1
         lijlij=lij*lij
         call rzero(cr(lg1),lij)
         call rzero(cr(lag2),lijlij)
         do 511 iavg=1 ,navg
            if (citype.eq.'m902') then
               call mn902(cr(lstor2),cr(lstor2),ncor2,
     $                    'transition density',c0(ia),b(nm),
     $                    'mcscr','mcscr','mcscf',junk,'mc')
               call mn830(cr(lstor2),cr(lstor2),ncor2,'mcscr','mcscr',
     $                    'mcscf','guga transition density matrix',
     $                    'mo 1ptdm','mo 2ptdm')
            else if (citype.eq.'m903') then
               call mn903(cr(lstor2),cr(lstor2),ncor2,
     $                    'transition density',c0(ia),b(nm),
     $                    'mcscr','mcscr',' ',junk,'mc')
           end if
c
           ia=ia+ncsf
           nm=nm+ncsf
c
c  ----- first .. read transition density
c
           call iosys('read real "mo 1ptdm" from mcscr',lij,cr(lag1),
     $                0,' ')
           call iosys('read real "mo 2ptdm" from mcscr',
     $                 lij**2,cr(lg2),0,' ')
c
c  lijlij = lij*lij
c
            do 512 i=1,lij
               cr(lg1-1+i)=cr(lg1-1+i)+wt(iavg)*cr(lag1-1+i)
 512        continue 
            do 513 i=1,lijlij
               cr(lag2-1+i)=cr(lag2-1+i)+wt(iavg)*cr(lg2-1+i)
 513        continue 
 511     continue
      end if 
c---------------------------------------------
c     call cg with transition density matrices
c---------------------------------------------
c
      call dcan(junk,cr(lag2),cr(lg2),naob(1),lij,(lij+1)*lij/2)
c
c     1pdm from guga needs to be scaled
c
      call scdm1(cr(lg1),naob(1),lij)
c
c
      call mccga(nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,lok,mix,cr(lg1),cr(lg2),
     $     nda1,lda1,nf35,nf36,
     $     t,cr(lstor2),cr(lstor2),ncor2,bufix,lbufso,
     $     rabcx,rabix,raibx,incor,sg,tg)
c
      lda1=mmda1
c
c--------------------------------------------------c
c    scale the results of the matrix times vector
c    with energy weights
c--------------------------------------------------c
      if(navg.gt.1) then
         nm1=nmix+1
         mx=nm1
         int=nmix+nlcsf+1
         do 20010 j=1,navg
            call sscal(ncsf,wt(j),t(mx),1)
            mx=mx+ncsf
20010    continue
      end if
      if(lmixt.ne.0) then
c        internal-external ci block
         call ainex(b(nm1),b(int),t(nm1),t(int),thc,wt,iptci,
     1              ncsfm,nwtm,lmixm,nspinm)
c        internal ci portion of the hessian
         call aint(t(int),b(int),wt,wener,iptci,lmixm,nwtm,nspinm)
c        internal ci coupling block
         call cint(t,t(int),b,b(int),tgrad,wt,iptci,nmix,lmixm,
     $             nspinm)
      endif
c
c
      return
      end
