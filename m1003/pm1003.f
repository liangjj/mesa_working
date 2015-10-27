*deck @(#)pm1003.f	5.1 11/6/94
      subroutine pm1003(cr,icr)
c***begin prologue     pm1003
c***date written       860501  (yymmdd)  
c***revision date      920210  (yymmdd)
c  10 february 1992    rlm at lanl
c      implementing bhl state-averaged code.
c  18 june     1991    rlm at lanl
c      creating separate disk files for derivative integrals,etc.
c   4 december 1986    pws at lanl
c      1. changing 'namint' to character variable.
c      2. changing iosys opens to character unit names
c      3. adding family-member-size to mcscr.
c***keywords           
c***author             lengsfield, byron 
c***source             @(#)pm1003.f	5.1 11/6/94
c***purpose            solves coupled perturbed mcscf equations. 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       pm1003
c
      implicit real*8(a-h,o-z)
      real*8 cr(*)
      integer icr(*)
c
      common /io/ inp,iout
c
      integer nobb(2),ncob(2),naob(2),nfob(2)
      character*8 mcscf,citype,prtflg,prtsav
      character*128 dints
      character*4 itoc
      character*4096 ops
      logical forc,debug,logkey,prtrhs,prtsol,prttag,prttah
      logical prtena,prtum,prtall,prtcp,prtfor
      integer outym,outle,outta
      integer wpadti
c..samc
      integer lmx(10,10,10),kmx(10,10,10),iptci(2,10,10),lmixm(10)
      integer nwtm(10),ncsfm(10)
      real*8 wt(10,10),wtaao(10)
c..samc
      real*8 zero,one
      parameter (zero=0.0d+00, one=1.0d+00)
c
  100 format(1x,'m1003:state-averaged coupled-perturbed mcscf',
     $     /,5x,'byron lengsfield  llnl may 1991')
  110 format(5x,'maximum # of cycles :',i8)
  120 format(5x,'convergence criteria:',e8.1)
  130 format(5x,'nograd              :',i8)
  140 format(5x,'force constant calculation')
  150 format(5x,'mcscf force constants stored on rwf')
c
c
      call getscm(0,cr,maxint,'size',0)
      ncore=iadtwp(maxint)
c
c     ----- temporarily turn off printing as much as possible -----
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtsav)
      prtflg='minimum'
      call iosys('write character "print flag" to rwf',0,0,0,prtflg)
c
      call iosys('read character options from rwf',-1,0,0,ops)
      istopj=intkey(ops,'cpmc=linear-equations=cycles',30,' ')
      lnbuf=intkey(ops,'cpmc=buffer-size',40000,' ')
      outym=intkey(ops,'cpmc=out-of-core-ymatrix',0,' ')
      outta=intkey(ops,'cpmc=out-of-core-taints',0,' ')
      outle=intkey(ops,'cpmc=out-of-core-linear-equations',0,' ')
      stopj=fpkey(ops,'cpmc=linear-equations=convergence',
     $            1.0d-06,' ')
      stopj=fpkey(ops,'cpmc=stopj',stopj,' ')
      istopj=intkey(ops,'cpmc=istopj',istopj,' ')
      debug=logkey(ops,'cpmc=debug',.false.,' ')
      prtrhs=logkey(ops,'cpmc=print=rhs',.false.,' ')
      prtsol=logkey(ops,'cpmc=print=solution',.false.,' ')
      prtena=logkey(ops,'cpmc=print=energy',.false.,' ')
      prttah=logkey(ops,'cpmc=print=ta_h',.false.,' ')
      prttag=logkey(ops,'cpmc=print=ta_g',.false.,' ')
      prtum=logkey(ops,'cpmc=print=u-matrix',.false.,' ')
      prtall=logkey(ops,'cpmc=print=all',.false.,' ')
      prtcp=logkey(ops,'cpmc=print=cpmc',.false.,' ')
      prtfor=logkey(ops,'cpmc=print=force',.false.,' ')
c
      if(logkey(ops,'cpmc=force',.false.,' ')) forc=.true.
      if(logkey(ops,'force-constants',.false.,' ')) then
         if(logkey(ops,'force-constants=numerical',.false.,' ')) then
c           numerical force constants, do nothing
         else
c           analytic force constants
            forc=.true.
         endif 
      endif
c
      if(prtall) prtcp=.true.
      if(prtcp)  prtfor=.true.
c
      write(iout,100)
      write(iout,110) istopj
      write(iout,120) stopj
c
      if(forc) then
         write(iout,140)
      end if
c
      call iosys('read integer mc_incorh from rwf',1,incorh,0,' ')
      call iosys('read integer mc_nbasis from rwf',1,nbf,0,' ')
      call iosys('read integer mc_ncore from rwf',1,nco,0,' ')
      call iosys('read integer mc_nactive from rwf',1,nao,0,' ')
      call iosys('read integer mc_norbs from rwf',1,nob,0,' ')
      call iosys('read integer mc_lbufso from rwf',1,lbufso,0,' ')
      call iosys('read integer mc_ncsf from rwf',1,ncsf,0,' ')
      call iosys('read integer mc_nmix from rwf',1,nmix,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $           1,natoms,0,' ')
c
      ndf=3*natoms
      if (debug) then
         write(iout,9001) nbf,nco,nao,nob,lbufso,ncsf,nmix,natoms
9001     format(' nbf nco nao nob lbufso ncsf nmix natoms',/,10(1x,i6))
      endif
c
c..samc
c
      navg=intkey(ops,'mcscf=navg',1,' ')
      mavg=intkey(ops,'cpmc=navg',1,' ')
      mcroot=intkey(ops,'mcscf=nroot',1,' ')
      mcroot=intkey(ops,'cpmc=nroot',mcroot,' ')
c
      noci=0
      if(logkey(ops,'mcscf=hf',.false.,' ')) then
         noci=1
         ncsf=0
      end if
c
      if(debug) then
         write(iout,*)' noci ',noci
      endif
c
      if(navg.gt.1) then
         call fparr(ops,'mcscf=weight',wt,navg,' ')
         mavg=0
      end if
      if(mavg.gt.1) then
         navg=mavg
         call fparr(ops,'cpmc=weight',wt,navg,' ')
      end if
c
      nspinm=1
      ncsfm(1)=ncsf
      nwtm(1)=navg
      if(navg.eq.1) then
         wt(1,1)=one
      else
         xx=zero
         do 11222 i=1,navg
            xx=xx+wt(i,1)
11222    continue
         xx=one/xx
         call sscal(navg,xx,wt,1)
         write(iout,11223)(wt(i,1),i=1,navg)
11223    format(/,' ** state-averaged weights ** ',/,
     $          5(2x,f10.7))
      end if
c
c     setup ci mixings pointers
      lmixt=0
      nlcsf=ncsf
      if(navg.gt.1) then
         write(iout,*)' setting-up state-averaged ci pointers '
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
      lncsf=navg*ncsf+lmixt
      if(debug) then
         write(iout,*)' lncsf nlcsf lmixt ',lncsf,nlcsf,lmixt      
      endif
c
c     allocate core
      if(navg.gt.1) then
         nroots=navg
      else
         nroots=mcroot
      end if
      noc=nco+nao
      nocnob=noc*nob
      nobnob=nob*nob
      nnao=nao*(nao+1)/2
      naonao=nao*nao
      ntnao=nnao*nnao
c
      nnob=nob*(nob+1)/2
c
      imixh=1
      icm=iadtwp(imixh+nocnob)
      ici=icm+nob*nob
      ibuf=ici+nlcsf
      iabcx=ibuf+lbufso
      iabix=iabcx+nnao*nao*nob
      iscr=iabix+nnao*nco*nob
      ifpq=iscr+nob*nob
      ir=ifpq+nob*nob
      isg=ir+nob*nob
      ilast=isg+nnao*naonao
c
      ilast=max(ilast,iscr+lbufso)
c
      incor=1
c
      if(ilast.gt.ncore.or.outta.ne.0) then
c
         incor=0
c
         iabcx=ibuf+lbufso
         iabix=iabcx+1
         iscr=iabix+1
         ifpq=iscr+nob*nob
         ir=ifpq+nob*nob
         isg=ir+nob*nob
         ilast=isg+nnao*naonao
c
         ilast=max(ilast,iscr+lbufso)
c
         if(ilast.gt.ncore) then
            call lnkerr(' m1003: ilast gt ncore ')
         endif
      end if
c
      nleft=ncore-ilast-1
      lblock=nob*noc+nnao+ntnao+2*lncsf+nnao*naonao
      igint=nnao+ntnao
      mdf=nleft/lblock
      mdf=min(mdf,ndf)
c
      if(mdf.lt.ndf) then
         write(iout,*)' *** multi-pass option invoked *** '
         write(iout,*)'  nleft  need ',nleft,lblock*ndf
      end if
c
      if(mdf.lt.1) then
         write(iout,*)'  lblock  nleft ',lblock,nleft
         call lnkerr(' increase core for m100: taint mdf lt 1 ')
      end if
c
      junk=ndf*(3*ncsf+nobnob+nocnob+2*nmix+naonao+nao)+250000
      call iosys('open mcscr as scratch on ssd',junk,0,0,' ')
      if(noci.eq.0) then
         call iosys('create real der_cigrad on mcscr',ndf*lncsf,0,0,' ')
         call iosys('create real der_civec on mcscr',ndf*lncsf,0,0,' ')
      end if
      call iosys('create real mo_ta_lag on mcscr',
     1            ndf*nocnob,0,0,' ')
      call iosys('create real der_grad on mcscr',ndf*nmix,0,0,' ')
c
      call iosys('read character "derivative integral filename"'
     1            //' from rwf',0,0,0,dints)
      call iosys('open dints as old',0,0,0,dints)
      call iosys('rewind all on dints',0,0,0,' ')
c
      call iosys('read real mcscf_ao_core_fock from rwf',
     1           nobnob,cr(iscr),0,' ')
c..bhl.bug
      call sqtotr(cr(ir),cr(iscr),nob,nnob)
c..bhl.bug
      call iosys('write real mc_core_fock to mcscr',
     1            nnob,cr(ir),0,' ')
c
      ncota=ilast
      naota=ncota+mdf*nco*nob
      ig1=naota+mdf*nob*nao
      ig2=ig1+mdf*nnao
      itaci=ig2+mdf*ntnao
      itbci=itaci+mdf*lncsf
      itg=itbci+mdf*lncsf
      ilast4=itg+mdf*nnao*naonao
c
      ilast4=max(ilast4,iscr+lbufso)
      need4=ilast4+2*lncsf+nnao+ntnao
      call getscm(need4,cr,ncor,'taint',0)
      mcore=iadtwp(ncor)-ilast4
c
      call iosys('read real "mcscf vector" from rwf',nob*nob,
     1           cr(icm),0,' ')
c
      if(incor.eq.1) then
         call getrab(cr(iabcx),cr(iscr),lbufso,'abcx',nao,nob,nao)
         if(nco.ne.0) then
            call getrab(cr(iabix),cr(iscr),lbufso,'abix',nco,nob,nao)
         end if
      endif
c
      call iosys('read character mcscf_type from rwf',-1,0,0,mcscf)
      call iosys('write character mcscf_type to mcscr',0,0,0,mcscf)
c
      if(noci.eq.0) then
         call iosys('read character "mcscf: ci used" from rwf',
     $              -1,0,0,citype)
         kx=0
         do 22133 i=1,navg
            call iosys('read real "mc root '//itoc(i)//'" from rwf',
     $                  ncsf,cr(ici+kx),0,' ')
         kx=kx+ncsf
22133    continue
      endif
c
      if (debug) then
          write(iout,*)' building ta_lagrangian'
          write(iout,*)' before gttamo: mo orbitals '
          call matout(cr(icm),nob,nob,nob,nob,iout)
      end if
c
c
      if(debug) then
         write(iout,*)'  incor for taint ',incor
      endif
c
      iindf=0
      do 10 i=1,ndf,mdf
         ldf=min(mdf,ndf-i+1)
c
         if(debug) then
            write(iout,*)' building ta_lagrangian: gttamo '
         endif
c
         call gttamo(cr(ncota),cr(naota),cr(iscr),cr(icm),cr(ibuf),
     1               nco,nao,nob,lbufso,ldf)
c
         if(debug) then
            write(iout,*)' building ta_lagrangian: taint '
         endif
c
         call taint(nob,nob,nco,nao,cr(icm),cr(ifpq),
     1              cr(ibuf),lbufso,cr(iabcx),cr(iabix),cr(isg),cr(itg),
     2              cr(iscr),cr(ir),cr(ncota),cr(naota),cr(ig1),cr(ig2),
     3              ldf,incor)
c
         if(debug) then
            write(iout,*)' building ta_lagrangian: derint'
            if(i.eq.1) then
               if(prttah) then
                  write(iout,*)' '
                  write(iout,*)' der-ta_g1 ndf=1 '
                  write(iout,17)(cr(ig1-1+ml),ml=1,nnao)
               endif
               if(prttag) then
                  write(iout,*)' '
                  write(iout,*)' der-ta_g2 ndf=1 '
                  write(iout,17)(cr(ig2-1+ml),ml=1,ntnao)
               endif
            endif
         endif
c
         call derint(cr(ig1),cr(ig2),cr(ibuf),lbufso,nao,ldf)
c
         if(debug) then
            if(i.eq.1) then
               write(iout,*)' '
               write(iout,*)' der-ta-ao_g1 ndf=1 '
               write(iout,17)(cr(ig1-1+ml),ml=1,nnao)
               write(iout,*)' der-ta-ao_g2 ndf=1 '
               write(iout,17)(cr(ig2-1+ml),ml=1,ntnao)
            endif
         endif
c
         ik1=0
         ik2=0
         ikci=0
c
         do 5 j=1,ldf
c
            call iosys('write real "mo one-electron integrals" '
     1                 //' to mcscr',nnao,cr(ig1+ik1),0,' ')
            call iosys('write real "mo two-electron integrals" '
     1                 //' to mcscr',ntnao,cr(ig2+ik2),0,' ')
c
            iindf=iindf+1
            kx=0
            do 4 k=1,navg
               junk=0
               if (citype.eq.'m902') then
                  call mn820(cr(ilast4),cr(ilast4),mcore,'mcscr',
     $                    'mcscr','mcscf','ci')
                  call mn902(cr(ilast4),cr(ilast4),mcore,'h times c',
     $                       cr(ici+kx),cr(itaci+ikci+kx),'mcscr',
     $                       'mcscr','ci',junk,'mc')
               else if (citype.eq.'m903') then
                  call mn903(cr(ilast4),cr(ilast4),mcore,'h times c',
     $                       cr(ici+kx),cr(itaci+ikci+kx),'mcscr',
     $                       'mcscr','ci',junk,'mc')
               end if
c
               ena=sdot(ncsf,cr(ici+kx),1,cr(itaci+ikci+kx),1)
               wtaao(k)=ena
c   
c              write(iout,*)' der. energy ndf ',ena,iindf
c
               call scopy(ncsf,cr(itaci+ikci+kx),1,cr(itbci+ikci+kx),1)
c
c              call sscal(ncsf,wt(k,1),cr(itaci+ikci+kx),1)
c
               kx=kx+ncsf
 4          continue
c
            call projct(cr(itaci+ikci),cr(ici),ncsfm,nwtm,nspinm)
            kx=0
            do 41 k=1,navg 
               call sscal(ncsf,wt(k,1),cr(itaci+ikci+kx),1)
               kx=kx+ncsf
 41         continue

c
c           general case derivative-ci gradient in the ci-vector basis
c
            if(lmixt.ne.0) then
               kx=nlcsf-1
               do 44 kk=1,nspinm
                  lmix0=lmixm(kk)
                  ncsf0=ncsfm(kk)
                  do 43 ll=1,lmix0
                     icc=iptci(1,ll,kk)
                     jcc=iptci(2,ll,kk)
                     ipt=(icc-1)*ncsf0
                     jpt=(jcc-1)*ncsf0
                     ena=sdot(ncsf0,cr(ici+jpt),1,cr(itbci+ikci+ipt),1)
c 
                     kx=kx+1
                     cr(itaci+ikci+kx)=(wt(icc,kk)-wt(jcc,kk))*ena
 43               continue
 44            continue
            end if
c
c..samc
c           call saxpy(ncsf,-ena,cr(ici),1,cr(itaci+ikci),1)
c
            ik1=ik1+nnao
            ik2=ik2+ntnao
            ikci=ikci+lncsf
  5      continue
c
         if(noci.eq.0) then
            call iosys('write real der_cigrad to mcscr'
     1                 //' without rewinding',ldf*lncsf,cr(itaci),0,' ')
            call iosys('write real der_civec to mcscr'
     1                 //' without rewinding',ldf*lncsf,cr(itbci),0,' ')
         end if
c..samc
c
  10  continue
c
      if(noci.eq.0) then 
c
c        move active integrals to mcscr for solution to lin. eqns.
c
         call iosys('read real "mo one-electron integrals" from rwf',
     1               nnao,cr(ig1),0,' ')
         call iosys('write real "mo one-electron integrals" to mcscr',
     1               nnao,cr(ig1),0,' ')
         call iosys('read real "mo two-electron integrals" from rwf',
     1               ntnao,cr(ig2),0,' ')
         call iosys('write real "mo two-electron integrals" to mcscr',
     1               ntnao,cr(ig2),0,' ')
      end if
c
      
      ncore=iadtwp(maxint)
c
      noc=nco+nao
      nnoc=noc*(noc+1)/2
      nobnob=nob*nob
      iscr=iabcx
      ilag=iscr+nobnob
      iymat=ilag+noc*nob
c
      lblock=noc*nob+nob*nob+nmix
c
      ilast=iymat+nobnob
      needm=ilast+lblock
      mxread=nobnob
c
c     test for minimum memory requirements
c
      if(needm.gt.ncore) then
         call lnkerr(' m1003: min. storage for ymatrix not avail. ')
      endif
c
c     does the mc_ao_ymatrix fits in core
c
      incor=1
c
      lnymat=nnoc*nobnob
c
      ilast=iymat+lnymat
      nleft=ncore-ilast
      mdf=nleft/lblock
      mdf=min(mdf,ndf)
      mdf=max(mdf,1)
c
      need0=ilast+mdf*lblock
      if(need0.lt.ncore.and.outym.eq.0)goto 90002
c
c     buffer-in the ymatrix
c
      mxread=(lnbuf/nobnob)*nobnob
      mxread=max(mxread,nobnob)
      mxread=min(mxread,lnymat)
c
      if(outym.ne.0) then
         mxread=2*nobnob
      endif
c
      ilast=iymat+mxread
      left=ncore-ilast
      mdf=left/lblock
      mdf=min(mdf,ndf)
c
      if(mdf.lt.1) then
         left=ncore-(iymat+lblock)
         mdf=1
         mxread=(left/nobnob)*nobnob
         ilast=iymat+mxread
         if(mxread.lt.0) then
            call lnkerr(' m1003: bug in core allocation ')
         endif
      endif
c
c
      incor=0
c
90002 continue
c
      if(debug) then
         write(iout,*)' m1003:incor and incorh for y-matrix ',
     $                incor,incorh
         if(mdf.lt.ndf) then
            write(iout,*)' *** multi-pass option invoked at 20 *** '
            write(iout,*)'  mdf  ndf ',mdf,ndf
         end if
      endif
c
      ita=ilast
      italag=ita+mdf*nob*nob
      igrad=italag+mdf*noc*nob
      ilast5=igrad+mdf*nmix
c
      if(ilast5.gt.ilast4) then
         call getscm(ilast5,cr,ncor,'ymatrix',0)
         ncore=iadtwp(ncor)-ilast5
      end if
c
      if(debug) then
          write(iout,*)' reading lagrangian and ymatrix '
      endif
c
      call iosys('read real mcscf_ao_lagrangian from rwf',
     1            noc*nob,cr(ilag),0,' ')
      if(incor.ne.0.and.incorh.ne.0) then
         call iosys('read real mc_ao_ymatrix from rwf',
     1               nnoc*nobnob,cr(iymat),0,' ')
      endif
      call iosys('read integer mcscf_mixhes from rwf',
     1            nocnob,icr(imixh),0,' ')
c
      if(debug) then
          write(iout,*)'  building ta-ci rhs '
      endif
c
      call iosys('rewind mo_der_overlap on dints',0,0,0,' ')
c
      do 20 i=1,ndf,mdf
         ldf=min(mdf,ndf-i+1)
c
         call gttaao(cr(ita),cr(iscr),cr(icm),nob,ldf)
c
         if(incor.ne.0.and.incorh.ne.0) then
c
            call xtalag(cr(italag),cr(ita),cr(iymat),cr(icm),cr(ilag),
     1                  cr(iscr),nco,nao,noc,nob,ldf)
c
         else
c
            call ytalag(cr(italag),cr(ita),cr(iymat),cr(icm),cr(ilag),
     1                  cr(iscr),nco,nao,noc,nob,ldf,mxread,incorh)
c
         endif
c
         call iosys('write real mo_ta_lag to mcscr without rewinding',
     1               ldf*nocnob,cr(italag),0,' ')
c
         if(debug) then
            if(i.eq.1) then
               write(iout,*)' ta-lag. for ndf=1 '
               call matout(cr(italag),nob,noc,nob,noc,iout)
            endif
         endif
c
         call derlag(cr(italag),cr(iscr),cr(igrad),cr(ibuf),lbufso,
     1               icr(imixh),nco,nao,nob,nmix,ldf)
c
         call iosys('write real der_grad to mcscr without rewinding',
     1               ldf*nmix,cr(igrad),0,' ')
c
  20  continue
c

      nobb(1)=nob
      ncob(1)=nco
      naob(1)=nao
      nfob(1)=0
      nsym=1
c..samc
      mdim=nmix+navg*ncsf+lmixt
c..samc
c
c..note
c     lncsf=navg*ncsf+lmixt
c..note

      nograd=intkey(ops,'cpmc=nograd',1,' ')
      if(debug) then
         write(iout,130) nograd
      endif
c
      kbuf=100000
c..samc
      kbuf=max(kbuf,lncsf,nmix)
      if(noci.eq.0) then
         mdf1=kbuf/lncsf
      else
         mdf1=ndf
      end if
c..samc
      mdf2=kbuf/nmix
      mdf1=min(ndf,mdf1)
      mdf2=min(ndf,mdf2)
      mdf=min(mdf1,mdf2)
c
c
      if(debug) then
         if(mdf.lt.ndf) then
             maxie=max(ncsf,nmix)
             write(iout,*)' *** multi-pass option invoked at 1000 *** '
             write(iout,*)'  nleft  need ',mdf*maxie,ndf*maxie
             write(iout,*)'  mdf  ndf ',mdf,ndf
         end if
      endif
c
c
      ilocs=wpadti(ibuf)
      ilocp=ilocs+nsym+1
      ilen=ilocp+noc
      imixo=ilen+noc
      igci=iadtwp(imixo+nocnob)
      igo=igci+mdf*(ncsf*navg+lmixt)
      ibci=igo+mdf*nmix
      ibo=ibci+mdf*(ncsf*navg+lmixt)
      ib=ibo+mdf*mdim
      ineed=ib+mdim
c
      call getscm(ineed,cr,ncor,'cpledr',0)
      ncore=iadtwp(ncor)-ineed
c
      call iosys('create real cpmc_soln on mcscr',ndf*mdim,0,0,' ')
c
      if(noci.eq.0) then
         call iosys('rewind der_cigrad on mcscr',0,0,0,' ')
      end if
c
      call iosys('rewind der_grad on mcscr',0,0,0,' ')
c
      call iosys('read integer mc_locsym from rwf',
     1            nsym+1,icr(ilocs),0,' ')
      call iosys('read integer mc_locp from rwf',
     1            noc,icr(ilocp),0,' ')
      call iosys('read integer mc_len from rwf',
     1            noc,icr(ilen),0,' ')
      call iosys('read integer mc_mixo from rwf',
     1            nocnob,icr(imixo),0,' ')
c
      kndf=0
      do 1000 i=1,ndf,mdf
c
         ldf=min(mdf,ndf-i+1)
c
         if(noci.eq.0) then
            call iosys('read real der_cigrad from mcscr'
     1                 //' without rewinding',ldf*lncsf,cr(igci),0,' ')
         end if
         call iosys('read real der_grad from mcscr without rewinding',
     1               ldf*nmix,cr(igo),0,' ')
c
         ikbo=0
         ikci=0
         iko=0
         if(debug) then
            write(iout,*)'  solve linear equations '
         endif
         do 1010 j=1,ldf
c
            kndf=kndf+1
            if(debug) then
               write(iout,*)' '
               write(iout,*)' perturbation ',kndf
            endif
c
            call scopy(nmix,cr(igo+iko),1,cr(ib),1)
            if(noci.eq.0) then
               call scopy(lncsf,cr(igci+ikci),1,cr(ib+nmix),1)
            end if
c
c
            if(prtrhs) then
               write(iout,*)' '
               write(iout,*)' cpmcscf-rhs  ndf= ',i-1+j
               write(iout,17)(cr(ib-1+mm),mm=1,mdim)
17             format(5(1x,f12.8))
            end if
c
            call cpledr(nsym,nobb,nobb,nfob,ncob,naob,cr(icm),
     #           nmix,ncsf,mdim,noci,istopj,stopj,cr(ib),sqcdf,
     #           cr(ineed),cr(ineed),ncore,
     #           lbufso,nograd,1,icr(ilocs),icr(ilocp),icr(ilen),
     #           icr(imixo),prtena,mcroot,
     #           nwtm,wt,ncsfm,nspinm,navg,iptci,lmixm,lmixt,nlcsf)
c
            call scopy(mdim,cr(ib),1,cr(ibo+ikbo),1)
c
            ikbo=ikbo+mdim
            ikci=ikci+lncsf
            iko=iko+nmix
c
 1010    continue
c
         call iosys('write real cpmc_soln to mcscr without rewinding',
     1              ldf*mdim,cr(ibo),0,' ')
c
 1000 continue
c
      iforce=ibuf
      ita=iforce+ndf*ndf
      ium=ita+ndf*nobnob
      igama=ium+ndf*nobnob
      ibvec=igama+ndf*mdim
      ilagt=ibvec+ndf*lncsf
      ilaga=ilagt+ndf*nocnob
      imclag=ilaga+ndf*nocnob
      itfor=imclag+nocnob
      ineed=itfor+ndf*ndf
      call getscm(ineed,cr,ncor,'mcforce',0)
cccc
      len=ndf*nobnob
      call rzero(cr(ium),len)
cccc
      call iosys('read real cpmc_soln from mcscr',
     1            ndf*mdim,cr(igama),0,' ')
      call iosys('read real mo_der_overlap from dints',
     1            ndf*nobnob,cr(ita),0,' ')
c
      if(noci.eq.0) then
         call iosys('read real der_civec from mcscr',
     1            ndf*lncsf,cr(ibvec),0,' ')
      end if
c
      call iosys('read real mo_ta_lag from mcscr',
     1            ndf*nocnob,cr(ilagt),0,' ')
      call iosys('read real mo_der_lagrangian from rwf',
     1            ndf*nocnob,cr(ilaga),0,' ')
      call iosys('read real mcscf_mo_lagrangian from rwf',
     1            nocnob,cr(imclag),0,' ')
cccc
      igam=igama
      iu=ium
      it=ita
      ilt=ilagt
      ila=ilaga
      do 7000 i=1,ndf
c
         call setta(cr(it),nob)
         call umat(cr(iu),cr(igam),icr(imixh),nco,nao,nob)
         if(debug) then
            write(iout,*)' ilagt before fixlag '
            call matout(cr(ilt),nob,noc,nob,noc,6)
            write(iout,*)' ilaga before fixlag '
            call matout(cr(ila),nob,noc,nob,noc,6)
         endif
         call fixlag(cr(ilt),cr(ila),nocnob)
c
         iu=iu+nobnob
         it=it+nobnob
         igam=igam+mdim
         ilt=ilt+nocnob
         ila=ila+nocnob
 7000 continue
c
c------------------------------------------------------------c
c     read integral contribution to the force constant matrix
c------------------------------------------------------------c
c
      if(forc) then
         itrfor=ndf*(ndf+1)/2
         if(debug) then
            call rzero(cr(itfor),itrfor)
         else
            call iosys('read real "integral force constants" from rwf',
     #                  itrfor,cr(itfor),0,' ')
         endif
c
         if(prtcp) then
            write(iout,*)' '
            write (iout,*) ' mcscf integral force constants '
            call print(cr(itfor),itrfor,ndf,iout)
         end if
c
         iu=ium
         iga=igama+nmix
         ibv=ibvec
         ilt=ilagt
         ila=ilaga
         it=ita
cc
         numdf=0
cc
         do 8010 i=1,ndf
cc
            ju=ium
            jga=igama+nmix
            jbv=ibvec
            jlt=ilagt
            jla=ilaga
            jt=ita
cc
cc
            do 8000 j=1,i
c
               call force(cr(iu),cr(iga),cr(ibv),
     1                    cr(ilt),cr(ila),cr(it),
     2                    cr(ju),cr(jga),cr(jbv),
     3                    cr(jlt),cr(jla),cr(jt),
     4                    cr(imclag), noc,nob, ncsf, xforc, prtall)
c
               numdf=numdf+1
               cr(iforce+numdf-1)=xforc
c
c
               ju=ju+nobnob
               jga=jga+mdim
               jbv=jbv+ncsf
               jlt=jlt+nocnob
               jla=jla+nocnob
               jt=jt+nobnob
c
 8000       continue
cc
            iu=iu+nobnob
            iga=iga+mdim
            ibv=ibv+ncsf
            ilt=ilt+nocnob
            ila=ila+nocnob
            it=it+nobnob
cc
 8010    continue
cc
         if(prtcp) then
            write(iout,*)'  '
            write(iout,*)
     $          'cpmcscf contributions to the force constants '
            call print(cr(iforce),itrfor,ndf,iout)
         end if
c
c        ----- add to cpmcscf contributions to the integral contributions
c
         call vadd(cr(iforce),cr(iforce),cr(itfor),itrfor)
c
         if(prtfor) then
            write(iout,150)
            call print(cr(iforce),itrfor,ndf,iout)
         end if
c
          call iosys('write real "cartesian second derivatives" to rwf',
     $               itrfor,cr(iforce),0,' ')
      endif
c
c     ----- form cpmcscf u-matrix -----
c       u(i,j) contains cpmcscf solutions  u(i,j)=-u(j,i)
c       t(i,j) contains the overlap contribution
      iu0=ium-1
      it0=ita-1
      do 9010 i=1,ndf*nobnob
         cr(iu0+i)=cr(iu0+i)+cr(it0+i)
 9010 continue
c
      if(prtum) then
         ix=ium
         do 9011 i=1,ndf
            write(iout,*)' '
            write(iout,*)'   u-matrix  ndf = ',i
            call matout(cr(iu+ix),nob,nob,nob,nob,iout)
            ix=ix+nobnob
 9011    continue
      endif
c
      write(iout,9012)
 9012 format(5x,'cpmcscf solutions stored on rwf ')
c
      call iosys('write real "cphf solutions" to rwf',
     #            ndf*nobnob,cr(ium),0,' ')
c
      call iosys('write character "print flag" to rwf',0,0,0,prtsav)
c
c     ----- exit this link -----
      call iosys('close dints',0,0,0,' ')
      call chainx(0)
c
      stop
      end
