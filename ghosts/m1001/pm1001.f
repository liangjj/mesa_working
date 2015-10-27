*deck %W%  %G%
      subroutine pm1001(cr,icr)
c***begin prologue     pm1001
c***date written       860501  (yymmdd)  
c***revision date      910618  (yymmdd)
c  18 june     1991    rlm at lanl
c      creating separate disk files for derivative integrals,etc.
c   4 december 1986    pws at lanl
c      1. changing 'namint' to character variable.
c      2. changing iosys opens to character unit names
c      3. adding family-member-size to mcscr.
c***keywords           
c***author             
c***source             %W%   %G%
c***purpose            solves coupled perturbed mcscf equations. 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       pm1001
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
      if(logkey(ops,'cpmc=force',.false.,' ').or.
     $   logkey(ops,'opt=force',.false.,' ')) forc=.true.
c
      if(prtall) prtcp=.true.
      if(prtcp)  prtfor=.true.
c
      write(iout,*)' '
      write(iout,*)' m1001:    coupled-perturbed mcscf  '
      write(iout,*)'           byron lengsfield  brl dec 1987 '
      write(iout,*)' '
      write(iout,*)' maximum number of cycles  ',istopj
      write(iout,*)' convergence criteria      ',stopj
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
c..bhl
c      write(iout,19001) nbf,nco,nao,nob,lbufso,ncsf,nmix,natoms
c19001 format(' nbf nco nao nob lbufso ncsf nmix natoms',/,10(1x,i6))
c..bhl
c     allocate core
      nroots=1
c
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
      ibuf=ici+ncsf
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
        incor=0
        iabcx=ibuf+lbufso
        iabix=iabcx+1
        iscr=iabix+1
        ifpq=iscr+nob*nob
        ir=ifpq+nob*nob
        isg=ir+nob*nob
        ilast=isg+nnao*naonao
        ilast=max(ilast,iscr+lbufso)
        if(ilast.gt.ncore) then
           call lnkerr(' m1001: ilast gt ncore ')
         endif
      end if
c
      nleft=ncore-ilast-1
      lblock=nob*noc+nnao+ntnao+2*ncsf+nnao*naonao
      igint=nnao+ntnao
      mdf=nleft/lblock
      mdf=min(mdf,ndf)
c
      if(mdf.lt.ndf) then
         write(iout,*)' ***** multi-pass option invoked ***** '
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
      call iosys('create real der_cigrad on mcscr',ndf*ncsf,0,0,' ')
      call iosys('create real der_civec on mcscr',ndf*ncsf,0,0,' ')
      call iosys('create real mo_ta_lag on mcscr',
     1            ndf*nocnob,0,0,' ')
      call iosys('create real der_grad on mcscr',ndf*nmix,0,0,' ')
c
c     ----- open integral units -----
      call iosys('read character "derivative integral filename"'
     $            //' from rwf',0,0,0,dints)
      call iosys('open dints as old',0,0,0,dints)
c
      call iosys('read real mcscf_ao_core_fock from rwf',
     1            nobnob,cr(iscr),0,' ')
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
      itbci=itaci+mdf*ncsf
      itg=itbci+mdf*ncsf
      ilast4=itg+mdf*nnao*naonao
c
      ilast4=max(ilast4,iscr+lbufso)
      need4=ilast4+2*ncsf+nnao+ntnao
      call getscm(need4,cr,ncor,'taint',0)
      mcore=iadtwp(ncor)-ilast4
c
      call iosys('read real "mcscf vector" from rwf',nob*nob,
     1           cr(icm),0,' ')
c
      if(incor.eq.1) then
         call getrab(cr(iabcx),cr(iscr),lbufso,'abcx',nao,nob,nao)
         call getrab(cr(iabix),cr(iscr),lbufso,'abix',nco,nob,nao)
      endif
c
      call iosys('read character mcscf_type from rwf',-1,0,0,mcscf)
      call iosys('write character mcscf_type to mcscr',0,0,0,mcscf)
c
c..rlm see m1003
c     this probably won't work for scf, but will be replaced by m1003
c     soon
      noci=0
      if(logkey(ops,'mcscf=hf',.false.,' ')) then
         noci=1
         ncsf=0
      endif
c..rlm
      if(noci.eq.0) then
         call iosys('read character "mcscf: ci used" from rwf',
     $              -1,0,0,citype)
         call iosys('read real "mc root '//itoc(nroots)//'" from rwf',
     1               ncsf,cr(ici),0,' ')
      endif
c..bhl
c      write(iout,*)' building ta_lagrangian'
c      endfile iout
c      backspace iout
c..bhl
      if (debug) then
          write(iout,*)' before gttamo: mo orbitals '
          call matout(cr(icm),nob,nob,nob,nob,iout)
      end if
c..bhl
c.debug
      write(iout,*)'  incor for taint ',incor
c.debug
      do 10 i=1,ndf,mdf
         ldf=min(mdf,ndf-i+1)
c..bhl
c        write(iout,*)' building ta_lagrangian: gttamo '
c        endfile iout
c        backspace iout
c..bhl
         call gttamo(cr(ncota),cr(naota),cr(iscr),cr(icm),cr(ibuf),
     1               nco,nao,nob,lbufso,ldf)
c..bhl
c        write(iout,*)' building ta_lagrangian: taint '
c        endfile iout
c        backspace iout
c..bhl
         call taint(nob,nob,nco,nao,cr(icm),cr(ifpq),
     1              cr(ibuf),lbufso,cr(iabcx),cr(iabix),cr(isg),cr(itg),
     2              cr(iscr),cr(ir),cr(ncota),cr(naota),cr(ig1),cr(ig2),
     3              ldf,incor)
c..bhl
c        write(iout,*)' building ta_lagrangian: derint'
c        endfile iout
c        backspace iout
c..bhl
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
c..bhl
         call derint(cr(ig1),cr(ig2),cr(ibuf),lbufso,nao,ldf)
c..bhl
c        if(i.eq.1) then
c           write(iout,*)' '
c           write(iout,*)' der-ta-ao_g1 ndf=1 '
c           write(iout,17)(cr(ig1-1+ml),ml=1,nnao)
c           write(iout,*)' der-ta-ao_g2 ndf=1 '
c           write(iout,17)(cr(ig2-1+ml),ml=1,ntnao)
c        endif
c..bhl
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
            zero=0.d0
            junk=0
            if (citype.eq.'m902') then
               call mn820(cr(ilast4),cr(ilast4),mcore,'mcscr','mcscr',
     $                    'mcscf','ci')
               call mn902(cr(ilast4),cr(ilast4),mcore,'h times c',
     $                    cr(ici),cr(itaci+ikci),'mcscr','mcscr',
     $                    'ci',junk,'mc')
            else if (citype.eq.'m903') then
               call mn903(cr(ilast4),cr(ilast4),mcore,'h times c',
     $                    cr(ici),cr(itaci+ikci),'mcscr','mcscr',
     $                    'ci',junk,'mc')
            end if
c
            ena=sdot(ncsf,cr(ici),1,cr(itaci+ikci),1)
c
c..bhl
            if(prtena) then
               write(iout,70002) i-1+j,ena
70002          format(' degree of freedom ',i5,'  e(ta+ao)=  ',f12.8)
            endif
c..bhl
            call scopy(ncsf,cr(itaci+ikci),1,cr(itbci+ikci),1)
c
            call saxpy(ncsf,-ena,cr(ici),1,cr(itaci+ikci),1)
c
            ik1=ik1+nnao
            ik2=ik2+ntnao
            ikci=ikci+ncsf
c
  5      continue
c
         call iosys('write real der_cigrad to mcscr without rewinding',
     1               ldf*ncsf,cr(itaci),0,' ')
         call iosys('write real der_civec to mcscr without rewinding',
     1               ldf*ncsf,cr(itbci),0,' ')
c
  10  continue
c
c  move active integrals to mcscr for solution to lin. eqns.
c
      call iosys('read real "mo one-electron integrals" from rwf',
     $            nnao,cr(ig1),0,' ')
      call iosys('write real "mo one-electron integrals" to mcscr',
     $            nnao,cr(ig1),0,' ')
      call iosys('read real "mo two-electron integrals" from rwf',
     $            ntnao,cr(ig2),0,' ')
      call iosys('write real "mo two-electron integrals" to mcscr',
     $            ntnao,cr(ig2),0,' ')
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
c  test for minimum memory requirements
c
      if(needm.gt.ncore) then
         call lnkerr(' m1001: min. storage for ymatrix not avail. ')
      endif
c
c  does the mc_ao_ymatrix fits in core
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
c   buffer-in the ymatrix
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
            call lnkerr(' m1001: bug in core allocation ')
         endif
      endif
c
c
      incor=0
c
90002 continue
c
      if(debug) then
         write(iout,*)' m1001:incor and incorh for y-matrix ',
     $                incor,incorh
         if(mdf.lt.ndf) then
            write(iout,*)
     $            ' ***** multi-pass option invoked at 20 *****'
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
         endfile iout
         backspace iout
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
         endfile iout
         backspace iout
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
c..bhl
c        if(i.eq.1) then
c           write(iout,*)' ta-lag. for ndf=1 '
c           call matout(cr(italag),nob,noc,nob,noc,iout)
c        endif
c..bhl
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
      noci=0
      mdim=ncsf+nmix
      nograd=1
c
      kbuf=100000
      kbuf=max(kbuf,ncsf,nmix)
      mdf1=kbuf/ncsf
      mdf2=kbuf/nmix
      mdf1=min(ndf,mdf1)
      mdf2=min(ndf,mdf2)
      mdf=min(mdf1,mdf2)
c
      if(debug) then
         if(mdf.lt.ndf) then
             maxie=max(ncsf,nmix)
             write(iout,*)
     $          ' ***** multi-pass option invoked at 1000 ***** '
             write(iout,*)'  nleft  need ',mdf*maxie,ndf*maxie
             write(iout,*)'  mdf  ndf ',mdf,ndf
         end if
      endif
c
      ilocs=wpadti(ibuf)
      ilocp=ilocs+21
      ilen=ilocp+noc
      imixo=ilen+noc
      igci=iadtwp(imixo+nocnob)
      igo=igci+mdf*ncsf
      ibci=igo+mdf*nmix
      ibo=ibci+mdf*ncsf
      ib=ibo+mdf*mdim
      ineed=ib+mdim
c
      call getscm(ineed,cr,ncor,'cpledr',0)
      ncore=iadtwp(ncor)-ineed
c
      call iosys('create real cpmc_soln on mcscr',ndf*mdim,0,0,' ')
c
      call iosys('rewind der_cigrad on mcscr',0,0,0,' ')
      call iosys('rewind der_grad on mcscr',0,0,0,' ')
c
      call iosys('read integer mc_locsym from rwf',
     1            nsym+1,icr(ilocs),0,' ')
c    1            21,icr(ilocs),0,' ')
      call iosys('read integer mc_locp from rwf',
     1            noc,icr(ilocp),0,' ')
      call iosys('read integer mc_len from rwf',
     1            noc,icr(ilen),0,' ')
      call iosys('read integer mc_mixo from rwf',
     1            nocnob,icr(imixo),0,' ')
c
      do 1000 i=1,ndf,mdf
c
         ldf=min(mdf,ndf-i+1)
c
         call iosys('read real der_cigrad from mcscr without rewinding',
     1               ldf*ncsf,cr(igci),0,' ')
         call iosys('read real der_grad from mcscr without rewinding',
     1               ldf*nmix,cr(igo),0,' ')
c
         ikbo=0
         ikci=0
         iko=0
c
         if(debug) then
            write(iout,*)'  solve linear equations '
            endfile iout
            backspace iout
         endif
c
         do 1010 j=1,ldf
c
            call scopy(nmix,cr(igo+iko),1,cr(ib),1)
            call scopy(ncsf,cr(igci+ikci),1,cr(ib+nmix),1)
c
            if(prtrhs) then
               write(iout,*)' '
               write(iout,*)' cpmcscf-rhs  ndf= ',i-1+j
               write(iout,17)(cr(ib-1+mm),mm=1,mdim)
17             format(5(1x,f12.8))
            end if
c
            call cpledr(nsym,nobb,nobb,nfob,ncob,naob,cr(icm),
     #                  nmix,ncsf,mdim,noci,istopj,stopj,cr(ib),sqcdf,
     #                  cr(ineed),cr(ineed),ncore,
     #                  lbufso,nograd,1,icr(ilocs),icr(ilocp),icr(ilen),
     #                  icr(imixo),prtena)
c
            if(prtsol) then
               write(iout,*)' solution to lin.-eqn. for ndf= ',i-1+j
               write(iout,17)(cr(ib-1+im),im=1,mdim)
            endif
c
            call scopy(mdim,cr(ib),1,cr(ibo+ikbo),1)
c
            ikbo=ikbo+mdim
            ikci=ikci+ncsf
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
      ilagt=ibvec+ndf*ncsf
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
      call iosys('read real der_civec from mcscr',
     1            ndf*ncsf,cr(ibvec),0,' ')
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
c..bhl
c1       write(iout,*)' ilagt before fixlag '
c1       call matout(cr(ilt),nob,noc,nob,noc,6)
c1       write(iout,*)' ilaga before fixlag '
c1       call matout(cr(ila),nob,noc,nob,noc,6)
c..bhl
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
      if(.not.forc) goto 8011
c
      itrfor=ndf*(ndf+1)/2
      if(debug) then
         call rzero(cr(itfor),itrfor)
      else
         call iosys('read real "integral force constants" from rwf',
     #            itrfor,cr(itfor),0,' ')
      endif
c
      if(prtcp) then
         write(iout,*)' '
         write (iout,*) ' mcscf integral force constants '
         call print(cr(itfor),itrfor,ndf,iout)
      end if
c
cc
c  store cpmcscf solution on unit it23
c    .. first add overlap and variational contributions
cc
c      call addv(cr(iscr),cr(iu),cr(it),nobnob)
c      call swrit(it23,cr(iga),intowp(ncsf))
cc
cccc
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
     1                 cr(ilt),cr(ila),cr(it),
     2                 cr(ju),cr(jga),cr(jbv),
     3                 cr(jlt),cr(jla),cr(jt),
     4                 cr(imclag), noc,nob, ncsf, xforc, prtall)
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
 8000    continue
cc
         iu=iu+nobnob
         iga=iga+mdim
         ibv=ibv+ncsf
         ilt=ilt+nocnob
         ila=ila+nocnob
         it=it+nobnob
cc
 8010 continue
cc
      if(prtcp) then
         write(iout,*)'  '
         write(iout,*) ' cpmcscf contributions to the force constants '
         call print(cr(iforce),itrfor,ndf,iout)
      end if
c
c  add to cpmcscf contributions to the integral contributions
c
      call vadd(cr(iforce),cr(iforce),cr(itfor),itrfor)
c
      if(prtfor) then
         write(iout,*)'  '
         write(iout,*) ' mcscf force constants '
         call print(cr(iforce),itrfor,ndf,iout)
      end if
c
       call iosys('write real "cartesian second derivatives" to rwf',
     $            itrfor,cr(iforce),0,' ')
ccc
cc
ccc
 8011 continue
c
c  form cpmcscf u-matrix
c       u(i,j) contains cpmcscf solutions  u(i,j)=-u(j,i)
c       t(i,j) contains the overlap contribution
c
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
 9012 format(/,' cpmcscf solutions stored on rwf ',//)
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
