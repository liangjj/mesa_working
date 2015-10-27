*deck @(#)outtrn.f	5.1  11/6/94
      subroutine outtrn(z,a,nwint,nco,nao,nob,nbf,nblock,
     $                  niter,ecore,enuc,ncsf,lsect,opnscf)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      871116   (yymmdd)
c   16 november 1987   bhl      lanl
c   modified to skip over transformations if nao=0
c   changes denoted by c..bhl..close shell scf
c    and to call getscm before call to ktran
c***revision date      871216   (yymmdd)
c   16 december 1987   bhl      brl	
c   code to transform the one-electron integrals to the
c   active-mo basis and write them to rwf
c   changes denoted by c..bhl..h_int
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)outtrn.f	5.1   11/6/94
c
c***purpose            to transform the one- and two-electron
c                      integrals form the atomic orbital to
c                      the molecular orbital basis.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit integer (a-z)
c
      real*8 eonel
      real*8 etwoel
      real*8 z,enuc,ecore,sdot
      integer a(1)
c
c...
c...  equivalence (z,a)
c...  integer a(nwint)
c...
      character*3 ians
      dimension z(2),ihd(20)
c
c
      common /io/ inp,iout
c
c
c     data ncorsx/40000/
c     save ncorsx
c
c
      nocc=nco+nao
c
c     cdir$ fastmd
c
      need=27
      nsym=1
      numij=(nao*(nao+1))/2
c
c
c     ----- now get the scf vector -----
c
c
      hmo=iadtwp(need)
      c=hmo+numij
      tempc=c+nbf*nbf
      eigval=c
      i10=wpadti(tempc+nbf**2)
      need=i10+200
      call getscm(need,a,maxcor,'mctrans8',0)
c
c     skip set 1
c
c
      ntot=nbf*nbf
      nbf22=ntot
c
      call iosys('read real mo_new from mcscr',ntot,z(c),0,' ')
c
c      if(niter.eq.1) then
c       write(iout,90011)
c90011  format(//'  input orbitals',/)
c       call vecout(z(c),nbf,nob)
c      endif
c
c     ----- core
c           density matrix, if necessary
c
      nnp=(nbf+1)*nbf/2
      nbf2=nnp
c
      cdens=tempc+nbf**2
      need=wpadti(cdens+nnp)
      call getscm(need,a,maxcor,'mctrans9',0)
c
c
c     ----- open the integral file and check the orthonormality of
c                   the orbitals
cc
      ihd(1)=nco
      ihd(2)=nao
      ihd(3)=nob
      ihd(4)=nbf
cc
      ihd(5)=nblock
c     ihd(6)=nbufj
c     ihd(7)=nbufk
      ihd(6)=nbf22
      ihd(7)=nbf22
cc
cc
      s=iadtwp(need)
      h=s+nnp
      v=h+nnp
      need=wpadti(v+nnp)
      call getscm(need,a,maxcor,'mctransa',0)
c
      call iosys('read real "overlap integrals" from rwf',
     $           nnp,z(s),0,' ')
      call iosys('read real "kinetic integrals" from rwf',
     $            nnp,z(h),0,' ')
      call iosys('read real "potential integrals" from rwf',
     $            nnp,z(v),0,' ')
c
      call vadd(z(h),z(h),z(v),nnp)
c
      t1=iadtwp(need)
      t2=t1+nbf**2
      t3=t2+nbf**2
      need=wpadti(t3+nbf**2)
      call getscm(need,a,maxcor,'mctransb',0)
c
c     schmidt orthogonalize the vectors
c
      call schvec(z(c),z(s),z(t1),z(t2),z(t3),nbf,nnp)
c
      call chknrm(z(c),z(s),z(t1),z(t2),z(t3),nbf,nnp)
c
      call iosys('write real mo_old on mcscr',ntot,z(c),0,' ')
c
c
c..bhl..h_int
c
c     transform one-electron integrals
c
      if(nao.ne.0) then
c
c     ----- step over core orbitals
c
       iof=nco*nbf
       call trn1e(z(h),z(c+iof),nbf,nao,nnp,numij,z(t1),z(t2),z(hmo))
c
       call iosys('write real mc_h_int to rwf',numij,z(hmo),0,' ')
c
      endif
c
c..bhl..h_int
c
c
      if(nco.ne.0) then
         call corden(z(c),z(cdens),nbf,nco,nnp)
      else
         call iosys('write real mc_core_fock to mcscr',nbf2,z(h),0,' ')
      endif
c
c     ----- form the frozen-core fock-matrix -----
c
      ecore=0.0d+00
      if (nco.gt.0) then
         cf=iadtwp(need)
         values=cf+nnp
         need=wpadti(values+min(60000,nnp**2))
         call getscm(need,a,maxcor,'mctransc',0)
         ntriang=(iadtwp(need)-values)/nnp
c
c
         call fock(z(values),z(cdens),z(cf),nnp,itap44,nbf,z(t1),z(t2),
     $        ntriang)
c
c     ----- form the frozen core energy -----
c
cps         call vadd(z(cf),z(cf),z(h),nnp)
         call trtosq(z(t1),z(cdens),nbf,nnp)
         call trtosq(z(t2),z(cf),nbf,nnp)
         call trtosq(z(t3),z(h),nbf,nnp)
         eonel=sdot(nbf**2,z(t1),1,z(t3),1)
         etwoel=0.5d+00*sdot(nbf**2,z(t1),1,z(t2),1)
         ecore=eonel+etwoel
         call iosys('write real "mcscf core 1e energy" to rwf',
     $        1,eonel,0,' ')
         call iosys('write real "mcscf core 2e energy" to rwf',
     $        1,etwoel,0,' ')
         call iosys('write real "mcscf one-electron energy" to rwf',
     $        1,eonel,0,' ')
         call iosys('write real "mcscf two-electron energy" to rwf',
     $        1,etwoel,0,' ')
         call iosys('write real "frozen core energy" to rwf',
     $        1,ecore,0,' ')
cps         call vmove(z(h),z(cf),nnp)
         call vadd(z(h),z(h),z(cf),nnp)
c
c     ----- store the core fock matrix and core density matrix ----
c
         nbf2=(nbf*(nbf+1))/2
         call iosys('write real mc_core_fock to mcscr',
     $               nbf2,z(h),0,' ')
         call iosys('write real "mcscf ao core density" to rwf',
     $               nbf2,z(cdens),0,' ')
c
      else
c
         eonel=0.d+00
         etwoel=0.d+00
c
         call iosys('write real "mcscf core 1e energy" to rwf',
     $        1,eonel,0,' ')
         call iosys('write real "mcscf core 2e energy" to rwf',
     $        1,etwoel,0,' ')
         call iosys('write real "frozen core energy" to rwf',
     $        1,ecore,0,' ')
      end if
c
c     ----- get nuclear repulsion energy from tape10 -----
c
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,enuc,0,' ')
c
c     ..... save one electron ints ( fock operator ) on tape16
c
      nbf2=(nbf*(nbf+1))/2
c
c     ----- transform the one-electron integrals -----
c
c
c..bhl..closed shell scf
c
      if(nao.ne.0) then

c     ----- step over core orbitals
c
       iof=nco*nbf
       call trn1e(z(h),z(c+iof),nbf,nao,nnp,numij,z(t1),z(t2),z(hmo))
c
       call iosys('does mc_tint1 exist on mcscr',0,0,0,ians)
c
       if(ians.eq.'no') then
          call iosys('create real mc_tint1 on mcscr',numij,0,0,' ')
       endif
c
       call iosys('write real mc_tint1 to mcscr',numij,z(hmo),0,' ')
       call iosys('write real "mo one-electron integrals" to mcscr',
     $      numij,z(hmo),0,' ')
c
      endif
c
c     ----- reallocate core for the two-electron transformation -----
c
      nnb=(nocc*(nocc+1))/2
      noc2=nnb
      nao2=(nao*(nao+1))/2
      nao4=nao2*nao2
      nocnbf=nocc*nbf
c
      t1=c+nbf*nbf
      t2=t1+nbf**2
      values=t2+nbf*nbf
c
      ntriang=min(nnp,60000/nnp)
      ntriang=max(ntriang,nao2**2/nnp+1)
c
      need=wpadti(values+ntriang*nnp)
c
      call getscm(0,a,maxcor,'mctransd',0)
c
      nwwp=iadtwp(maxcor)
c
c..bhl.big.start
c
      call iosys('does mc_j exist on mcscr',0,0,0,ians)
c
      lenjda=noc2*nbf2
      lenkda=nocnbf*nbf2
      lenkkda=noc2*nbf22
      lenjkda=lenjda+lenkda+2
      if(ians.eq.'no') then
         mtot=lenjkda
         mstot=lenkkda+2
         call iosys('create real mc_j on mcscr',mtot,0,0,' ')
         call iosys('create real mc_k on mcscr',mstot,0,0,' ')
      endif
c
c     ihd(8)=jnij
c     ihd(9)=jmij
c     ihd(10)=jnkl
c     ihd(11)=jmkl
c     ihd(12)=kknij
c     ihd(13)=kkmij
c     ihd(14)=kknkl
c     ihd(15)=kkmkl
      ihd(16)=0
      ihd(17)=0
      ihd(18)=1
c
cc
c  ..... put header record of j - k integral file
cc
      call iosys('does mc_jk_hdr exist on mcscr',0,0,0,ians)
      if(ians.eq.'no') then
         call iosys('create integer mc_jk_hdr on mcscr',20,0,0,' ')
      endif
c
      call iosys('write integer mc_jk_hdr to mcscr',20,ihd,0,' ')
c
c
c ---- core allocation for jtran
c
      nbufj=nbf22
      nbufk=nbf22
c
      last=iadtwp(need)
c
      ncors=nwwp-nbf22-3*nbufj-noc2-3*nbufk-nocnbf
      ncors=ncors-nbf*nao*nao2-nao2*nao*nao-last
c
      ncors=min(ncors,lenjkda)
c552      if(ncors.lt.0) call lnkerr(' ncors lt 0 in jtran ')
c
      asort=last
      t3=asort+ncors
      xbj=t3+nbf*nbf
      labj=xbj+nbufj
      labjj=labj+nbufj
      jlab=labjj+nbufj
      xbk=jlab+noc2
      labk=xbk+nbufk
      labkk=labk+nbufk
      klab=labkk+nbufk
      valt=klab+nocnbf
      xval=valt+nbf*nao*nao2
      need=wpadti(xval+nao2*nao*nao)
      call getscm(need,a,maxcor,'mctranse',0)
c
c
c
cbigmc
c
c..bhl..closed shell scf
c
      if(nao.ne.0) then
       call reordc(z(c),z(t1),nco,nao,nbf)
      endif
c
c
      call jtran(z(values),nnp,ntriang,z(c),nbf,nco,nao,nocc,z(t1),
     $     z(t2),z(t3),z(asort),ncors, z(valt),z(xval), lsortj,
     $     z(xbj),z(labj),z(labjj),z(jlab),nbufj,
     $     z(xbk),z(labk),z(labkk),z(klab),nbufk)
c
c       write(iout,*)'  jtran completed '
c
c
c..bhl..closed shell scf
c
      if(nao.ne.0) then
c
       ntot=numij*numij
       call iosys('does mc_tint2 exist on mcscr',0,0,0,ians)
       if(ians.eq.'no') then
          call iosys('create real mc_tint2 on mcscr',ntot,0,0,' ')
       endif
c
       call iosys('write real mc_tint2 to mcscr',ntot,z(values),0,' ')
       call iosys('write real "mo two-electron integrals" to mcscr',
     $      ntot,z(values),0,' ')
c
      endif
c
c
c ---- core allocation for final transform
c
c
      xbk=xbj
      labk=xbk+nbufk
      labkk=labk+nbufk
      klab=labkk+nbufk
      need=wpadti(klab)
c
      call getscm(need,a,maxcor,'ktran',0)
c
      call ktran(z(values),nnp,ntriang,z(c),nbf,nco,nao,nocc,z(t1),
     $     z(asort),ncors,  lsortj,
     $     z(xbk),z(labk),z(labkk),z(klab),nbufk)
c
c
      return
      end
