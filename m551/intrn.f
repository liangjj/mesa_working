*deck  @(#)intrn.f	1.3 7/30/91
      subroutine intrn(z,a,nwint,nco,nao,nob,nbf,nblock,
     $                 niter,ecore,enuc,ncsf,lsect,opnscf)
c
c***begin prologue     intrn
c***date written       871022   (yymmdd)
c***revision date      871116   (yymmdd)
c   16 november 1987   bhl lanl
c   modified to skip over transformations if nao=0
c   changes denoted by c..bhl..close shell scf
c***revision date      871216   (yymmdd)
c   16 december 1987   bhl      brl	
c   code to transform the one-electron integrals to the
c   active-mo basis and write them to rwf
c   changes denoted by c..bhl..h_int
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)intrn.f	1.1   11/30/90
c
c***purpose            to transform the one- and two-electron integrals
c                      from the atomic-orbital to the molecular orbital
c                      basis.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       intrn
c
      implicit integer (a-z)
c
      real*8 eonel
      real*8 etwoel
      real*8 z,enuc,ecore,sdot
      character*3 ians
      integer a(1)
c
      dimension z(2),ihd(20)
c
      common /io/ inp,iout
c...
c...  equivalence (cr,icr)
c...  integer a(nwint)
c...
c
      nocc=nco+nao
c
      need=27
      nsym=1
      numij=(nao*(nao+1))/2
c
c
c     ----- now get the scf vector -----
c
      ntot=nbf*nbf
      hmo=iadtwp(need)
      c=hmo+numij
      tempc=c+nbf*nbf
      eigval=c
      i10=wpadti(tempc+nbf**2)
      need=i10+200
      call getscm(need,a,maxcor,'intrn: mctrans1',0)
c
      call iosys('read real mo_new from mcscr',ntot,z(c),0,' ')
c
c
c     ----- core
c           density matrix, if necessary
c
      nnp=(nbf+1)*nbf/2
      cdens=tempc+nbf**2
      need=wpadti(cdens+nnp)
      call getscm(need,a,maxcor,'intrn: mctrans2',0)
c
c
c     ----- open the integral file and check the orthonormality of
c                   the orbitals
c
      lsect1=lsect+1
      nnpp=nbf*nbf
      npj=nblock/nnp
      nbufj=npj*nnp
      npk=nblock/nnpp
      nbufk=npk*nnpp
cc
      if(npk.le.0) then
         write(iout,*)'  buffer size for jk integrals is to small'
         write(iout,*)'  nnpp  nblock ',nnpp,nblock
         call lnkerr(' intrn: buffer size for jk integrals too small')
      endif
c
      s=iadtwp(need)
      h=s+nnp
      v=h+nnp
      need=wpadti(v+nnp)
      call getscm(need,a,maxcor,'intrn: mctrans3',0)
c
      call iosys('read real "overlap integrals" from rwf',
     $           nnp,z(s),0,' ')
      call iosys('read real "kinetic integrals" from rwf',
     $           nnp,z(h),0,' ')
      call iosys('read real "potential integrals" from rwf',
     $           nnp,z(v),0,' ')

c
      call vadd(z(h),z(h),z(v),nnp)
c
      t1=iadtwp(need)
      t2=t1+nbf**2
      t3=t2+nbf**2
      need=wpadti(t3+nbf**2)
      call getscm(need,a,maxcor,'intrn: mctrans4',0)
c
c     schmidt orthogonalize the vectors
c
      call schvec(z(c),z(s),z(t1),z(t2),z(t3),nbf,nnp)
c
      call chknrm(z(c),z(s),z(t1),z(t2),z(t3),nbf,nnp)
c
      call iosys('write real mo_old to mcscr',ntot,z(c),0,' ')
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
      if(nco.ne.0) then
         call corden(z(c),z(cdens),nbf,nco,nnp)
      endif
c
c     ----- form the frozen-core fock-matrix -----
c
      ecore=0.0d+00
      if (nco.gt.0) then
         cf=iadtwp(need)
         values=cf+nnp
         need=wpadti(values+min(60000,nnp**2))
c
         call getscm(need,a,maxcor,'intrn: mctrans5',0)
c
c        maxcor=maxcor-10
         ntriang=(iadtwp(need)-values)/nnp
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
         call iosys('write real "frozen core energy" to rwf',
     $        1,ecore,0,' ')
cps         call vmove(z(h),z(cf),nnp)
         call vadd(z(h),z(h),z(cf),nnp)
c
c     ----- store the core fock matrix and core density matrix ----
c
         nbf2=(nbf*(nbf+1))/2
         call iosys('write real mc_core_fock to mcscr',
     $              nbf2,z(h),0,' ')
         call iosys('write real "mcscf ao core density" to rwf',
     $               nbf2,z(cdens),0,' ')
c
      else
c
        eonel=0.d+00
        etwoel=0.d+00
        call iosys('write real "mcscf core 1e energy" to rwf',
     $              1,eonel,0,' ')
        call iosys('write real "mcscf core 2e energy" to rwf',
     $              1,etwoel,0,' ')
        call iosys('write real "frozen core energy" to rwf',
     $              1,ecore,0,' ')
      end if
c
c     ----- get nuclear repulsion energy from tape10 -----
c
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,enuc,0,' ')
c
c     ----- transform the one-electron integrals -----
c
c..bhl..closed shell scf
c
      nbf2=(nbf*(nbf+1))/2
c
      if(nao.ne.0) then
c
c     ----- step over core orbitals
c
         iof=nco*nbf
         call trn1e(z(h),z(c+iof),nbf,nao,nnp,numij,z(t1),z(t2),z(hmo))
         call iosys('does mc_tint1 exist on mcscr',0,0,0,ians)
         if(ians.eq.'no') then
            call iosys('create real mc_tint1 on mcscr',numij,0,0,' ')
         endif
c
         call iosys('write real mc_tint1 to mcscr',numij,z(hmo),0,' ')
         call iosys('write real "mo one-electron integrals" to mcscr',
     $               numij,z(hmo),0,' ')
c
      endif
c
c     ----- reallocate core for the two-electron transformation -----
c
      nnb=(nocc*(nocc+1))/2
      t1=c+nbf*nbf
      t2=t1+nbf**2
      valj=t2+nocc*nbf
      valk=valj+nnb*nnp
      values=valk+nnb*nbf*nbf
      ntriang=min(nnp,20000/nnp)
      nao2=nao*(nao+1)/2
      ntriang=max(ntriang,nao2**2/nnp+1)
      top=values+ntriang*nnp
c
c      val=values+ntriang*nnp
c      lenbin=min(640,nao2**2)
c      lab=val+lenbin
c      bin=lab+lenbin
c      maxipt=lenbin
c      asort=bin+lenbin
c
c      call getscm(0,a,maxr,'mctran55',0)
c      maxcor=iadtwp(maxr-512)
c      lnsort=maxcor-asort
c      lnsort=max(lnsort,10000)
c      need=wpadti(asort+lnsort)
c
      need=wpadti(top)
      if (ntriang.le.0) call lnkerr(' intrn: not enough core in bltr')
      call getscm(need,a,maxcor,'intrn: mctrans6',0)
c
      lst16=1
      jstrt=1
c
      if (nao.ne.0) then
         call reordc(z(c),z(t1),nco,nao,nbf)
      end if
c
      call trn2e(z(values),nnp,ntriang,z(c),nbf,nco,nocc,z(t1),z(t2),
     $     numij,z(valj),z(valk), nbufj,nbufk,z(values),z(t1),nao,
     $     lst16,kstrt)
c
      ihd(1)=nco
      ihd(2)=nao
      ihd(3)=nob
      ihd(4)=nbf
      ihd(5)=nblock
      ihd(6)=nbufj
      ihd(7)=nbufk
      ihd(8)=0
      ihd(16)=0
      ihd(17)=0
      ihd(18)=lsect1
cc
c  ..... put header record of j - k integral file
cc
      call iosys('does mc_jk_hdr exist on mcscr',0,0,0,ians)
      if(ians.eq.'no') then
         call iosys('create integer mc_jk_hdr on mcscr',20,0,0,' ')
      endif
c
      call iosys('write integer mc_jk_hdr on mcscr',20,ihd,0,' ')
c
c..bhl..closed shell scf
c
      if(nao.ne.0) then
c
       ntot2=numij*numij
       call iosys('does mc_tint2 exist on mcscr',0,0,0,ians)
       if(ians.eq.'no') then
          call iosys('create real mc_tint2 on mcscr',ntot2,0,0,' ')
       endif
c
       call iosys('write real mc_tint2 on mcscr',
     $             ntot2,z(values),0,' ')
       call iosys('write real "mo two-electron integrals"  on mcscr',
     $      ntot2,z(values),0,' ')
c
      endif
c
c
      return
      end
