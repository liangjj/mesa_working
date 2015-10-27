*deck @(#)mn814.f	5.1  11/6/94
      subroutine mn814(z,a,nder)
c***begin prologue     m814
c***date written       860819   (yymmdd)
c***revision date      891212   (yymmdd)
c
c   december 12, 1989   bhl at llnl
c      revised to process one degree of freedom at a time and
c      to read a chained derivative integral file
c
c   7 february 1988  bhl at brl
c      write derivative-core fock matrix to rwf for use
c      by link 1021 to construct dependent u-matrix terms
c
c   4 december 1986  pws at lanl
c      changing 'namint' and iosys open to character.
c
c***keywords           m814, link 814, derivative mcscf lagrangian
c***author             lengsfield, byron (brl)
c***source             @(#)mn814.f	5.1   11/6/94
c***purpose            to form the derivative mcscf lagrangian and
c                      derivative integrals for cpmcscf
c***description
c
c***references         b.h. lengsfield iii encylopedia galactica vol. i page 324
c
c
      implicit integer (a-z)
c
      character*4096 ops
      character*32 xform,ffnext
      character*3 ans
      logical logkey
      real*8 z(*)
      integer a(*)
c
      parameter (maxcor=20000)
c
      common /io/     inp,iout
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- find the name of the vector set to use in the transformation
c
c      write(iout,1234)
c 1234 format(/,'  m814: ')
c
      if(logkey(ops,'mcscf',.false.,'mcscf'))then
         xform='"mcscf vector"'
      else if(logkey(ops,'hf=quadratic',.false.,' '))then
         xform='"mcscf vector"'
      else
         xform='"scf vector"'
      endif
c
      pos=index(ops,'xform=')
      if (pos.gt.0) then
         if (ffnext(ops,pos+6,start,end).ne.'string') then
            call lnkerr('illegal transformation orbital set: '//
     #                   ops(pos:end))
         end if
         xform=ops(start:end)
      end if
c
c     ----- get the dimensions we need for core allocation -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $            1,nbf,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $            1,natoms,0,' ')
      call iosys('read integer mc_ncore from rwf',1,ncore,0,' ')
      call iosys('read integer mc_nactive from rwf',1,nactiv,0,' ')
c
      numder=natoms*3
      nnp=nbf*(nbf+1)/2
      nnpact=nactiv*(nactiv+1)/2
      nocc=nactiv+ncore
      nbfna=nbf*nactiv
      nbfnc=nbf*ncore
c
      if(nder.eq.1) then
         call iosys('rewind all on rwf',0,0,0,' ')
         call iosys('rewind all on rdints',0,0,0,' ')
c
         call iosys('does mo_der_lagrangian exist on rwf',
     #               0,0,0,ans)
         if(ans.eq.'no') then
            call iosys('create real mo_der_lagrangian on rwf',
     $                  nbf*nocc*numder,0,0,' ')
         end if
         call iosys('does half_der_lagrangian exist on rwf',
     #              0,0,0,ans)
         if(ans.eq.'no') then
            call iosys('create real half_der_lagrangian on rwf',
     $                  nbf*nocc*numder,0,0,' ')
         end if
c
         call iosys('does mo_der_overlap exist on dints',0,0,0,ans)
         if(ans.eq.'no') then
            call iosys('create real "ao derivative core lagrangian"'
     #                  //' on dints',numder*nnp,0,0,' ')
            call iosys('create real der_i(xbcd) on dints',
     #                 nbfna*nnpact*numder,0,0,' ')
            call iosys('create real der_g on dints',nnpact**2*numder,
     #                  0,0,' ')
            call iosys('create real der_h on dints',nnpact*numder,
     #                  0,0,' ')
            call iosys('create real der_f(x,b) on dints',nbfna*numder,
     #                  0,0,' ')
            call iosys('create real mo_der_overlap on dints',
     #                 numder*nbf*nbf,0,0,' ')
            call iosys('create real "ao derivative overlap integrals"'//
     #                 ' on dints',
     #                 numder*nnp,0,0,' ')
            if(ncore.gt.0) then
               call iosys('create real der_f(x,j) on dints',
     #                    nbfnc*numder,0,0,' ')
               call iosys('create real der_fab(x,k) on dints',
     #                     nbfnc*numder,0,0,' ')
            end if
         end if
c
         call iosys('rewind all on dints',0,0,0,' ')
c
      end if
c
c     ----- get  the transformation matrix (vector) and check for
c           orthonormality
c
      lag=1
      c=lag+nbf*nocc
      tmptop=c+nbf**2
      s=tmptop
      t1=s+nnp
      t2=t1+nbf**2
      t3=t2+nbf**2
      top=t3+nbf**2
c
      call getscm(top,z,junk,'m814 checking orthonormality of c',0)
c
      call iosys('read real '//xform//' from rwf',nbf**2,z(c),0,' ')
c
c      write(iout,*)'  orbitals read from ',xform
c      call matout(z(c),nbf,nbf,nbf,nbf,iout)
c
      call iosys('read real "overlap integrals" from rwf',
     $            nnp,z(s),0,' ')
c
      call chknrm(z(c),z(s),z(t1),z(t2),z(t3),nbf,nnp)
c
c
      ntriang=min(nnp,51200/nnp)
c
      ao=tmptop
c
c     ----- form the core fock matrix -----
c
      if (ncore.gt.0) then
         corden=ao+nnp
         values=corden+nnp
         t1=values+nnp*ntriang
         t2=t1+nbf**2
         top=t2+nbf**2
c
         call getscm(top,z,junk,'core-fock matrix construction',0)
c
         call cdens(z(corden),nnp,z(c),nbf,ncore,z(t1))
c
         call iosys('read real "ao derivative one-electron integrals"'//
     #              ' from rdints',nnp,z(ao),(nder-1)*nnp,' ')
c
c.bhl      call trtosq(z(t2),z(ao),nbf,nnp)
c.bhl      write(iout,*)' ndf=1 ao der. ints '
c.bhl      call matout(z(t2),nbf,nbf,nbf,nbf,iout)
c
         call dfock(z(values),z(corden),z(ao),nnp,nbf,z(t1),z(t2),
     #              ntriang,nder)
c
c..bhl
c
      call iosys('write real "ao derivative core lagrangian" to dints',
     #            nnp,z(ao),(nder-1)*nnp,' ')
c
c..bhl
c
      end if
c
c     ----- allocate memory for the one-electron transformations -----
c
      mo=ao
      fxj=ao+nnp
      fxb=fxj
      t1=fxj+max(ncore,nactiv)*nbf
      dab=t1+nbf**2
      top=dab+nactiv**2
c
      call getscm(top,z,junk,'m814, one-electron transformations',0)
c
      if (ncore.le.0) then
         call iosys('read real "ao derivative one-electron integrals"'//
     #' from rdints',nnp,z(ao),(nder-1)*nnp,' ')
      end if
c
c     ----- zero out the lagrangian, and start its formation as we go
c
      call rzero(z(lag),nbf*nocc)
c
      call trnd1e(z(c),z(ao),z(fxj),z(fxb),z(mo),z(t1),nbf,nnp,ncore,
     #            nactiv,nder,nnpact,nbf*ncore,nbf*nactiv,z(lag),
     #            z(dab),nocc)
c
c     ----- and the transformation of the overlap integrals -----
c
      mo=ao+nnp
      t1=mo+nbf**2
      t2=t1+nbf**2
      top=t2+nbf**2
c
      call getscm(top,z,junk,'derivative overlap transformations',0)
c
      call trnds(z(c),z(ao),z(mo),z(t1),z(t2),nbf,nnp,nder,nbf**2)
c
c     ----- and for active-orbital density section -----
c
      if (ncore.gt.0) then
         fxk=ao+nnp
         t1=fxk+nbf*ncore
         t2=t1+nbf**2
         d=t2+nbf**2
         values=d+nbf**2
         top=values+nnp*ntriang
c
         call getscm(top,z,junk,'active-orbital density',0)
c
c        ----- get density and transform to ao basis -----
c
         call fmdab(z(d),z(c),z(t1),z(t2),nactiv,nnpact,nbf,nnp,ncore)
c
c        ----- form fock matrices -----
c
         call rzero(z(ao),nnp)
         call dfock(z(values),z(d),z(ao),nnp,nbf,z(t1),z(t2),
     #              ntriang,nder)
c
c        ----- transform -----
c
         call fab(z(c),z(fxk),z(t1),z(ao),nbf,ncore,nbf*ncore,nnp,
     #            z(lag),nder,nocc)
c
      end if
c
c     ----- allocate space for two-electron section:
c           3 index transformation
c
c
      call getscm(0,z,canget,'maximum size?',0)
c
      t1=tmptop
      t2=t1+nbf**2
      third=t2+nbf**2
      mo=third+nbf*nactiv*nnpact
      lab=wpadti(mo)
      values=iadtwp(lab+nnp)
      asort=values+nnp*ntriang
      top=min(canget,asort+nnp**2)
      top=max(top,mo+nnpact**2)
c
      call getscm(top,z,junk,'derivative 3-index transformation',0)
c
c     ----- three-index transformations -----
c
      lnsort=top-asort+1
c
      call trnd2e(z(c),z(values),z(t1),z(t2),z(mo),a(lab),
     #            nbf,nnp,nactiv,ncore,nnpact,ntriang,nbf*nactiv,
     #            z(asort),lnsort,nder,z(third))
c
c     ----- add two-electron portion to the lagrangian -----
c
      sqdm=tmptop
      twopdm=sqdm
      t1=sqdm+nactiv**2*nnpact
      xbcd=t1
      top=max(t1+nnpact**2,xbcd+nbf*nactiv*nnpact)
c
      call getscm(top,z,junk,'lagrangian construction',0)
c
      call lagran(z(lag),z(xbcd),z(sqdm),z(t1),z(twopdm),nactiv,
     #            nbf,nnpact,nactiv**2,nbf*nactiv*nnpact,
     #            nnpact*(nnpact+1)/2,nder,nocc,ncore)
c
c     ----- transform the lagrangians to the mo basis -----
c
      molag=tmptop
      top=molag+nbf*nocc
c
      call getscm(top,z,junk,'tranformation of the lagrangians',0)
c
      call trnlag(z(c),z(lag),z(molag),nbf,nocc,nder)
c
c     ----- close integrals unit -----
c
c..bhl      call iosys('close ints',0,0,0,' ')
c
c     ----- time to chain to next link -----
c
c..bhl      call chainx(0)
c
c
      return
      end
