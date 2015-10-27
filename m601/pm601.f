*deck  @(#)pm601.f	5.1 11/6/94
      subroutine pm601(z,a)
c
c***begin prologue     m601
c***date written       900410   (yymmdd)
c***revision date               (yymmdd)
c
c***keywords           m601, symmetry, hartree-fock, scf,
c***author             martin, richard (lanl)
c***source              @(#)pm601.f	5.1 11/6/94
c***purpose            projects eigenvectors into pure irreducible
c                      representations.
c***description
c
c***references
c
c***routines called
c***end prologue       m601
c
      implicit integer (a-z)
c
      parameter (maxrep=14)
      character*4096 ops
      character*128 namchk
      character*8 lirrep(maxrep),prtflg
      character*32 chrkey,xform
c
      logical debug, prnt, logkey
c
      integer a(*)
      integer numso(maxrep), lambda(maxrep)
c
      real*8 z(*)
c
      common /io/ inp,iout
c
      data maxcor/1/, prnt/.false./
      save maxcor,prnt
      parameter (debug=.false.)
c
 1000 format(1x,'raw salcs')
 1010 format(1x,'salc overlap')
 1020 format(1x,'orthogonalized salcs')
 1025 format(1x,'projected vectors')
 1030 format(1x,'m601:')
 1040 format(5x,'orbital  index:',20i4,(/20x,20i4))
 1050 format(5x,'symmetry index:',20i4,(/20x,20i4))
 1060 format(5x,'symmetry label:',20(1x,a3),(/20x,20(1x,a3)))
c
c     ----- recover the options string -----
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     has printing been turned off externally?
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
c
c     ----- get the dimensions, symmetry information, etc -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbf,0,' ')
      nnp=(nbf+1)*nbf/2
      call iosys('read integer "number of irreducible representations"'
     $          //' from rwf',1,nirrep,0,' ')
      call iosys('read integer "number of symmetry orbitals" from rwf',
     $           nirrep,numso,0,' ')
      call iosys('read integer "degeneracies of irreducible'//
     $           ' representations" from rwf',nirrep,lambda,0,' ')
      call iosys('read character "labels of irreducible '//
     $           'representations" from rwf',nirrep*len(lirrep(1)),
     $            0,0,lirrep)
c
c
c     allocate core.
c
      s=1
      smhalf=s+nnp
      u=smhalf+nnp
      usym=u+nbf*nbf
      eigval=usym+nbf*nbf
      salc=eigval+nbf
      bftoso=salc+nbf*nbf
      t1=bftoso+nbf*nbf
      t2=t1+nbf*nbf
      t3=t2+nbf*nbf
      t4=t3+nbf*nbf
      norm=t4+nbf
      triang=norm+nbf*nirrep
      vecrep=wpadti(triang+nnp)
      end=vecrep+nbf
c
      call getscm(end,z,maxcor,'m601',0)
c
c     read in the overlap matrix.
      call iosys('read real "overlap integrals" from rwf',
     $            nnp,z(s),0,' ')
      call iosys('read real "salc transformation matrix" from rwf',
     $            -1,z(salc),0,' ')
      if(debug) then
         write(iout,1000)
         call matout(z(salc),nbf,nbf,nbf,nbf,iout)
      end if
c
c     symmetrically orthogonalize the basis function to salc matrix.
      call trtosq(z(t1),z(s),nbf,nnp)
      call ebc(z(t2),z(t1),z(salc),nbf,nbf,nbf)
      call ebtc(z(t1),z(salc),z(t2),nbf,nbf,nbf)
      if(debug) then
         write(iout,1010)
         call matout(z(t1),nbf,nbf,nbf,nbf,iout)
      end if
      call sqtotr(z(usym),z(t1),nbf,nnp)
c     usym now holds the salc x salc overlap matrix.
      iprint=0
      call sinv(z(usym),z(smhalf),z(u),z(eigval),z(t1),z(t2),
     $          nbf,nnp,z(triang),iprint)
      call trtosq(z(u),z(smhalf),nbf,nnp)
      call ebc(z(bftoso),z(salc),z(u),nbf,nbf,nbf)
c
      if(debug) then
         write(iout,1020)
         call matout(z(bftoso),nbf,nbf,nbf,nbf,iout)
      end if
c
c     project the vectors, good stuff comes back in usym.
c
c     read in the transformation vector, perhaps overriding the directive
c     on the rwf file.
      xform=chrkey(ops,'transformation=vector',' ',' ')
      if(xform.eq.' ') then
         call iosys('read character "transformation vector" from rwf',
     $              -1,0,0,xform)
      else
         call iosys('write character "transformation vector" to rwf',
     $              0,0,0,xform)
      end if
      if(logkey(ops,'transformation=chk',.false.,' ')) then
         call iosys('read character "checkpoint filename" from rwf',
     $               0,0,0,namchk)
         call iosys('open chk as unknown',0,0,0,namchk)
         call iosys('read real '//xform//' from chk',-1,z(u),0,' ')
         call iosys('close chk',0,0,0,' ')
      else
         call iosys('read real '//xform//' from rwf',-1,z(u),0,' ')
      end if
c
      call trtosq(z(salc),z(s),nbf,nnp)
      nmo=nbf
      call symprj(nirrep,nbf,nbf,numso,lambda,lirrep,z(u),nmo,z(bftoso),
     $            z(usym),a(vecrep),z(salc),z(t1),z(t2),z(t3),z(t4),
     $            z(norm))
c
c     save the projected vectors on the rwf, and store their
c     associated symmetry indices.
      call iosys('write real '//xform//' to rwf',nbf*nbf,z(usym),
     $           0,' ')
      call iosys('write integer "transformation vector symmetries" '
     $           //'to rwf',nbf,a(vecrep),0,' ')
c
      write(iout,1030)
      if(prnt) then
         write (iout,1025)
         call matout(z(usym),nbf,nbf,nbf,nbf,iout)
      end if
      prnt=.true.
      if (prnt) then
         write(iout,1040) (i,i=1,nbf)
         write(iout,1050) (a(vecrep+i-1),i=1,nbf)
         write(iout,1060) (lirrep(a(vecrep+i-1)),i=1,nbf)
      end if
c     ----- and exit gracefully -----
c
      call chainx(0)
c
c
      stop
      end

