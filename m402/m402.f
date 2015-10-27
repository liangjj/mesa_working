*deck @(#)pm402.f	5.1  11/6/94
      program m402
c
c***begin prologue     m402
c***date written       900410   (yymmdd)
c***revision date               (yymmdd)
c
c***keywords           m402, symmetry, hartree-fock, scf,
c***author             martin, richard (lanl)
c***source             @(#)pm402.f	5.1   11/6/94
c***purpose            projects eigenvectors into pure irreducible
c                      representations.
c***description
c
c***references
c
c***routines called
c***end prologue       m402
c
      implicit integer (a-z)
c
      parameter (maxrep=14)
      character*4096 ops
      character*8 lirrep(maxrep),prtflg
      character*32 xform
c
      logical debug,prnt,logkey,drop
c
      integer a
      integer numso(maxrep), lambda(maxrep), numsn(maxrep)
c
      logical kohn
c
      real*8 z
c
      common /io/ inp,iout
c
      data maxcor/1/, prnt/.false./, debug/.false./
      save maxcor,prnt,debug
      pointer(p,z(1)), (p,a(1))
c
 1000 format(1x,'raw salcs')
 1010 format(1x,'salc overlap')
 1020 format(1x,'orthogonalized salcs')
 1025 format(1x,'projected vectors')
 1030 format(1x,'m402:project symmetry vectors')
 1040 format(5x,'orbital  index:',20i4,(/20x,20i4))
 1050 format(5x,'symmetry index:',20i4,(/20x,20i4))
 1060 format(5x,'symmetry label:',20(1x,a3),(/20x,20(1x,a3)))
 1070 format(5x,'calculation performed in reduced basis')
      call drum
      write (iout,1030)
c
c     ----- recover the options string -----
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     has printing been turned off externally?
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
      if (logkey(ops,'print=m402=all',.false.,' ')) then
          prnt=.true.
          debug=.true.
      endif
c
c     ----- is this a kohn calculation? -----
      kohn=logkey(ops,'kohn',.false.,' ')
c
c     ----- get the dimensions, symmetry information, etc -----
c
c----------------------------------------------------------------------c
c         drop is a flag which made m330 drop function components      c
c         oldnbf is the original basis set with no functions deleted   c
c                nbf is the size of the deleted set                    c
c----------------------------------------------------------------------c
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbf,0,' ')
      drop=logkey(ops,'drop',.false.,' ')
      if (drop) then
         call iosys('read integer "old nbf" from rwf',1,oldnbf,0,' ')
         write(iout,1070)
      else
         oldnbf=nbf
      endif
      call iosys('read integer "number of irreducible representations"'
     $          //' from rwf',1,nirrep,0,' ')
      call iosys('read integer "number of symmetry orbitals" from rwf',
     $           nirrep,numso,0,' ')
      call icopy(numso,numsn,nirrep)
      call iosys('read integer "degeneracies of irreducible'//
     $           ' representations" from rwf',nirrep,lambda,0,' ')
      call iosys('read character "labels of irreducible '//
     $           'representations" from rwf',nirrep*len(lirrep(1)),
     $            0,0,lirrep)
c
c
c     allocate core.
c
      nnp=(nbf+1)*nbf/2
      s=1
      smhalf=s+nnp
      u=smhalf+nnp
      usym=u+nbf*nbf
      eigval=usym+nbf*nbf
      salc=eigval+nbf
      salct=salc
      if (drop) then
          salct=salc+nbf*nbf
      endif
      bftoso=salct+oldnbf*oldnbf
      t1=bftoso+nbf*nbf
      t2=t1+nbf*nbf
      t3=t2+nbf*nbf
      t4=t3+nbf*nbf
      t5=t4+nbf
      norm=t5+nbf*nbf
      triang=norm+nbf*nirrep
      vecrep=wpadti(triang+nnp)
      bfns=vecrep+nbf
      indx=bfns
      if (drop) then
          indx=bfns+nbf
      endif
      top=indx+oldnbf
c
c      call getscm(top,z,maxcor,'m402',0)
      call getmem(top,p,ngot,'m402',0)
c
c     read in the overlap matrix.
      call iosys('read real "overlap integrals" from rwf',
     $            nnp,z(s),0,' ')
      call iosys('read real "salc transformation matrix" from rwf',
     $            -1,z(salct),0,' ')
      if (drop) then
          call nsalc(z(salct),z(salc),a(indx),nirrep,numso,numsn,
     1               oldnbf,nbf,prnt)
      endif
      if (kohn) then
         call aosym(z(salc),numsn,a(bfns),nirrep,nbf)
      endif
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
      call iosys('write real "orthogonal salc transformation matrix"'
     $          //' to rwf',nbf*nbf,z(bftoso),0,' ')
c
c     project the vectors, good stuff comes back in usym.
c
c     read in the transformation vector.
      xform='"guess vector"'
      call iosys('read real '//xform//' from rwf',-1,z(u),0,' ')
c
      call trtosq(z(t5),z(s),nbf,nnp)
      nmo=nbf
      call symprj(nirrep,nbf,nbf,numsn,lambda,lirrep,z(u),nmo,z(bftoso),
     $            z(usym),a(vecrep),z(t5),z(t1),z(t2),z(t3),z(t4),
     $            z(norm),ops)
      if (drop) then
c----------------------------------------------------------------------c
c         overwrite original salc matrix with new salc matrix          c
c----------------------------------------------------------------------c
         call iosys('write integer "number of symmetry orbitals" to '
     1            //'rwf',nirrep,numsn,0,' ')
         call iosys('write real "salc transformation matrix" to rwf',
     1            nbf*nbf,z(salc),0,' ')
      end if
c
c     save the projected vectors on the rwf, and store their
c     associated symmetry indices.
      call iosys('write real '//xform//' to rwf',nbf*nbf,z(usym),
     $           0,' ')
      call iosys('write integer "guess vector symmetries" '
     $           //'to rwf',nbf,a(vecrep),0,' ')
c
      if(prnt) then
         write (iout,1025)
         call matout(z(usym),nbf,nbf,nbf,nbf,iout)
      end if
      prnt=logkey(ops,'print=m402=symmetry-information',.false.,' ')
      if (prnt) then
         write(iout,1040) (i,i=1,nbf)
         write(iout,1050) (a(vecrep+i-1),i=1,nbf)
         write(iout,1060) (lirrep(a(vecrep+i-1)),i=1,nbf)
      end if
c
c
      call iosys('write character ran402 to rwf',0,0,0,'yes')
c
c     ----- and exit gracefully -----
c
      call getmem(-ngot,p,idum,'m402',idum)
      call chainx(0)
c
c
      stop
      end
