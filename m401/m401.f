*deck @(#)pm401.f	5.2  4/17/95
      program m401
c***begin prologue     m401
c***date written       850810   (yymmdd)
c***revision date      931994   (yymmdd)
c
c  19 january 1994     rlm at lanl
c     adding symmetry projection routines from m402
c
c  01 august  1990     bhl at llnl
c     option to expand natural orbitals with set specification
c
c  22 july    1989     bhl at llnl
c     option guess=merge added to pickup two basis sets
c     and merge them. a sav unit was introduced and
c     an array ind(3,200) specifies the orbital pickup
c
c  19 january 1988     bhl at brl
c     option changed so guess=chk ..etc. is accepted
c     as well as scf=guess=chk
c
c  09 january 1988     bhl at brl
c     print option changed from print=scf=guess to print=guess
c
c  15 may 1987         pws at lanl
c     fixing up the core allocation for 32/64 bit machines
c
c***keywords           m401, link 401, guess, $alter, huckel, core,
c                      rdchk, rdrwf, rdinp, $vectors
c***author             martin, richard    (lanl)
c***source             @(#)pm401.f	5.2   4/17/95
c***purpose            prepares an initial guess for the scf link.
c***description
c     the type of guess desired is specified in the $route section of the
c     input deck as guess=(option).
c     m401 currently recognizes:
c       guess=core    in this case the guess comes from diagonalizing the
c                     core hamiltonian (h=t+v), where t is the kinetic
c                     energy and v is the nuclear attraction operator.
c                     this type of guess has the advantage that it will
c                     work for any basis set, but is generally a fairly
c                     poor starting point.
c       guess=huckel  an extended huckel calculation is performed to obtain
c                     starting vectors.  the huckel calculation is done
c                     in a mimimum(sto3g) basis set, and the occupied
c                     orbitals are projected onto the current basis using
c                     the corresponding orbital transformation.  this type
c                     of guess is also available for any basis set, and is
c                     usually fairly decent.
c       guess=rdchk   the initial vectors are read from the checkpoint file
c                     and projected onto the current basis.  this is
c                     usually an excellent guess if you have performed a
c                     previous calculation in a  different basis, on a
c                     different electronic state, or at a nearby geometry.
c       guess=rdrwf   the initial vectors are read from the read/write
c                     file, and projected onto the current basis.
c       guess=rdinp   the initial guess is read from the input file.
c                     the input is signaled through a $vectors section,
c                     followed by each of the occupied orbitals, read in
c                     a list directed format.
c       guess=merge   two vector sets are read and merged. the vector
c                     are specified by defining set1 and set2 in the
c                     option string. an ind array is used to determine
c                     which vectors will be used. the first element of ind
c                     species the vector set, the next element the starting vector and
c                     the third element is the last vector, so ind is a
c                     series of triples.
c                     ex:  guess=(merge	,ind=( 1,1,3, 2,1,3, 1,7,20 ))
c                          set1=chk  set2=sav
c                     the first  3 vectors are orbitals 1-3 from set1
c                     the second 3 vectors are orbitals 1-3 from set2
c                     the remaining vectors are taken from set1
c
c     additional options:
c       alter         after the initial guess is generated, it may be
c                     modified by a $alter section in the input.
c                     this directive is followed by a list of orbital pairs.
c                     for example, the input section
c                       $alter
c                        6,7
c                       $end
c                     will interchange orbital 6 with 7 before storing
c                     them for the scf link.
c       print=guess   print the initial orbitals.
c       timing        print timing statistics for this link.
c
c***references
c
c***routines called
c     m401
c       drum(mdutil)
c       iosys(io)
c       traktm(mdutil)
c       getscm(mdutil)
c       sinv(util)
c       prjges(local)
c         sizmb(local)
c         rdbas(local)
c         rdmo(local)
c         minbas(local)
c           replic(util)
c           sto3g(local)
c             coefs(local)
c             scalef(local)
c             defsc(local)
c           togenc(util)
c           fmstrt(util)
c           ehuckl(local)
c             eneg(local)
c           normal(util)
c         ovrlap(local)
c           stwod(util)
c           fmonel(util)
c           trans1(util)
c           putrec(local)
c         hukmat(local)
c           getmo(local)
c           alter(local)
c       rwfges(local)
c         rdmo(local)
c         alter(local)
c       corges(local)
c         getmo(local)
c         alter(local)
c       schmdt(math)
c       wvec(util)
c       chainx(mdutil)
c
c***end prologue       m401
      implicit integer(a-z)
      parameter (maxnbf=2000)
      parameter (maxrep=14)
      integer a, isym, ngot(2)
      real*8 z, sym
      real*8 maxerr, fpkey
      character*4096 ops
      character*8 gestyp,prtflg,filenm
      character*16 xform
      character*16 bflabl(maxnbf)
      character*16 symlabl(maxnbf)
      character*3 answer
      character*8 lirrep(maxrep)
      logical altges,prnt
      logical logkey
      logical usesym
c
      common/io/inp,iout
      pointer(p,z(1)), (p,a(1))
      pointer(psym,sym(1)), (psym,isym(1))
c
      data maxerr/1.0d-06/, altges/.false./, prnt/.true./
      save maxerr,altges,prnt
c
 1000 format(1x,'m401:initial guess')
 1010 format(5x,'projected checkpoint.')
 1020 format(5x,'read from the input file.')
 1030 format(5x,'projected huckel.')
 1031 format(5x,'expanding vectors from chk')
 1032 format(5x,'merging vectors ')
 1040 format(5x,'read from the rwf.')
 1045 format(5x,'read from the chk.')
 1050 format(5x,'core hamiltonian.')
 1100 format(5x,'initial guess vector.')
c
c
c     retrieve the options string.
      call iosys('read character options from rwf',-1,0,0,ops)
      gestyp='core'
      if(logkey(ops,'scf=guess=rdchk',.false.,' ')) gestyp='rdchk'
      if(logkey(ops,'scf=guess=rwf',.false.,' ')) gestyp='rwf'
      if(logkey(ops,'scf=guess=chk',.false.,' ')) gestyp='chk'
      if(logkey(ops,'scf=guess=rdinp',.false.,' ')) gestyp='rdinp'
      if(logkey(ops,'scf=guess=core',.false.,' ')) gestyp='core'
      if(logkey(ops,'scf=guess=huckel',.false.,' ')) gestyp='huckel'
      altges=logkey(ops,'scf=guess=alter',.false.,' ')
c
      if(logkey(ops,'guess=rdchk',.false.,' ')) gestyp='rdchk'
      if(logkey(ops,'guess=rwf',.false.,' ')) gestyp='rwf'
      if(logkey(ops,'guess=chk',.false.,' ')) gestyp='chk'
      if(logkey(ops,'guess=rdinp',.false.,' ')) gestyp='rdinp'
      if(logkey(ops,'guess=core',.false.,' ')) gestyp='core'
      if(logkey(ops,'guess=huckel',.false.,' ')) gestyp='huckel'
      if(logkey(ops,'guess=expand',.false.,' ')) gestyp='expand'
      if(logkey(ops,'guess=merge',.false.,' ')) gestyp='merge'
      if(logkey(ops,'guess=fixchk',.false.,' ')) gestyp='fixchk'
      altges=logkey(ops,'guess=alter',altges,' ')
      maxerr=fpkey(ops,'maxerr',1.d-06,' ')
      usesym=.not.logkey(ops,'sym=off',.false.,' ')
      if(logkey(ops,'drop',.false.,' ')) then
         usesym=.false.
      endif
c
c     modify the local options if we're in the course of an optimization.
c     specifically, get the guess from the rwf and don't alter it.
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
c
      if(prtflg.eq.'minimum') then
         gestyp='rwf'
         prnt=.false.
         altges=.false.
      endif
c
c     signal our presence.
      if(prnt) then
         write(iout,1000)
         if(gestyp.eq.'rdchk') write(iout,1010)
         if(gestyp.eq.'rdinp') write(iout,1020)
         if(gestyp.eq.'huckel') write(iout,1030)
         if(gestyp.eq.'expand') write(iout,1031)
         if(gestyp.eq.'merge') write(iout,1032)
         if(gestyp.eq.'rwf') write(iout,1040)
         if(gestyp.eq.'chk') write(iout,1045)
         if(gestyp.eq.'fixchk') write(iout,1045)         
         if(gestyp.eq.'core') write(iout,1050)
      endif
c
c     retrieve the basis set information.
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbasis,0,' ')
      nnp=nbasis*(nbasis+1)/2
      nbsq=nbasis*nbasis
c
c     allocate core needed by all options.
      c=1
      eigval=c+nbsq
      s=eigval+nbasis
      smhalf=s+nnp
c     scratch.
      u=smhalf+nnp
      iu=wpadti(u)
      t1=u+nbsq
      t2=t1+nbsq
      t3=t2+nbsq
      hmat=t3+nbsq
      triang=hmat+nnp
      top=wpadti(triang+nnp)
c     core guess
      t=u
      v=t+nnp
      top=max(top,wpadti(v+nnp))
c     checkpoint guess.
c
c      call getscm(0,z,maxcor,'m401',0)
      call getmem(top,p,ngot(1),'m401',0)
c
c     form s**-1/2.
      call iosys('read real "overlap integrals" from rwf',
     $           nnp,z(s),0,' ')
      iprint=0
      call sinv(z(s),z(smhalf),z(u),z(eigval),z(t1),z(t2),
     $          nbasis,nnp,z(triang),iprint)
c
c     save a copy of s**-1/2
      call iosys('write real "overlap to minus half" to rwf',
     $            nnp,z(smhalf),0,' ') 
c
c     obtain the guess.
      if(gestyp.eq.'rdchk'.or.gestyp.eq.'huckel'
     $                    .or.gestyp.eq.'rdinp') then
         call prjges(gestyp,nbasis,nnp,z(c),z(eigval),
     $               z(smhalf),altges,ops)
      else if(gestyp.eq.'rwf'.or.gestyp.eq.'chk'.or.
     $        gestyp.eq.'fixchk') then
         filenm=gestyp
         if ( gestyp.eq.'fixchk') then
              filenm='chk'
         endif
         call rdges(nbasis,z(c),z(eigval),z(triang),filenm,altges,ops)
      else if(gestyp.eq.'core') then
         call corges(z(smhalf),z(t),z(v),z(hmat),z(c),z(eigval),
     $               z(triang),z(u),z(t1),nbasis,nnp,altges)
      else if(gestyp.eq.'expand') then
         if(logkey(ops,'guess=expand=no',.false.,' ')) then
            set=2
            rrtt=intkey(ops,'guess=expand=no',1,' ')
            write(iout,*)'  expanding natural orbital set ',rrtt
         else if(logkey(ops,'guess=expand=mcscf',.false.,' ')) then
            set=1
         else
            set=0
         end if
         call expand(z(c),z(t1),nbasis,nsmall,set,z(eigval),rrtt)
         if(logkey(ops,'print=eguess',.false.,' ')) then
            write(iout,1100)
            call iosys('read character "basis function labels"'
     $                 //' from rwf',len(bflabl(1))*nbasis,0,0,bflabl)
            call wvec(z(c),z(eigval),nbasis,nbasis,bflabl,' ')
         endif
      else if(gestyp.eq.'merge') then
         call merge(z(c),z(t1),z(t2),z(eigval),z(t3),z(u),nbasis,ops)
      else
         call lnkerr('unrecognized guess type in m401:'//gestyp)
      endif
c
c     normalize the guess.
      call schmdt(z(c),z(s),z(u),z(t1),z(t2),nbasis,nbasis,nnp,maxerr)
c
c
c     print the guess if requested.
      if(logkey(ops,'print=guess',.false.,' ')) then
         if(logkey(ops,'print=guess=all',.false.,' ')) then
            write(iout,1100)
            call iosys('read character "basis function labels" from'
     $               //' rwf',len(bflabl(1))*nbasis,0,0,bflabl)
            call wvec(z(c),z(eigval),nbasis,nbasis,bflabl,' ')
         else
            call iosys('read integer "number of alpha electrons" from'//
     $                 ' rwf',1,nummo,0,' ')
            nummo=min(nbasis,nummo+5)
            write(iout,1100)
            call iosys('read character "basis function labels" from'//
     $                 ' rwf',len(bflabl(1))*nbasis,0,0,bflabl)
            call wvec(z(c),z(eigval),nbasis,nummo,bflabl,' ')
         endif
      endif
c
c
      if(gestyp.eq.'expand') then
         if(.not.logkey(ops,'guess=norot',.false.,' ')) then
             nsmll=intkey(ops,'guess=fix',0,' ')
             if(nsmll.ne.0) then
                msmll=nbasis-nsmll
                write(iout,*)' rotating all orbitals greater than ',
     $                         nsmll
             else
                msmll=nsmall
             endif
             call expges(z(t),z(v),z(hmat),z(c),z(eigval),
     $                   z(triang),z(u),z(t1),z(t2),nbasis,nnp,msmll)
         else
             write(iout,*)' the expanded orbitals will not be rotated'
         endif
         if(logkey(ops,'symmetry',.false.,' ')) then
            call symorb(z(c),z(triang),z(t1),z(t2),z(t3),z(u),nbasis,
     $                  nnp,nsmall,ops)
         endif
         call schmdt(z(c),z(s),z(u),z(t1),z(t2),nbasis,nbasis,nnp,
     $               maxerr)
         if(logkey(ops,'guess=alter',.false.,' ')) then
            call alter(z(c),z(eigval),z(triang),nbasis)
         end if
         if(logkey(ops,'print=guess',.false.,' ')) then
            write(iout,1100)
            call wvec(z(c),z(eigval),nbasis,nbasis,bflabl,' ')
         endif
      elseif ( gestyp.eq.'fixchk') then
         if(.not.logkey(ops,'guess=norot',.false.,' ')) then
            nsmll=intkey(ops,'guess=fix',0,' ')
            if(nsmll.ne.0) then
               msmll=nbasis-nsmll
               write(iout,*)' rotating all orbitals greater than ',nsmll
            else
               msmll=nsmall
            endif
            call expges(z(t),z(v),z(hmat),z(c),z(eigval),
     $                  z(triang),z(u),z(t1),z(t2),nbasis,nnp,msmll)
         else
            write(iout,*)' the expanded orbitals will not be rotated'
         endif
         if(logkey(ops,'print=guess',.false.,' ')) then
            write(iout,1100)
            call wvec(z(c),z(eigval),nbasis,nbasis,bflabl,' ')
         endif
      endif
c
c     --- project pure symmetry vectors if using symmetry ---
      if(usesym) then
         salc=1
         bftoso=salc+nbasis*nbasis
         lambda=wpadti(bftoso+nbasis*nbasis)
         vecrep=lambda+maxrep
         numso=vecrep+nbasis
         top=numso+maxrep
         call getmem(top,psym,ngot(2),'sym',0)
c
c        --- get the dimensions, symmetry information, etc ---
         call iosys('read integer "number of irreducible'
     $              //' representations" from rwf',1,nirrep,0,' ')
         call iosys('read integer "number of symmetry orbitals"'
     $              //' from rwf',nirrep,isym(numso),0,' ')
         call iosys('read integer "degeneracies of irreducible'//
     $              ' representations" from rwf',
     $              nirrep,isym(lambda),0,' ')
         call iosys('read character "labels of irreducible'//
     $              ' representations" from rwf',nirrep*len(lirrep(1)),
     $               0,0,lirrep)
         call iosys('read real "salc transformation matrix" from rwf',
     $               -1,sym(salc),0,' ')
         call symtriz(ops,.false.,nbasis,nnp,z(s),z(smhalf),z(c),
     $                sym(bftoso),isym(numso),sym(salc),isym(vecrep),
     $                nirrep,isym(lambda),lirrep)
c
c        --- save the projected vector symmetry indices and the 
c            orthogonalized salc transformation matrix ---
         call iosys('write real "orthogonal salc transformation '//
     $              'matrix" to rwf',nbasis*nbasis,sym(bftoso),0,' ')
         call iosys('write integer "guess vector symmetries"'
     $              //' to rwf',nbasis,isym(vecrep),0,' ')
         if (logkey(ops,'print=guess=salc',.false.,' ')) then
            write(iout,*)'Symmetrized guess vector'
            do 10 i=1,nbasis
               symlabl(i)=lirrep(isym(vecrep+(i-1)))
 10         continue 
            call matprt(z(c),nbasis,nbasis,nbasis,nbasis,1,2,bflabl,
     $           symlabl,0,z(eigval),.true.)
         endif
         call getmem(-ngot(2),psym,idum,'sym',idum)
      endif
c
c     --- write the new guess to the rwf ---
      call iosys('write real "guess vector" to rwf',nbsq,z(c),0,' ')
      xform='"guess vector"'
      call iosys('write character "transformation vector" to rwf',
     $           0,0,0,xform)
      call iosys('write real "guess orbital energies" to rwf',nbasis,
     $           z(eigval),0,' ')
      call iosys('write real "orbital energies" to rwf',nbasis,
     $           z(eigval),0,' ')
      call iosys('write real "scf vector" to rwf',nbsq,z(c),0,' ')
c
c
      call getmem(-ngot(1),p,idum,'m401',idum)
      call chainx(0)
c
c
      stop
      end
