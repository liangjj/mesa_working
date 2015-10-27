*deck @(#)rdges.f	5.1  11/6/94
      subroutine rdges(nbasis,c,eigval,scrtch,file,altges,ops)
c***begin prologue     rdges
c***date written       850601  yymmdd
c***revision date      880127  yymmdd
c
c  27 january  1988    bhl at brl
c       if file=rwf code modified for mcscf optimization
c
c  19 january  1988    bhl at brl
c       modifying for guess=chk=no type options
c
c   2 november 1987    pws at lanl
c       modifying for scf=guess=chk=mcscf type options
c
c***keywords           initial guess, molecular orbitals, rwf, chk, i/o
c***author             martin, richard (lanl)
c***source             @(#)rdges.f	5.1   11/6/94
c***purpose            reads molecular orbitals/eigenvalues from 'file'.
c***description
c     call rdges(nbasis,c,eigval,scrtch,file,altges)
c       nbasis  the basis set dimension.
c       c       the molecular orbital array(nbasis,nbasis).
c       eigval  the eigenvalue vector(nbasis).
c       scrtch  scratch vector(nbasis).
c       file    character string; 'rwf' or 'chk'.
c       altges  logical flag (.true., alter the guess).
c***references         (none)
c***routines called    rdmo(m401), alter(m401)
c***end prologue       rdges
c
      implicit integer(a-z)
c
      character*(*) ops
      logical logkey
      real*8 c(nbasis,nbasis),eigval(nbasis),scrtch(nbasis)
      character*(*) file
      character temp*4, namchk*128, itoc*4
      logical altges
c
c.. bhl.. read from rwf changed for mcscf optimization
c
      nnbf=nbasis*nbasis
      nbf=nbasis
c..bhl 7/10/90
      if(file.eq.'rwf') then
         if (logkey(ops,'scf=guess=rwf=mcscf',.false.,' ').or.
     $       logkey(ops,'guess=rwf=mcscf',.false.,' ').or.
     $       logkey(ops,'mcscf',.false.,' ')) then
            call iosys('read real "mcscf vector" from rwf',nnbf,c,0,' ')
            call iosys('read real "mcscf orbital energies" from rwf',
     $           nbf,eigval,0,' ')
         else
            call iosys('read real "scf vector" from rwf',nnbf,c,0,' ')
            call iosys('read real "orbital energies" from rwf',
     $           nbf,eigval,0,' ')
         end if
      else if(file.eq.'chk') then
c
c        prepare the checkpoint file for action.
         call iosys('read character "checkpoint filename" from rwf',
     $        0,0,0,namchk)
         call iosys('open chk as old',0,0,0,namchk)
c
         if (logkey(ops,'scf=guess=chk=mcscf',.false.,' ').or.
     $       logkey(ops,'guess=chk=mcscf',.false.,' ')) then
            call iosys('read real "mcscf vector" from chk',nnbf,c,0,' ')
            call iosys('read real "mcscf orbital energies" from chk',
     $           nbf,eigval,0,' ')
c..bhl
         else if (logkey(ops,'guess=chk=no',.false.,' ')) then
c..bhl
            root=intkey(ops,'guess=chk=no',1,' ')
            call iosys('read real "no vector '//itoc(root)//'"'
     $                 //' from chk',nnbf,c,0,' ')
            call iosys('read real "no occ '//itoc(root)//'" from chk',
     $                 nbf,eigval,0,' ')
         else if (logkey(ops,'guess=chk=unconverged-scf',
     $            .false.,' ')) then
            call iosys('read real "scf unconverged vector" from chk',
     $                 nnbf,c,0,' ')
            call iosys('read real "scf unconverged orbital energies"'
     $                 //' from chk',nbf,eigval,0,' ')
         else if (logkey(ops,'guess=chk=unconverged-mcscf',
     $            .false.,' ')) then
            call iosys('read real "mcscf unconverged vector" from chk',
     $                 nnbf,c,0,' ')
            call iosys('read real "mcscf unconverged orbital energies"'
     $                 //' from chk',nbf,eigval,0,' ')
         else
            call iosys('read real "scf vector" from chk',nnbf,c,0,' ')
            call iosys('read real "orbital energies" from chk',
     $                 nbf,eigval,0,' ')
         end if
c
         call iosys('close chk',namchk,0,0,' ')
      else
         temp=file
         call lnkerr('unrecognized file name in rdges:'//temp)
      endif
c
c     possibly alter the guess.
      if(altges) call alter(c,eigval,scrtch,nbasis)
c
c
      return
      end
