*deck %W%  %G%
      subroutine pm105(z,a)
c
c***begin prologue     m105
c***date written       850601   (yymmdd)
c***revision date      920417   (yymmdd)
c
c   17 april    1992   rlm at lanl
c        modifying to read the auxiliary basis set information.
c   28 july     1988   bhl at llnl
c        passing sympt and mompt to genbas so we can handle linear molecules
c
c   16 november 1987   pws at lanl
c        adding call to 'grppt' to form the pointer and
c        other associated arrays needed to reorder the
c        two-particle density matrix for the gradient codes
c        for ci or mcscf gradients.
c
c    2 december 1986   pws at lanl
c        changing 'namdat' to a character  variable.
c
c***keywords           m105, link 105, basis, ecp, $basis, input,
c                      print
c***author             saxe, paul and martin, richard    (lanl)
c***source             %W% %G% 
c***purpose            prepares the basis set information.
c***description
c     m105 currently recognizes the option strings:
c       print_basis   print the basis set information with normalized
c                     primitives.
c       print_ecp     print the effective core potential information.
c       basis=string  the default basis set is denoted by the string.
c       q=n           molecular charge.
c       2s+1=n        molecular spin multiplicity..
c       timing        print timing statistics for this link.
c
c***references
c
c***routines called
c     m105
c       drum(mdutil)
c       iosys(io)
c       traktm(mdutil)
c       fillel(util)
c       getscm(mdutil)
c       genbas(local)
c       chainx(mdutil)
c
c***end prologue       m105
      implicit integer (a-z)
      integer a(*)
      real*8 z(*)
c
c
      parameter (maxnbf=2000,maxatm=2000)
c
      character*2 atsymb(0:104),basnam(maxatm)*16
      character*16 atname(maxatm)
      character*16 bflabl(maxnbf)
      character*4096 ops
      character*128 namdat
      integer sympt(maxnbf),mompt(maxnbf)
      logical logkey,usesym
      logical aux
c
      common /io/     inp,iout
c
      data aux/.false./
c
c     ----- get the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- blank out bflabl -----
c
      do 321 i=1,maxnbf
         bflabl(i)='                '
 321  continue
c
c     ----- get dimensions needed -----
c
      call iosys('read character "data filename" from rwf',
     $     0,0,0,namdat)
      call iosys('read integer "number of atoms" from rwf',
     $     1,natoms,0,' ')
c
c     ----- get atomic symbols -----
c
      call fillel(0,104,atsymb)
c
c     ----- are we to read in an auxiliary basis set -----
c           if so, the auxiliary bases are prefixed with "aux-"
      aux=logkey(ops,'scf=aux',.false.,' ')
c
c     ----- divvy up core, get it, and then read in basis type info -----
c
      atomz=1
      atomno=wpadti(atomz+natoms)
      top=atomno+natoms
c
      call getscm(top+120000,a(1),ngot,'m102 main',0)
c
      call genbas(a,z,top,ngot,atsymb,a(atomno),a(atomz),basnam,
     #            atname,natoms,namdat,inp,iout,bflabl,
     $            sympt,ops,mompt)
      if(aux) then
c
c     ----- blank out bflabl again -----
c
         do 322 i=1,maxnbf
            bflabl(i)='                '
 322     continue
         call genaux(a,z,top,ngot,atsymb,a(atomno),a(atomz),basnam,
     #               atname,natoms,namdat,inp,iout,bflabl,
     $               sympt,ops,mompt)
      endif
c
c     ----- re-allocate core, then work out arrays for gradients -----
c
      call iosys('read integer "number basis types" from rwf',
     $     1,ntypes,0,' ')
c
      ncont=1
      nobf=ncont+natoms*ntypes
      need=nobf+ntypes
c
      call getscm(need,a,ngot,'derivative',0)
c
      call iosys('read integer "number of contraction coefficients"'//
     $     ' from rwf',natoms*ntypes,a(ncont),0,' ')
      call iosys('read integer "number of pure functions" from rwf',
     $     ntypes,a(nobf),0,' ')
c
      call grppt(a(ncont),a(nobf),natoms,ntypes,a(need),ops)
c
c     ----- write a symmetry flag to the rwf -----
      usesym=.not.logkey(ops,'sym=off',.false.,' ')
      call iosys('write integer usesym to rwf',1,usesym,0,' ')
c
c     ----- exit gracefully -----
c
      call chainx(0)
c
c
      stop
      end
