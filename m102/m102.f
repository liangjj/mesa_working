
*deck @(#)pm102.f	5.1 11/6/94
      program m102
c***begin prologue     pm102.f
c***date written       850601   (yymmdd)
c***revision date      11/6/94
c
c    7 july     1993   rlm at lanl
c        removing obsolete symmetry stuff for linear molecules.
c   28 july     1988   bhl at llnl
c        passing sympt and mompt to genbas so we can handle linear molecules
c   16 november 1987   pws at lanl
c        adding call to 'grppt' to form the pointer and
c        other associated arrays needed to reorder the
c        two-particle density matrix for the gradient codes
c        for ci or mcscf gradients.
c    2 december 1986   pws at lanl
c        changing 'namdat' to a character  variable.
c
c***keywords           m102, link 102, basis, ecp, $basis, input,
c                      print
c***author             saxe, paul and martin, richard    (lanl)
c***source             @(#)pm102.f	5.1 11/6/94
c***purpose            prepares the basis set information.
c***description
c     m102 currently recognizes the option strings:
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
c***end prologue       pm102.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a
      real*8 z
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer maxnbf,maxatm
c
      parameter (maxnbf=2000,maxatm=2000)
c
      integer inp,iout
      integer i,atomz,atomno,ncont,top
      integer natoms,ngot,need,nobf,ntypes
      integer wpadti
      integer idum
      character*2 atsymb(0:104),basnam(maxatm)*16
      character*16 atname(maxatm)
      character*16 bflabl(maxnbf)
      character*4096 ops
      character*128 namdat
      logical logkey,usesym
c
      common /io/     inp,iout
      pointer (p,z(1)), (p,a(1))
c
      call drum
c     --- get the options string ---
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- blank out bflabl ---
      do 321 i=1,maxnbf
         bflabl(i)='                '
 321  continue
c
c     --- get dimensions needed ---
      call iosys('read character "data filename" from rwf',
     $           0,0,0,namdat)
      call iosys('read integer "number of atoms" from rwf',
     $           1,natoms,0,' ')
c
c     --- get atomic symbols ---
      call fillel(0,104,atsymb)
c
c     --- divvy up core, get it, and then read in basis type info ---
c
      atomz=1
      atomno=wpadti(atomz+natoms)
      top=atomno+natoms
c
c      call getscm(top+120000,a(1),ngot,'m102 main',0)
c      call getmem(top+120000,p,ngot,'m102 main',0)
      call getmem(top,p,ngot,'m102 main',0)
c
      call genbas(atsymb,a(atomno),a(atomz),basnam,
     $            atname,natoms,namdat,inp,iout,bflabl,ops)
      call getmem(-ngot,p,idum,'m102 main',idum)
c
c     --- re-allocate core, then work out arrays for gradients ---
      call iosys('read integer "number basis types" from rwf',
     $            1,ntypes,0,' ')
c
      ncont=1
      nobf=ncont+natoms*ntypes
      need=nobf+ntypes
c
c      call getscm(need,a,ngot,'derivative',0)
      call getmem(need,p,ngot,'derivative',0)
c
      call iosys('read integer "number of contraction coefficients"'//
     $           ' from rwf',natoms*ntypes,a(ncont),0,' ')
      call iosys('read integer "number of pure functions" from rwf',
     $            ntypes,a(nobf),0,' ')
c
      call grppt(a(ncont),a(nobf),natoms,ntypes,ops)
      call getmem(-ngot,p,idum,'derivative',idum)
c
c     --- write a symmetry flag to the rwf ---
      usesym=.not.logkey(ops,'sym=off',.false.,' ')
      call iosys('write integer usesym to rwf',1,usesym,0,' ')
c
c     --- exit gracefully ---
      call chainx(0)
c
c
      stop
      end
