*deck @(#)phyfil.f	5.1 11/6/94
      subroutine phyfil(n,phycon)
c***begin prologue     phyfil.f
c***date written       921016  yymmdd  
c***revision date      11/6/94
c
c   february 2, 1993   rlm at lanl
c      replacing the definition of the constants (one,four)
c      with a data statement, rather than a parameter statement
c      because for some unknown reason the sun stores the constants
c      in double precision with the former, and only single precision
c      in the latter -- a very strange compiler bug. this results in
c      the value of pi being slightly off, which propagates into
c      slight errors in the integrals.
c 
c***keywords           physical constants
c***author             frisch, et al., gaussian92 
c***source             @(#)phyfil.f	5.1   11/6/94
c***purpose            loads the values of several physical constants 
c***description
c     call phyfil(n,phycon)
c       n       the dimension of phycon.
c       phycon  receives the values of the physical constants.
c
c***references
c
c     these data were taken from tables in:
c
c     physics today, august 1987, page bg11, which in turn came from:
c
c     e. r. cohen, b. n. taylor, "the 1986 adjustment of the
c     fundamental physical constants," report of the codata task
c     group on fundamental constants, codata bulletin 63, pergamon,
c     elmsford, ny (1986).
c
c***routines called
c
c***end prologue       phyfil.f
      implicit none
c     --- input variables -----
      integer n
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 phycon(n)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      real*8 one,four,atan
      real*8 toang,tokg,toe,planck,avog,jpcal,hartre,slight,boltz
      real*8 finei,volmol
c
      data one    /1.0d+00/
      data four   /4.0d+00/
c
c                                         angstroms per bohr
      data toang  /0.529177249 d 00/
c                                         kilograms per atomic mass unit
      data tokg   /1.6605402  d-27/
c                                         electrostatic units (esu)
c                                            per electron charge
c                                            pure appl. chem.,
c                                            2, 717 (1973)
      data toe    /4.803242   d-10/
c                                         planck constant, joule-seconds
      data planck /6.6260755   d-34/
c                                         avogadro constant
      data avog   /6.0221367   d+23/
c                                         joules per calorie
      data jpcal  /4.184      d 00/
c                                         joules per hartre
      data hartre /4.3597482   d-18/
c                                         speed of light, cm sec(-1)
      data slight /2.99792458 d+10/
c                                         boltzman constant, joules per
c                                            kelvin
      data boltz  /1.380658   d-23/
c
c                                         inverse fine structure constant.
      data finei  /137.0359895  d+00/
c
c     molar volume of ideal gas in m**3 at 273.15 k
c     j. phys. chem. ref. data 2 717 (1973).
      data volmol/22.41410d-3/
      save toang, tokg, toe, planck, avog, jpcal, hartre, slight, boltz
      save finei, volmol
c
      common/io/inp,iout
c
c     --- load the array
      call rzero(phycon,n)
      phycon( 1) = toang
      phycon( 2) = tokg
      phycon( 3) = toe
      phycon( 4) = planck
      phycon( 5) = avog
      phycon( 6) = jpcal
c
c     --- meters per bohr.
      phycon( 7) = toang / float(10)**10
      phycon( 8) = hartre
      phycon( 9) = slight
      phycon(10) = boltz
      phycon(11) = one / finei
c
c     --- pi
      phycon(12)=four*atan(one)
c
c     --- electron mass in kg.
      phycon(13) = hartre*finei*finei*float(10000) / (slight*slight)
      phycon(14) = volmol
c
c
      return
      end
