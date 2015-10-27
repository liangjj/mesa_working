*deck @(#)phyold.f	5.1  11/6/94
      subroutine phyold(n,phycon)
c***begin prologue     phyold.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           physical constants
c***author             binkley, et. al., gaussian 82.
c***source             @(#)phyold.f	5.1   11/6/94
c***purpose            loads the values of several physical constants.
c***description
c     call phyold(n,phycon)
c       n       the dimension of phycon.
c       phycon  receives the values of the physical constants.
c
c***references         pure and applied chemistry, 51, 1 (1979).
c***routines called    (none)
c***end prologue       phyold.f
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
      real*8 one,four
      real*8 toang,tokg,toe,planck,avog,jpcal,tomet,hartre
      real*8 slight,boltz,finei
c
      parameter (one=1.0d+00,four=4.0d+00)
c
c     --- these data were taken from:
c     pure and applied chemistry, 51, 1 (1979) unless noted otherwise.
c
c                                         angstroms per bohr
      data toang  /0.52917706 d+00/
c                                         kilograms per atomic mass unit
      data tokg   /1.6605655  d-27/
c                                         electrostatic units (esu)
c                                            per electron charge
c                                            pure appl. chem.,
c                                            2, 717 (1973)
      data toe    /4.803242   d-10/
c                                         planck constant, joule-seconds
      data planck /6.626176   d-34/
c                                         avogadro constant
      data avog   /6.022045   d+23/
c                                         joules per calorie
      data jpcal  /4.184      d+00/
c                                         metres per bohr
      data tomet  /5.2917706  d-11/
c                                         joules per hartree
      data hartre /4.359814   d-18/
c                                         speed of light, cm sec(-1)
      data slight /2.99792458 d+10/
c                                         boltzman constant, joules per
c                                            kelvin
      data boltz  /1.380662   d-23/
c
c                                         inverse fine structure constant.
c                                         rev. mod. phys. 41 375 (1969).
      data finei  /137.03602  d+00/
      save toang,tokg,toe,planck,avog,jpcal,tomet,hartre
      save slight,boltz,finei
c
      common/io/inp,iout
c
c     --- load the constants
      call rzero(phycon,n)
      phycon( 1) = toang
      phycon( 2) = tokg
      phycon( 3) = toe
      phycon( 4) = planck
      phycon( 5) = avog
      phycon( 6) = jpcal
      phycon( 7) = tomet
      phycon( 8) = hartre
      phycon( 9) = slight
      phycon(10) = boltz
      phycon(11) = one / finei
c     --- pi.
      phycon(12)=four*atan(one)
c
c
      return
      end
