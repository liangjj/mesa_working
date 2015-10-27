*deck dmacon
      subroutine dmacon
c***begin prologue  dmacon
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (macon-s, dmacon-d)
c***author  (unknown)
c***see also  dbvsup
c***routines called  d1mach
c***common blocks    dml5mc
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  dmacon
      double precision d1mach
      integer ke, lpar
      double precision dd, eps, fouru, sqovfl, sru, twou, uro
      common /dml5mc/ uro,sru,eps,sqovfl,twou,fouru,lpar
c***first executable statement  dmacon
      uro = d1mach(4)
      sru = sqrt(uro)
      dd = -log10(uro)
      lpar = 0.5d0*dd
      ke = 0.5d0 + 0.75d0*dd
      eps = 10.0d0**(-2*ke)
      sqovfl = sqrt(d1mach(2))
      twou = 2.0d0*uro
      fouru = 4.0d0*uro
      return
      end
