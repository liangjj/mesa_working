*deck macon
      subroutine macon
c***begin prologue  macon
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (macon-s, dmacon-d)
c***author  (unknown)
c***description
c
c    sets up machine constants using r1mach
c
c***see also  bvsup
c***routines called  r1mach
c***common blocks    ml5mco
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  macon
      common /ml5mco/ uro,sru,eps,sqovfl,twou,fouru,lpar
c***first executable statement  macon
      uro=r1mach(4)
      sru=sqrt(uro)
      dd=-log10(uro)
      lpar=0.5*dd
      ke=0.5+0.75*dd
      eps=10.**(-2*ke)
      sqovfl=sqrt(r1mach(2))
      twou=2.0*uro
      fouru=4.0*uro
      return
      end
