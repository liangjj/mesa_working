*deck chksn4
      subroutine chksn4 (mbdcnd, nbdcnd, alpha, beta, cofx, singlr)
c***begin prologue  chksn4
c***subsidiary
c***purpose  subsidiary to sepx4
c***library   slatec
c***type      single precision (chksn4-s)
c***author  (unknown)
c***description
c
c     this subroutine checks if the pde sepx4
c     must solve is a singular operator.
c
c***see also  sepx4
c***routines called  (none)
c***common blocks    spl4
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  chksn4
c
      common /spl4/   kswx       ,kswy       ,k          ,l          ,
     1                ait        ,bit        ,cit        ,dit        ,
     2                mit        ,nit        ,is         ,ms         ,
     3                js         ,ns         ,dlx        ,dly        ,
     4                tdlx3      ,tdly3      ,dlx4       ,dly4
      logical         singlr
      external cofx
c***first executable statement  chksn4
      singlr = .false.
c
c     check if the boundary conditions are
c     entirely periodic and/or mixed
c
      if ((mbdcnd.ne.0 .and. mbdcnd.ne.3) .or.
     1    (nbdcnd.ne.0 .and. nbdcnd.ne.3)) return
c
c     check that mixed conditions are pure neuman
c
      if (mbdcnd .ne. 3) go to  10
      if (alpha.ne.0.0 .or. beta.ne.0.0) return
   10 continue
c
c     check that non-derivative coefficient functions
c     are zero
c
      do  30 i=is,ms
         xi = ait+(i-1)*dlx
         call cofx (xi,ai,bi,ci)
         if (ci .ne. 0.0) return
   30 continue
c
c     the operator must be singular if this point is reached
c
      singlr = .true.
      return
      end
