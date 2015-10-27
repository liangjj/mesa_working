*deck chksng
      subroutine chksng (mbdcnd, nbdcnd, alpha, beta, gama, xnu, cofx,
     +   cofy, singlr)
c***begin prologue  chksng
c***subsidiary
c***purpose  subsidiary to sepeli
c***library   slatec
c***type      single precision (chksng-s)
c***author  (unknown)
c***description
c
c     this subroutine checks if the pde sepeli
c     must solve is a singular operator.
c
c***see also  sepeli
c***routines called  (none)
c***common blocks    splpcm
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  chksng
c
      common /splpcm/ kswx       ,kswy       ,k          ,l          ,
     1                ait        ,bit        ,cit        ,dit        ,
     2                mit        ,nit        ,is         ,ms         ,
     3                js         ,ns         ,dlx        ,dly        ,
     4                tdlx3      ,tdly3      ,dlx4       ,dly4
      logical         singlr
c***first executable statement  chksng
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
   10 if (nbdcnd .ne. 3) go to  20
      if (gama.ne.0.0 .or. xnu.ne.0.0) return
   20 continue
c
c     check that non-derivative coefficient functions
c     are zero
c
      do  30 i=is,ms
         xi = ait+(i-1)*dlx
         call cofx (xi,ai,bi,ci)
         if (ci .ne. 0.0) return
   30 continue
      do  40 j=js,ns
         yj = cit+(j-1)*dly
         call cofy (yj,dj,ej,fj)
         if (fj .ne. 0.0) return
   40 continue
c
c     the operator must be singular if this point is reached
c
      singlr = .true.
      return
      end
