*deck cpevlr
      subroutine cpevlr (n, m, a, x, c)
c***begin prologue  cpevlr
c***subsidiary
c***purpose  subsidiary to cpzero
c***library   slatec
c***type      single precision (cpevlr-s)
c***author  (unknown)
c***see also  cpzero
c***routines called  (none)
c***revision history  (yymmdd)
c   810223  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  cpevlr
      real a(*),c(*)
c***first executable statement  cpevlr
      np1=n+1
      do 1 j=1,np1
            ci=0.0
            cim1=a(j)
            mini=min(m+1,n+2-j)
            do 1 i=1,mini
               if(j .ne. 1) ci=c(i)
               if(i .ne. 1) cim1=c(i-1)
               c(i)=cim1+x*ci
    1 continue
      return
      end
