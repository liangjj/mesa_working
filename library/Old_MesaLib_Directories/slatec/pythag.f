*deck pythag
      real function pythag (a, b)
c***begin prologue  pythag
c***subsidiary
c***purpose  compute the complex square root of a complex number without
c            destructive overflow or underflow.
c***library   slatec
c***type      single precision (pythag-s)
c***author  (unknown)
c***description
c
c     finds sqrt(a**2+b**2) without overflow or destructive underflow
c
c***see also  eisdoc
c***routines called  (none)
c***revision history  (yymmdd)
c   811101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  pythag
      real a,b
c
      real p,q,r,s,t
c***first executable statement  pythag
      p = max(abs(a),abs(b))
      q = min(abs(a),abs(b))
      if (q .eq. 0.0e0) go to 20
   10 continue
         r = (q/p)**2
         t = 4.0e0 + r
         if (t .eq. 4.0e0) go to 20
         s = r/t
         p = p + 2.0e0*p*s
         q = q*s
      go to 10
   20 pythag = p
      return
      end
