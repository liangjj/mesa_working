*deck dxpsi
      double precision function dxpsi (a, ipsik, ipsix)
c***begin prologue  dxpsi
c***subsidiary
c***purpose  to compute values of the psi function for dxlegf.
c***library   slatec
c***category  c7c
c***type      double precision (xpsi-s, dxpsi-d)
c***keywords  psi function
c***author  smith, john m., (nbs and george mason university)
c***routines called  (none)
c***revision history  (yymmdd)
c   820728  date written
c   890126  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  dxpsi
      double precision a,b,c,cnum,cdenom
      dimension cnum(12),cdenom(12)
      save cnum, cdenom
c
c        cnum(i) and cdenom(i) are the ( reduced ) numerator
c        and 2*i*denominator respectively of the 2*i th bernoulli
c        number.
c
      data cnum(1),cnum(2),cnum(3),cnum(4),cnum(5),cnum(6),cnum(7),
     1cnum(8),cnum(9),cnum(10),cnum(11),cnum(12)
     2    / 1.d0,     -1.d0,    1.d0,     -1.d0, 1.d0,
     3   -691.d0,  1.d0,     -3617.d0, 43867.d0, -174611.d0, 77683.d0,
     4   -236364091.d0/
      data cdenom(1),cdenom(2),cdenom(3),cdenom(4),cdenom(5),cdenom(6),
     1 cdenom(7),cdenom(8),cdenom(9),cdenom(10),cdenom(11),cdenom(12)
     2/12.d0,120.d0,   252.d0,   240.d0,132.d0,
     3  32760.d0, 12.d0,  8160.d0, 14364.d0, 6600.d0, 276.d0, 65520.d0/
c***first executable statement  dxpsi
      n=max(0,ipsix-int(a))
      b=n+a
      k1=ipsik-1
c
c        series expansion for a .gt. ipsix using ipsik-1 terms.
c
      c=0.d0
      do 12 i=1,k1
      k=ipsik-i
   12 c=(c+cnum(k)/cdenom(k))/b**2
      dxpsi=log(b)-(c+.5d0/b)
      if(n.eq.0) go to 20
      b=0.d0
c
c        recurrence for a .le. ipsix.
c
      do 15 m=1,n
   15 b=b+1.d0/(n-m+a)
      dxpsi=dxpsi-b
   20 return
      end
