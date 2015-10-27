*deck xpsi
      real function xpsi (a, ipsik, ipsix)
c***begin prologue  xpsi
c***subsidiary
c***purpose  to compute values of the psi function for xlegf.
c***library   slatec
c***category  c7c
c***type      single precision (xpsi-s, dxpsi-d)
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
c***end prologue  xpsi
      real a,b,c,cnum,cdenom
      dimension cnum(12),cdenom(12)
      save cnum, cdenom
c
c        cnum(i) and cdenom(i) are the ( reduced ) numerator
c        and 2*i*denominator respectively of the 2*i th bernoulli
c        number.
c
      data cnum(1),cnum(2),cnum(3),cnum(4),cnum(5),cnum(6),cnum(7),
     1cnum(8),cnum(9),cnum(10),cnum(11),cnum(12)
     2    / 1.,     -1.,    1.,     -1., 1.,
     3   -691.,  1.,     -3617., 43867., -174611., 77683.,
     4   -236364091./
      data cdenom(1),cdenom(2),cdenom(3),cdenom(4),cdenom(5),cdenom(6),
     1 cdenom(7),cdenom(8),cdenom(9),cdenom(10),cdenom(11),cdenom(12)
     2/12.,120.,   252.,   240.,132.,
     3  32760., 12.,  8160., 14364., 6600., 276., 65520./
c***first executable statement  xpsi
      n=max(0,ipsix-int(a))
      b=n+a
      k1=ipsik-1
c
c        series expansion for a .gt. ipsix using ipsik-1 terms.
c
      c=0.
      do 12 i=1,k1
      k=ipsik-i
   12 c=(c+cnum(k)/cdenom(k))/b**2
      xpsi=log(b)-(c+.5/b)
      if(n.eq.0) go to 20
      b=0.
c
c        recurrence for a .le. ipsix.
c
      do 15 m=1,n
   15 b=b+1./(n-m+a)
      xpsi=xpsi-b
   20 return
      end
