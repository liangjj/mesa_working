*deck dfspvd
      subroutine dfspvd (t, k, x, ileft, vnikx, nderiv)
c***begin prologue  dfspvd
c***subsidiary
c***purpose  subsidiary to dfc
c***library   slatec
c***type      double precision (bsplvd-s, dfspvd-d)
c***author  (unknown)
c***description
c
c   **** double precision version of bsplvd ****
c calculates value and deriv.s of all b-splines which do not vanish at x
c
c  fill vnikx(j,ideriv), j=ideriv, ... ,k  with nonzero values of
c  b-splines of order k+1-ideriv , ideriv=nderiv, ... ,1, by repeated
c  calls to dfspvn
c
c***see also  dfc
c***routines called  dfspvn
c***revision history  (yymmdd)
c   780801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  dfspvd
      implicit double precision (a-h,o-z)
      dimension t(*),vnikx(k,*)
      dimension a(20,20)
c***first executable statement  dfspvd
      call dfspvn(t,k+1-nderiv,1,x,ileft,vnikx(nderiv,nderiv))
      if (nderiv .le. 1)               go to 99
      ideriv = nderiv
      do 15 i=2,nderiv
         idervm = ideriv-1
         do 11 j=ideriv,k
   11       vnikx(j-1,idervm) = vnikx(j,ideriv)
         ideriv = idervm
         call dfspvn(t,0,2,x,ileft,vnikx(ideriv,ideriv))
   15    continue
c
      do 20 i=1,k
         do 19 j=1,k
   19       a(i,j) = 0.d0
   20    a(i,i) = 1.d0
      kmd = k
      do 40 m=2,nderiv
         kmd = kmd-1
         fkmd = kmd
         i = ileft
         j = k
   21       jm1 = j-1
            ipkmd = i + kmd
            diff = t(ipkmd) - t(i)
            if (jm1 .eq. 0)            go to 26
            if (diff .eq. 0.d0)          go to 25
            do 24 l=1,j
   24          a(l,j) = (a(l,j) - a(l,j-1))/diff*fkmd
   25       j = jm1
            i = i - 1
                                       go to 21
   26    if (diff .eq. 0.)             go to 30
         a(1,1) = a(1,1)/diff*fkmd
c
   30    do 40 i=1,k
            v = 0.d0
            jlow = max(i,m)
            do 35 j=jlow,k
   35          v = a(i,j)*vnikx(j,m) + v
   40       vnikx(i,m) = v
   99                                  return
      end
