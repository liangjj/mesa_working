*deck ccmpb
      subroutine ccmpb (n, ierror, an, bn, cn, b, ah, bh)
c***begin prologue  ccmpb
c***subsidiary
c***purpose  subsidiary to cblktr
c***library   slatec
c***type      complex (compb-s, ccmpb-c)
c***author  (unknown)
c***description
c
c     ccmpb computes the roots of the b polynomials using subroutine
c     tevlc which is a modification the eispack program tqlrat.
c     ierror is set to 4 if either tevlc fails or if a(j+1)*c(j) is
c     less than zero for some j.  ah,bh are temporary work arrays.
c
c***see also  cblktr
c***routines called  cpadd, inxcb, r1mach, tevlc
c***common blocks    ccblk
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  ccmpb
c
      dimension       an(*)      ,bn(*)      ,cn(*)      ,b(*)       ,
     1                ah(*)      ,bh(*)
      common /ccblk/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
c***first executable statement  ccmpb
      eps = r1mach(4)
      bnorm = abs(bn(1))
      do 102 j=2,nm
         bnorm = max(bnorm,abs(bn(j)))
         arg = an(j)*cn(j-1)
         if (arg) 119,101,101
  101    b(j) = sign(sqrt(arg),an(j))
  102 continue
      cnv = eps*bnorm
      if = 2**k
      kdo = k-1
      do 108 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i4 = i2+i2
         ipl = i4-1
         ifd = if-i4
         do 107 i=i4,ifd,i4
            call inxcb (i,l,ib,nb)
            if (nb) 108,108,103
  103       js = i-ipl
            jf = js+nb-1
            ls = 0
            do 104 j=js,jf
               ls = ls+1
               bh(ls) = bn(j)
               ah(ls) = b(j)
  104       continue
            call tevlc (nb,bh,ah,ierror)
            if (ierror) 118,105,118
  105       lh = ib-1
            do 106 j=1,nb
               lh = lh+1
               b(lh) = -bh(j)
  106       continue
  107    continue
  108 continue
      do 109 j=1,nm
         b(j) = -bn(j)
  109 continue
      if (npp) 117,110,117
  110 nmp = nm+1
      nb = nm+nmp
      do 112 j=1,nb
         l1 = mod(j-1,nmp)+1
         l2 = mod(j+nm-1,nmp)+1
         arg = an(l1)*cn(l2)
         if (arg) 119,111,111
  111    bh(j) = sign(sqrt(arg),-an(l1))
         ah(j) = -bn(l1)
  112 continue
      call tevlc (nb,ah,bh,ierror)
      if (ierror) 118,113,118
  113 call inxcb (if,k-1,j2,lh)
      call inxcb (if/2,k-1,j1,lh)
      j2 = j2+1
      lh = j2
      n2m2 = j2+nm+nm-2
  114 d1 = abs(b(j1)-b(j2-1))
      d2 = abs(b(j1)-b(j2))
      d3 = abs(b(j1)-b(j2+1))
      if ((d2 .lt. d1) .and. (d2 .lt. d3)) go to 115
      b(lh) = b(j2)
      j2 = j2+1
      lh = lh+1
      if (j2-n2m2) 114,114,116
  115 j2 = j2+1
      j1 = j1+1
      if (j2-n2m2) 114,114,116
  116 b(lh) = b(n2m2+1)
      call inxcb (if,k-1,j1,j2)
      j2 = j1+nmp+nmp
      call cpadd (nm+1,ierror,an,cn,b(j1),b(j1),b(j2))
  117 return
  118 ierror = 4
      return
  119 ierror = 5
      return
      end
