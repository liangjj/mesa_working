*deck cpadd
      subroutine cpadd (n, ierror, a, c, cbp, bp, bh)
c***begin prologue  cpadd
c***subsidiary
c***purpose  subsidiary to cblktr
c***library   slatec
c***type      single precision (cpadd-s)
c***author  (unknown)
c***description
c
c   cpadd computes the eigenvalues of the periodic tridiagonal matrix
c   with coefficients an,bn,cn.
c
c   n    is the order of the bh and bp polynomials.
c   bp   contains the eigenvalues on output.
c   cbp  is the same as bp except type complex.
c   bh   is used to temporarily store the roots of the b hat polynomial
c        which enters through bp.
c
c***see also  cblktr
c***routines called  bcrh, pgsf, ppgsf, pppsf
c***common blocks    ccblk
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  cpadd
c
      complex         cx         ,fsg        ,hsg        ,
     1                dd         ,f          ,fp         ,fpp        ,
     2                cdis       ,r1         ,r2         ,r3         ,
     3                cbp
      dimension       a(*)       ,c(*)       ,bp(*)      ,bh(*)      ,
     1                cbp(*)
      common /ccblk/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
      external        pgsf       ,pppsf      ,ppgsf
c***first executable statement  cpadd
      scnv = sqrt(cnv)
      iz = n
      if (bp(n)-bp(1)) 101,142,103
  101 do 102 j=1,n
         nt = n-j
         bh(j) = bp(nt+1)
  102 continue
      go to 105
  103 do 104 j=1,n
         bh(j) = bp(j)
  104 continue
  105 ncmplx = 0
      modiz = mod(iz,2)
      is = 1
      if (modiz) 106,107,106
  106 if (a(1)) 110,142,107
  107 xl = bh(1)
      db = bh(3)-bh(1)
  108 xl = xl-db
      if (pgsf(xl,iz,c,a,bh)) 108,108,109
  109 sgn = -1.
      cbp(1) = cmplx(bcrh(xl,bh(1),iz,c,a,bh,pgsf,sgn),0.)
      is = 2
  110 if = iz-1
      if (modiz) 111,112,111
  111 if (a(1)) 112,142,115
  112 xr = bh(iz)
      db = bh(iz)-bh(iz-2)
  113 xr = xr+db
      if (pgsf(xr,iz,c,a,bh)) 113,114,114
  114 sgn = 1.
      cbp(iz) = cmplx(bcrh(bh(iz),xr,iz,c,a,bh,pgsf,sgn),0.)
      if = iz-2
  115 do 136 ig=is,if,2
         xl = bh(ig)
         xr = bh(ig+1)
         sgn = -1.
         xm = bcrh(xl,xr,iz,c,a,bh,pppsf,sgn)
         psg = pgsf(xm,iz,c,a,bh)
         if (abs(psg)-eps) 118,118,116
  116    if (psg*ppgsf(xm,iz,c,a,bh)) 117,118,119
c
c     case of a real zero
c
  117    sgn = 1.
         cbp(ig) = cmplx(bcrh(bh(ig),xm,iz,c,a,bh,pgsf,sgn),0.)
         sgn = -1.
         cbp(ig+1) = cmplx(bcrh(xm,bh(ig+1),iz,c,a,bh,pgsf,sgn),0.)
         go to 136
c
c     case of a multiple zero
c
  118    cbp(ig) = cmplx(xm,0.)
         cbp(ig+1) = cmplx(xm,0.)
         go to 136
c
c     case of a complex zero
c
  119    it = 0
         icv = 0
         cx = cmplx(xm,0.)
  120    fsg = (1.,0.)
         hsg = (1.,0.)
         fp = (0.,0.)
         fpp = (0.,0.)
         do 121 j=1,iz
            dd = 1./(cx-bh(j))
            fsg = fsg*a(j)*dd
            hsg = hsg*c(j)*dd
            fp = fp+dd
            fpp = fpp-dd*dd
  121    continue
         if (modiz) 123,122,123
  122    f = (1.,0.)-fsg-hsg
         go to 124
  123    f = (1.,0.)+fsg+hsg
  124    i3 = 0
         if (abs(fp)) 126,126,125
  125    i3 = 1
         r3 = -f/fp
  126    if (abs(fpp)) 132,132,127
  127    cdis = sqrt(fp**2-2.*f*fpp)
         r1 = cdis-fp
         r2 = -fp-cdis
         if (abs(r1)-abs(r2)) 129,129,128
  128    r1 = r1/fpp
         go to 130
  129    r1 = r2/fpp
  130    r2 = 2.*f/fpp/r1
         if (abs(r2) .lt. abs(r1)) r1 = r2
         if (i3) 133,133,131
  131    if (abs(r3) .lt. abs(r1)) r1 = r3
         go to 133
  132    r1 = r3
  133    cx = cx+r1
         it = it+1
         if (it .gt. 50) go to 142
         if (abs(r1) .gt. scnv) go to 120
         if (icv) 134,134,135
  134    icv = 1
         go to 120
  135    cbp(ig) = cx
         cbp(ig+1) = conjg(cx)
  136 continue
      if (abs(cbp(n))-abs(cbp(1))) 137,142,139
  137 nhalf = n/2
      do 138 j=1,nhalf
         nt = n-j
         cx = cbp(j)
         cbp(j) = cbp(nt+1)
         cbp(nt+1) = cx
  138 continue
  139 ncmplx = 1
      do 140 j=2,iz
         if (aimag(cbp(j))) 143,140,143
  140 continue
      ncmplx = 0
      do 141 j=2,iz
         bp(j) = real(cbp(j))
  141 continue
      go to 143
  142 ierror = 4
  143 continue
      return
      end
