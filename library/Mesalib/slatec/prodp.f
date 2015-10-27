*deck prodp
      subroutine prodp (nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a,
     +   b, c, d, u, w)
c***begin prologue  prodp
c***subsidiary
c***purpose  subsidiary to blktri
c***library   slatec
c***type      single precision (prodp-s, procp-c)
c***author  (unknown)
c***description
c
c prodp applies a sequence of matrix operations to the vector x and
c stores the result in y (periodic boundary conditions).
c
c bd,bm1,bm2 are arrays containing roots of certain b polynomials.
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively.
c aa         array containing scalar multipliers of the vector x.
c na         is the length of the array aa.
c x,y        the matrix operations are applied to x and the result is y.
c a,b,c      are arrays which contain the tridiagonal matrix.
c m          is the order of the matrix.
c d,w,u      are working arrays.
c is         determines whether or not a change in sign is made.
c
c***see also  blktri
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  prodp
c
      dimension       a(*)       ,b(*)       ,c(*)       ,x(*)       ,
     1                y(*)       ,d(*)       ,u(*)       ,bd(*)      ,
     2                bm1(*)     ,bm2(*)     ,aa(*)      ,w(*)
c***first executable statement  prodp
      do 101 j=1,m
         y(j) = x(j)
         w(j) = y(j)
  101 continue
      mm = m-1
      mm2 = m-2
      id = nd
      ibr = 0
      m1 = nm1
      m2 = nm2
      ia = na
  102 if (ia) 105,105,103
  103 rt = aa(ia)
      if (nd .eq. 0) rt = -rt
      ia = ia-1
      do 104 j=1,m
         y(j) = rt*w(j)
  104 continue
  105 if (id) 128,128,106
  106 rt = bd(id)
      id = id-1
      if (id .eq. 0) ibr = 1
c
c begin solution to system
c
      bh = b(m)-rt
      ym = y(m)
      den = b(1)-rt
      d(1) = c(1)/den
      u(1) = a(1)/den
      w(1) = y(1)/den
      v = c(m)
      if (mm2-2) 109,107,107
  107 do 108 j=2,mm2
         den = b(j)-rt-a(j)*d(j-1)
         d(j) = c(j)/den
         u(j) = -a(j)*u(j-1)/den
         w(j) = (y(j)-a(j)*w(j-1))/den
         bh = bh-v*u(j-1)
         ym = ym-v*w(j-1)
         v = -v*d(j-1)
  108 continue
  109 den = b(m-1)-rt-a(m-1)*d(m-2)
      d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
      w(m-1) = (y(m-1)-a(m-1)*w(m-2))/den
      am = a(m)-v*d(m-2)
      bh = bh-v*u(m-2)
      ym = ym-v*w(m-2)
      den = bh-am*d(m-1)
      if (den) 110,111,110
  110 w(m) = (ym-am*w(m-1))/den
      go to 112
  111 w(m) = 1.
  112 w(m-1) = w(m-1)-d(m-1)*w(m)
      do 113 j=2,mm
         k = m-j
         w(k) = w(k)-d(k)*w(k+1)-u(k)*w(m)
  113 continue
      if (na) 116,116,102
  114 do 115 j=1,m
         y(j) = w(j)
  115 continue
      ibr = 1
      go to 102
  116 if (m1) 117,117,118
  117 if (m2) 114,114,123
  118 if (m2) 120,120,119
  119 if (abs(bm1(m1))-abs(bm2(m2))) 123,123,120
  120 if (ibr) 121,121,122
  121 if (abs(bm1(m1)-bd(id))-abs(bm1(m1)-rt)) 114,122,122
  122 rt = rt-bm1(m1)
      m1 = m1-1
      go to 126
  123 if (ibr) 124,124,125
  124 if (abs(bm2(m2)-bd(id))-abs(bm2(m2)-rt)) 114,125,125
  125 rt = rt-bm2(m2)
      m2 = m2-1
  126 do 127 j=1,m
         y(j) = y(j)+rt*w(j)
  127 continue
      go to 102
  128 return
      end
