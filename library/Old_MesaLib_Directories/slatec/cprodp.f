*deck cprodp
      subroutine cprodp (nd, bd, nm1, bm1, nm2, bm2, na, aa, x, yy, m,
     +   a, b, c, d, u, y)
c***begin prologue  cprodp
c***subsidiary
c***purpose  subsidiary to blktri
c***library   slatec
c***type      single precision (cprodp-s, cprocp-c)
c***author  (unknown)
c***description
c
c prodp applies a sequence of matrix operations to the vector x and
c stores the result in yy. (periodic boundary conditions and complex
c case)
c
c bd,bm1,bm2     are arrays containing roots of certain b polynomials.
c nd,nm1,nm2     are the lengths of the arrays bd,bm1,bm2 respectively.
c aa             array containing scalar multipliers of the vector x.
c na             is the length of the array aa.
c x,yy      the matrix operations are applied to x and the result is yy.
c a,b,c          are arrays which contain the tridiagonal matrix.
c m              is the order of the matrix.
c d,u,y          are working arrays.
c isgn           determines whether or not a change in sign is made.
c
c***see also  blktri
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  cprodp
c
      complex         y          ,d          ,u          ,v          ,
     1                den        ,bh         ,ym         ,am         ,
     2                y1         ,y2         ,yh         ,bd         ,
     3                crt
      dimension       a(*)       ,b(*)       ,c(*)       ,x(*)       ,
     1                y(*)       ,d(*)       ,u(*)       ,bd(*)      ,
     2                bm1(*)     ,bm2(*)     ,aa(*)      ,yy(*)
c***first executable statement  cprodp
      do 101 j=1,m
         y(j) = cmplx(x(j),0.)
  101 continue
      mm = m-1
      mm2 = m-2
      id = nd
      m1 = nm1
      m2 = nm2
      ia = na
  102 iflg = 0
      if (id) 111,111,103
  103 crt = bd(id)
      id = id-1
      iflg = 1
c
c begin solution to system
c
      bh = b(m)-crt
      ym = y(m)
      den = b(1)-crt
      d(1) = c(1)/den
      u(1) = a(1)/den
      y(1) = y(1)/den
      v = cmplx(c(m),0.)
      if (mm2-2) 106,104,104
  104 do 105 j=2,mm2
         den = b(j)-crt-a(j)*d(j-1)
         d(j) = c(j)/den
         u(j) = -a(j)*u(j-1)/den
         y(j) = (y(j)-a(j)*y(j-1))/den
         bh = bh-v*u(j-1)
         ym = ym-v*y(j-1)
         v = -v*d(j-1)
  105 continue
  106 den = b(m-1)-crt-a(m-1)*d(m-2)
      d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
      y(m-1) = (y(m-1)-a(m-1)*y(m-2))/den
      am = a(m)-v*d(m-2)
      bh = bh-v*u(m-2)
      ym = ym-v*y(m-2)
      den = bh-am*d(m-1)
      if (abs(den)) 107,108,107
  107 y(m) = (ym-am*y(m-1))/den
      go to 109
  108 y(m) = (1.,0.)
  109 y(m-1) = y(m-1)-d(m-1)*y(m)
      do 110 j=2,mm
         k = m-j
         y(k) = y(k)-d(k)*y(k+1)-u(k)*y(m)
  110 continue
  111 if (m1) 112,112,114
  112 if (m2) 123,123,113
  113 rt = bm2(m2)
      m2 = m2-1
      go to 119
  114 if (m2) 115,115,116
  115 rt = bm1(m1)
      m1 = m1-1
      go to 119
  116 if (abs(bm1(m1))-abs(bm2(m2))) 118,118,117
  117 rt = bm1(m1)
      m1 = m1-1
      go to 119
  118 rt = bm2(m2)
      m2 = m2-1
c
c matrix multiplication
c
  119 yh = y(1)
      y1 = (b(1)-rt)*y(1)+c(1)*y(2)+a(1)*y(m)
      if (mm-2) 122,120,120
  120 do 121 j=2,mm
         y2 = a(j)*y(j-1)+(b(j)-rt)*y(j)+c(j)*y(j+1)
         y(j-1) = y1
         y1 = y2
  121 continue
  122 y(m) = a(m)*y(m-1)+(b(m)-rt)*y(m)+c(m)*yh
      y(m-1) = y1
      iflg = 1
      go to 102
  123 if (ia) 126,126,124
  124 rt = aa(ia)
      ia = ia-1
      iflg = 1
c
c scalar multiplication
c
      do 125 j=1,m
         y(j) = rt*y(j)
  125 continue
  126 if (iflg) 127,127,102
  127 do 128 j=1,m
         yy(j) = real(y(j))
  128 continue
      return
      end
