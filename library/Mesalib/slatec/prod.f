*deck prod
      subroutine prod (nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a,
     +   b, c, d, w, u)
c***begin prologue  prod
c***subsidiary
c***purpose  subsidiary to blktri
c***library   slatec
c***type      single precision (prod-s, proc-c)
c***author  (unknown)
c***description
c
c prod applies a sequence of matrix operations to the vector x and
c stores the result in y.
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
c***end prologue  prod
c
      dimension       a(*)       ,b(*)       ,c(*)       ,x(*)       ,
     1                y(*)       ,d(*)       ,w(*)       ,bd(*)      ,
     2                bm1(*)     ,bm2(*)     ,aa(*)      ,u(*)
c***first executable statement  prod
      do 101 j=1,m
         w(j) = x(j)
         y(j) = w(j)
  101 continue
      mm = m-1
      id = nd
      ibr = 0
      m1 = nm1
      m2 = nm2
      ia = na
  102 if (ia) 105,105,103
  103 rt = aa(ia)
      if (nd .eq. 0) rt = -rt
      ia = ia-1
c
c scalar multiplication
c
      do 104 j=1,m
         y(j) = rt*w(j)
  104 continue
  105 if (id) 125,125,106
  106 rt = bd(id)
      id = id-1
      if (id .eq. 0) ibr = 1
c
c begin solution to system
c
      d(m) = a(m)/(b(m)-rt)
      w(m) = y(m)/(b(m)-rt)
      do 107 j=2,mm
         k = m-j
         den = b(k+1)-rt-c(k+1)*d(k+2)
         d(k+1) = a(k+1)/den
         w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
  107 continue
      den = b(1)-rt-c(1)*d(2)
      w(1) = 1.
      if (den) 108,109,108
  108 w(1) = (y(1)-c(1)*w(2))/den
  109 do 110 j=2,m
         w(j) = w(j)-d(j)*w(j-1)
  110 continue
      if (na) 113,113,102
  111 do 112 j=1,m
         y(j) = w(j)
  112 continue
      ibr = 1
      go to 102
  113 if (m1) 114,114,115
  114 if (m2) 111,111,120
  115 if (m2) 117,117,116
  116 if (abs(bm1(m1))-abs(bm2(m2))) 120,120,117
  117 if (ibr) 118,118,119
  118 if (abs(bm1(m1)-bd(id))-abs(bm1(m1)-rt)) 111,119,119
  119 rt = rt-bm1(m1)
      m1 = m1-1
      go to 123
  120 if (ibr) 121,121,122
  121 if (abs(bm2(m2)-bd(id))-abs(bm2(m2)-rt)) 111,122,122
  122 rt = rt-bm2(m2)
      m2 = m2-1
  123 do 124 j=1,m
         y(j) = y(j)+rt*w(j)
  124 continue
      go to 102
  125 return
      end
