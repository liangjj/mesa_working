*deck cproc
      subroutine cproc (nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a,
     +   b, c, d, w, yy)
c***begin prologue  cproc
c***subsidiary
c***purpose  subsidiary to cblktr
c***library   slatec
c***type      complex (cprod-s, cproc-c)
c***author  (unknown)
c***description
c
c proc applies a sequence of matrix operations to the vector x and
c stores the result in y.
c aa     array containing scalar multipliers of the vector x.
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively.
c bd,bm1,bm2 are arrays containing roots of certain b polynomials.
c na     is the length of the array aa.
c x,y    the matrix operations are applied to x and the result is y.
c a,b,c  are arrays which contain the tridiagonal matrix.
c m      is the order of the matrix.
c d,w    are work arrays.
c isgn   determines whether or not a change in sign is made.
c
c***see also  cblktr
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  cproc
c
      complex         y          ,d          ,w          ,bd         ,
     1                crt        ,den        ,y1         ,y2         ,
     2                x          ,a          ,b          ,c
      dimension       a(*)       ,b(*)       ,c(*)       ,x(*)       ,
     1                y(*)       ,d(*)       ,w(*)       ,bd(*)      ,
     2                bm1(*)     ,bm2(*)     ,aa(*)      ,yy(*)
c***first executable statement  cproc
      do 101 j=1,m
         y(j) = x(j)
  101 continue
      mm = m-1
      id = nd
      m1 = nm1
      m2 = nm2
      ia = na
  102 iflg = 0
      if (id) 109,109,103
  103 crt = bd(id)
      id = id-1
c
c begin solution to system
c
      d(m) = a(m)/(b(m)-crt)
      w(m) = y(m)/(b(m)-crt)
      do 104 j=2,mm
         k = m-j
         den = b(k+1)-crt-c(k+1)*d(k+2)
         d(k+1) = a(k+1)/den
         w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
  104 continue
      den = b(1)-crt-c(1)*d(2)
      if (abs(den)) 105,106,105
  105 y(1) = (y(1)-c(1)*w(2))/den
      go to 107
  106 y(1) = (1.,0.)
  107 do 108 j=2,m
         y(j) = w(j)-d(j)*y(j-1)
  108 continue
  109 if (m1) 110,110,112
  110 if (m2) 121,121,111
  111 rt = bm2(m2)
      m2 = m2-1
      go to 117
  112 if (m2) 113,113,114
  113 rt = bm1(m1)
      m1 = m1-1
      go to 117
  114 if (abs(bm1(m1))-abs(bm2(m2))) 116,116,115
  115 rt = bm1(m1)
      m1 = m1-1
      go to 117
  116 rt = bm2(m2)
      m2 = m2-1
  117 y1 = (b(1)-rt)*y(1)+c(1)*y(2)
      if (mm-2) 120,118,118
c
c matrix multiplication
c
  118 do 119 j=2,mm
         y2 = a(j)*y(j-1)+(b(j)-rt)*y(j)+c(j)*y(j+1)
         y(j-1) = y1
         y1 = y2
  119 continue
  120 y(m) = a(m)*y(m-1)+(b(m)-rt)*y(m)
      y(m-1) = y1
      iflg = 1
      go to 102
  121 if (ia) 124,124,122
  122 rt = aa(ia)
      ia = ia-1
      iflg = 1
c
c scalar multiplication
c
      do 123 j=1,m
         y(j) = rt*y(j)
  123 continue
  124 if (iflg) 125,125,102
  125 return
      end
