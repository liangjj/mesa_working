*deck cmptr3
      subroutine cmptr3 (m, a, b, c, k, y1, y2, y3, tcos, d, w1, w2, w3)
c***begin prologue  cmptr3
c***subsidiary
c***purpose  subsidiary to cmgnbn
c***library   slatec
c***type      complex (tri3-s, cmptr3-c)
c***author  (unknown)
c***description
c
c     subroutine to solve tridiagonal systems.
c
c***see also  cmgnbn
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   890206  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  cmptr3
      complex         a          ,b          ,c          ,y1         ,
     1                y2         ,y3         ,tcos       ,d          ,
     2                w1         ,w2         ,w3         ,x          ,
     3                xx         ,z
      dimension       a(*)       ,b(*)       ,c(*)       ,k(4)       ,
     1                tcos(*)    ,y1(*)      ,y2(*)      ,y3(*)      ,
     2                d(*)       ,w1(*)      ,w2(*)      ,w3(*)
      integer k1p1, k2p1, k3p1, k4p1
c
c***first executable statement  cmptr3
      mm1 = m-1
      k1 = k(1)
      k2 = k(2)
      k3 = k(3)
      k4 = k(4)
      k1p1 = k1+1
      k2p1 = k2+1
      k3p1 = k3+1
      k4p1 = k4+1
      k2k3k4 = k2+k3+k4
      if (k2k3k4 .eq. 0) go to 101
      l1 = k1p1/k2p1
      l2 = k1p1/k3p1
      l3 = k1p1/k4p1
      lint1 = 1
      lint2 = 1
      lint3 = 1
      kint1 = k1
      kint2 = kint1+k2
      kint3 = kint2+k3
  101 continue
      do 115 n=1,k1
         x = tcos(n)
         if (k2k3k4 .eq. 0) go to 107
         if (n .ne. l1) go to 103
         do 102 i=1,m
            w1(i) = y1(i)
  102    continue
  103    if (n .ne. l2) go to 105
         do 104 i=1,m
            w2(i) = y2(i)
  104    continue
  105    if (n .ne. l3) go to 107
         do 106 i=1,m
            w3(i) = y3(i)
  106    continue
  107    continue
         z = 1./(b(1)-x)
         d(1) = c(1)*z
         y1(1) = y1(1)*z
         y2(1) = y2(1)*z
         y3(1) = y3(1)*z
         do 108 i=2,m
            z = 1./(b(i)-x-a(i)*d(i-1))
            d(i) = c(i)*z
            y1(i) = (y1(i)-a(i)*y1(i-1))*z
            y2(i) = (y2(i)-a(i)*y2(i-1))*z
            y3(i) = (y3(i)-a(i)*y3(i-1))*z
  108    continue
         do 109 ip=1,mm1
            i = m-ip
            y1(i) = y1(i)-d(i)*y1(i+1)
            y2(i) = y2(i)-d(i)*y2(i+1)
            y3(i) = y3(i)-d(i)*y3(i+1)
  109    continue
         if (k2k3k4 .eq. 0) go to 115
         if (n .ne. l1) go to 111
         i = lint1+kint1
         xx = x-tcos(i)
         do 110 i=1,m
            y1(i) = xx*y1(i)+w1(i)
  110    continue
         lint1 = lint1+1
         l1 = (lint1*k1p1)/k2p1
  111    if (n .ne. l2) go to 113
         i = lint2+kint2
         xx = x-tcos(i)
         do 112 i=1,m
            y2(i) = xx*y2(i)+w2(i)
  112    continue
         lint2 = lint2+1
         l2 = (lint2*k1p1)/k3p1
  113    if (n .ne. l3) go to 115
         i = lint3+kint3
         xx = x-tcos(i)
         do 114 i=1,m
            y3(i) = xx*y3(i)+w3(i)
  114    continue
         lint3 = lint3+1
         l3 = (lint3*k1p1)/k4p1
  115 continue
      return
      end
