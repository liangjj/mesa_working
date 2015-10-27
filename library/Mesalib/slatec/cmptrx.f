*deck cmptrx
      subroutine cmptrx (idegbr, idegcr, m, a, b, c, y, tcos, d, w)
c***begin prologue  cmptrx
c***subsidiary
c***purpose  subsidiary to cmgnbn
c***library   slatec
c***type      complex (trix-s, cmptrx-c)
c***author  (unknown)
c***description
c
c     subroutine to solve a system of linear equations where the
c     coefficient matrix is a rational function in the matrix given by
c     tridiagonal  ( . . . , a(i), b(i), c(i), . . . ).
c
c***see also  cmgnbn
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  cmptrx
c
      complex         a          ,b          ,c          ,y          ,
     1                tcos       ,d          ,w          ,x          ,
     2                xx         ,z
      dimension       a(*)       ,b(*)       ,c(*)       ,y(*)       ,
     1                tcos(*)    ,d(*)       ,w(*)
      integer kb, kc
c***first executable statement  cmptrx
      mm1 = m-1
      kb = idegbr+1
      kc = idegcr+1
      l = kb/kc
      lint = 1
      do 108 k=1,idegbr
         x = tcos(k)
         if (k .ne. l) go to 102
         i = idegbr+lint
         xx = x-tcos(i)
         do 101 i=1,m
            w(i) = y(i)
            y(i) = xx*y(i)
  101    continue
  102    continue
         z = 1./(b(1)-x)
         d(1) = c(1)*z
         y(1) = y(1)*z
         do 103 i=2,mm1
            z = 1./(b(i)-x-a(i)*d(i-1))
            d(i) = c(i)*z
            y(i) = (y(i)-a(i)*y(i-1))*z
  103    continue
         z = b(m)-x-a(m)*d(mm1)
         if (abs(z) .ne. 0.) go to 104
         y(m) = (0.,0.)
         go to 105
  104    y(m) = (y(m)-a(m)*y(mm1))/z
  105    continue
         do 106 ip=1,mm1
            i = m-ip
            y(i) = y(i)-d(i)*y(i+1)
  106    continue
         if (k .ne. l) go to 108
         do 107 i=1,m
            y(i) = y(i)+w(i)
  107    continue
         lint = lint+1
         l = (lint*kb)/kc
  108 continue
      return
      end
