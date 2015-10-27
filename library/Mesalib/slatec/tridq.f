*deck tridq
      subroutine tridq (mr, a, b, c, y, d)
c***begin prologue  tridq
c***subsidiary
c***purpose  subsidiary to pois3d
c***library   slatec
c***type      single precision (tridq-s)
c***author  (unknown)
c***see also  pois3d
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900308  renamed routine from trid to tridq.  (wrb)
c   900402  added type section.  (wrb)
c***end prologue  tridq
      dimension       a(*)       ,b(*)       ,c(*)       ,y(*)       ,
     1                d(*)
c***first executable statement  tridq
      m = mr
      mm1 = m-1
      z = 1./b(1)
      d(1) = c(1)*z
      y(1) = y(1)*z
      do 101 i=2,mm1
         z = 1./(b(i)-a(i)*d(i-1))
         d(i) = c(i)*z
         y(i) = (y(i)-a(i)*y(i-1))*z
  101 continue
      z = b(m)-a(m)*d(mm1)
      if (z .ne. 0.) go to 102
      y(m) = 0.
      go to 103
  102 y(m) = (y(m)-a(m)*y(mm1))/z
  103 continue
      do 104 ip=1,mm1
         i = m-ip
         y(i) = y(i)-d(i)*y(i+1)
  104 continue
      return
      end
