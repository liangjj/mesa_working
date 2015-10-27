*deck bsplvb
      subroutine bsplvb ( t, jhigh, index, x, ileft, biatx )
      implicit real *8 (a-h,o-z)
calculates the value of all possibly nonzero b-splines at  x  of
c  order max(jhigh,(j+1)(index-1)) on  t .
c     dimension (t(ileft+jhigh)
      save
      dimension t(1),biatx(jhigh), deltam(20),deltap(20)
      data j/1/
content of j, deltam, deltap is expected unchanged between calls.
                                        go to (10,20), index
   10 j = 1
      biatx(1) = 1.
      if (j .ge. jhigh)                 go to 99
c
   20    jp1 = j + 1
         deltap(j) = t(ileft+j) - x
         deltam(j) = x - t(ileft-j+1)
         saved = 0.d0
         do 26 i=1,j
            term = biatx(i)/(deltap(i) + deltam(jp1-i))
            biatx(i) = saved + deltap(i)*term
   26       saved = deltam(jp1-i)*term
         biatx(jp1) = saved
         j = jp1
         if (j .lt. jhigh)              go to 20
c
   99                                   return
      end
