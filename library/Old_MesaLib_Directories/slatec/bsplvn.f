*deck bsplvn
      subroutine bsplvn (t, jhigh, index, x, ileft, vnikx)
c***begin prologue  bsplvn
c***subsidiary
c***purpose  subsidiary to fc
c***library   slatec
c***type      single precision (bsplvn-s, dfspvn-d)
c***author  (unknown)
c***description
c
c calculates the value of all possibly nonzero b-splines at *x* of
c  order max(jhigh,(j+1)(index-1)) on *t*.
c
c***see also  fc
c***routines called  (none)
c***revision history  (yymmdd)
c   780801  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  bsplvn
      dimension t(*),vnikx(*)
      dimension deltam(20),deltap(20)
      save j, deltam, deltap
      data j/1/,(deltam(i),i=1,20),(deltap(i),i=1,20)/40*0./
c***first executable statement  bsplvn
                                       go to (10,20),index
   10 j = 1
      vnikx(1) = 1.
      if (j .ge. jhigh)                go to 99
c
   20    ipj = ileft+j
         deltap(j) = t(ipj) - x
         imjp1 = ileft-j+1
         deltam(j) = x - t(imjp1)
         vmprev = 0.
         jp1 = j+1
         do 26 l=1,j
            jp1ml = jp1-l
            vm = vnikx(l)/(deltap(l) + deltam(jp1ml))
            vnikx(l) = vm*deltap(l) + vmprev
   26       vmprev = vm*deltam(jp1ml)
         vnikx(jp1) = vmprev
         j = jp1
         if (j .lt. jhigh)             go to 20
c
   99                                  return
      end
