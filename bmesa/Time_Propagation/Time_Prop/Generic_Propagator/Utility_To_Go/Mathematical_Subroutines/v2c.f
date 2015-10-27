*deck v2c.f	
      subroutine v2c(vr,vi,v,n,num)
c***begin prologue     v2c
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, complex
c***author             schneider barry(nsf)
c***source             
c***purpose            v = vr + i*vi
c***description
c                      call v2c(vr,vi,v,n,num)
c                        v        complex output vector of length n.
c                        vr       real input vector of length n.
c                        vi       real input vector of length n.
c                        n        vector lengths.
c                        num      number of vectors
c
c***references
c***routines called    (none)
c***end prologue       v2c
      real*8 vr(n,num), vi(n,num)
      complex*16 v(n,num), eye
      data eye /(0.d0,1.d0)/
      do 1 i=1,num
         do 2 j=1,n
            v(j,i)= vr(j,i) + eye*vi(j,i)
    2    continue
    1 continue
      return
      end
