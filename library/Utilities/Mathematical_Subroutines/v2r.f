*deck v2r.f	
      subroutine v2r(vr,vi,v,n,num)
c***begin prologue     v2r
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, complex
c***author             schneider barry(nsf)
c***source             
c***purpose            vr = real(v)  vi =imag(v)
c***description
c                      call v2r(vr,vi,v,n,num)
c                        v        complex input vector of length n.
c                        vr       real output vector of length n.
c                        vi       real output vector of length n.
c                        n        vector lengths.
c                        num      number of vectors
c
c***references
c***routines called    (none)
c***end prologue       v2r
      real*8 vr(n,num), vi(n,num)
      complex*16 v(n,num)
      do 1 i=1,num
         do 2 j=1,n
            vr(j,i) = real(v(j,i))
            vi(j,i) = imag(v(j,i))
    2    continue
    1 continue
      return
      end
