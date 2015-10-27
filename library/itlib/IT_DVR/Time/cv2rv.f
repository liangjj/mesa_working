*deck cv2rv.f	
      subroutine cv2rv(vr,vc,n,num)
c***begin prologue     cv2rv
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, complex
c***author             schneider barry(nsf)
c***source             
c***purpose            vr = real(vc)  vi =imag(v)
c***description
c                      call cv2rv(vr,vc,n,num)
c                        v        complex input vector of length n.
c                        vr       real output vector of length n.
c                        vi       real output vector of length n.
c                        n        vector lengths.
c                        m        column lengths.
c                        num      number of vectors
c
c***references
c***routines called    (none)
c***end prologue       cv2rv
      real*8 vr(n,2,num)
      complex*16 vc(n,num)
      common /io/ inp,iout
      do 1 i=1,num
         do 2 j=1,n
            vr(j,1,i) = real(vc(j,i))
            vr(j,2,i) = imag(vc(j,i))
    2    continue
    1 continue
      return
      end
