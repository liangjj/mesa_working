*deck rv2cv.f	
      subroutine rv2cv(vr,vc,n,num)
c***begin prologue     rv2cv
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, complex
c***author             schneider barry(nsf)
c***source             
c***purpose            vc(i) = vr(i,1) + i*vi(i,2)
c***description
c                      call rv2cv(vr,vc,n,num)
c                        vc       complex output vector of length n.
c                        vr       real input vector of length n.
c                        n        vector lengths.
c                        num      number of vectors
c
c***references
c***routines called    (none)
c***end prologue       rv2cv
      real*8 vr(n,2,num)
      complex*16 vc(n,num), eye
      data eye /(0.d0,1.d0)/
      common /io/ inp,iout
      do 1 i=1,num
         do 2 j=1,n
            vc(j,i)= vr(j,1,i) + eye*vr(j,2,i)
    2    continue
    1 continue
      return
      end
