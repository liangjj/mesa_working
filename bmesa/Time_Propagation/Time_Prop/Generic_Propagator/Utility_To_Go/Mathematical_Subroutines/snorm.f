*deck snorm
c***begin prologue     snorm
c***date written       861108   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m1104, link 1104, wighted norm
c***author             schneider, barry (lanl)
c***source             m1104
c***purpose            scalar product in r space.
c***description        norm of two vectors including integration weight
c***                   in definition of scalar product.
c
c***references         none
c
c***routines called
c***end prologue       snorm
      function snorm (f1,f2,wt,n,nowgt)
      real *8 f1, f2, wt, sum, snorm
      logical nowgt
      dimension f1(n), f2(n), wt(n)
      sum=0.d+00
      if (.not.nowgt) then
         do 10 i=1,n
            sum=sum+f1(i)*wt(i)*f2(i)
   10    continue
      else
          do 20 i=1,n
             sum=sum+f1(i)*f2(i)
   20     continue
      endif
      snorm=sum
      return
      end
