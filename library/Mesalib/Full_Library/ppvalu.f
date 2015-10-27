*deck ppvalu
      function ppvalu (break, coef, l, k, x, jderiv )
calculates value at  x  of  jderiv-th derivative of spline from pp-repr
      implicit real *8 (a-h,o-z)
      save
      dimension break(l), coef(k,l)
      ppvalu = 0.d0
      fmmjdr = k - jderiv
c  derivatives of order  k  or higher are identically zero.
      if (fmmjdr .le. 0.d0)               go to 99
c  find index i of largest breakpoint to the left of  x .
      call interv ( break, l, x, i, ndummy )
c  evaluate  jderiv -th derivative of  i -th polynomial piece at  x .
      h = x - break(i)
      lim=k-jderiv
      do 10 mm=1,lim
         m=k-mm+1
         ppvalu = (ppvalu/fmmjdr)*h + coef(m,i)
   10    fmmjdr = fmmjdr - 1.d0
   99                                   return
      end
