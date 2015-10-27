*deck intrv.f
c***begin prologue     intrv
c***date written       930623   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           newton-raphson root finder
c***author             schneider, barry (nsf)
c***source             math
c***purpose            locate interval in which a function goes through zero
c***                   
c***                   
c
c***references         
c
c***routines called    
c***end prologue       intrv
      subroutine intrv(fun,zval,rl,del,n,m,cntzro,prnt)
c
      implicit integer (a-z)
      real*8 fun, zval, rl, del, val, x, prd
      logical prnt
      common/io/inp, iout 
      dimension fun(n), zval(m+1)
c
      cntzro=0
      zval(1)=rl
      x=rl
      val=fun(1)
      do 10 i=2,n
         if (cntzro.eq.m) go to 100
         x=x+del
         prd=val*fun(i)
         if (prd.lt.0.d0) then
              cntzro=cntzro+1
              zval(cntzro+1)=x
              val=fun(i)
         endif
   10 continue
  100 write(iout,1) cntzro
    1 format(/,5x,'there are ',i3,' zeros in the input interval')
      if (prnt) then
          write(iout,2) ( zval(i),i=1,cntzro+1) 
      endif
    2 format(/,5x,'the intervals for the zeros',(/,5(1x,e15.8)))
      return
      end       
