*deck factl
c***begin prologue     factl
c***date written       920405   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           factorials
c***author             schneider, barry (nsf)
c***source             mylib
c***purpose            calculate factorials from zero to n
c***description        
c***                   
c***                   
c
c***references         
c
c***routines called    
c***end prologue      factl
      subroutine factl(fact,n)
      implicit integer (a-z)
      real *8 fact
      dimension fact(0:n)
      data maxn /100/
      common/io/ inp, iout
      if(n.gt.100) then
         call lnkerr('factorial will overflow')
      endif
      fact(0)=1.d0
      if (n.eq.0) then
          return
      endif
      do 10 i=1,n
         fact(i) = i * fact(i-1) 
   10 continue
      return
      end
