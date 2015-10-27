*deck dfactl
c***begin prologue     dfactl
c***date written       920405   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           double factorials
c***author             schneider, barry (nsf)
c***source             mylib
c***purpose            calculate double factorials from zero to n
c***description        
c***                   
c***                   
c
c***references         
c
c***routines called    
c***end prologue      factl
      subroutine dfactl(dfact,n)
      implicit integer (a-z)
      real *8 dfact
      dimension dfact(0:n)
      dfact(0)=1.d0
      dfact(1)=1.d0
      dfact(2)=3.d0
      if (n.gt.2) then
          do 10 i=3,n
             dfact(i) = ( i + i - 1) * dfact(i-1) 
   10     continue
      endif 
      return
      end
