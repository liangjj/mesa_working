*deck vadex12.f
c***begin prologue     vadex12
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            two particle exponential
c***                   
c***references         
c
c***routines called    
c***end prologue       vadex12
      subroutine vadex12(v,q1,q2,a12,b12,n,prn)
      implicit integer (a-z)
      real*8 v, q1, q2, a12, b12
      logical prn
      character*80 title
      dimension v(n), q1(n), q2(n)
      common/io/inp, iout
      do 10 i=1,n
         v(i) = v(i) + a12 * exp ( - b12 * ( q2(i) + q1(i) ) )
 10   continue   
      if(prn) then
         title='two dimensional exponential potential'
         call prntrm(title,v,n,1,n,1,iout)
      endif
      return
      end       



