*deck vexp12.f
c***begin prologue     vexp12
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            exponential pair potential
c***                   
c***references         
c
c***routines called    
c***end prologue       vexp12
      subroutine vexp12(v,q1,q2,a12,b12,n1,n2,prn)
      implicit integer (a-z)
      real*8 v, q1, q2, a12, b12
      logical prn
      character*80 title
      dimension v(n2,n1), q1(n1), q2(n2)
      common/io/inp, iout
      do 10 i=1,n1
         do 20 j=1,n2
            v(j,i) = v(j,i) + a12 * exp ( - b12 * abs( q2(j) - q1(i) ) )
 20      continue   
 10   continue   
      if(prn) then
         title='two dimensional exponential interaction'
         call prntrm(title,v,n2,n1,n2,n1,iout)
      endif
      return
      end       



