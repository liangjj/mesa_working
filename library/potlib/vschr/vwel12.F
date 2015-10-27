*deck vwel12.f
c***begin prologue     vwel12
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            two dimensional potential well
c***                   
c***references         
c
c***routines called    
c***end prologue       vwel12
      subroutine vwel12(v,q1,q2,a12,d12,n1,n2,prn)
      implicit integer (a-z)
      real*8 v, q1, q2, a12, d12
      logical prn
      character*80 title
      dimension v(n2,n1), q1(n1), q2(n2)
      common/io/inp, iout
      do 10 i=1,n1
         do 20 j=1,n2
            if( abs(q1(i)-q2(j)).le.a12) then
                v(j,i) = v(j,i) + d12
            endif
 20      continue   
 10   continue   
      if(prn) then
         title='two dimensional well'
         call prntrm(title,v,n2,n1,n2,n1,iout)
      endif
      return
 1    format(/,9x,'coordinate',8x,'potential')
 2    format(5x,e15.8,3x,e15.8)
      end       



