*deck vwad12.f
c***begin prologue     vwad12
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            two dimensional potential well in adiabatic
c***                   coordinates
c***                   
c***references         
c
c***routines called    
c***end prologue       vwad12
      subroutine vwad12(v,q1,q2,depth,n,prn)
      implicit integer (a-z)
      real*8 v, q1, q2, depth
      logical prn
      character*80 title
      dimension v(n), q1(n), q2(n)
      common/io/inp, iout
      write(iout,1) depth
      do 10 i=1,n
         v(i) = v(i) + depth
 10   continue   
      if(prn) then
         title='two dimensional well'
         call prntrm(title,v,n,1,n,1,iout)
      endif
      return
 1    format(/,1x,'two dimensional potential well depth = ',e15.8)
      end       



