*deck preph.f
c***begin prologue     preph
c***date written       970503   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           precondition
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***routines called    
c***end prologue       preph
      subroutine preph(ham,v,d,n)
      implicit integer (a-z)
      real*8 ham, v, d
      dimension ham(n,n), v(n,n), d(n)
      common/io/inp, iout
      do 10 i=1,n
         d(i) = ham(i,i) +v(i,i)
   10 continue
      do 20 i=1,n
         ham(i,i) = 0.d0
 20   continue   
      return
      end       
