*deck vcad12.f
c***begin prologue     vcad12
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            one-d model coulomb interaction
c***                   
c***references         
c
c***routines called    
c***end prologue       vcad12
      subroutine vcad12(v,x,y,z,n,prn)
      implicit integer (a-z)
      real*8 v, x, y, z
      logical prn
      character*80 title
      dimension v(*), x(*), y(*)
      common/io/inp, iout
      write(iout,1) z
      do 10 i=1,n
         v(i) = v(i) + z / (x(i) + y(i))
 10   continue   
      if(prn) then
         title='potential'
         call prntrm(title,v,n,1,n,1,iout)
      endif
      return
 1    format(/,1x,'potential = z / ( x1 + x2 ) where',/,1x,
     1            '         z = ',e15.8)
      end       



