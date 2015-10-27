*deck vhmo.f
c***begin prologue     vhmo
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            potential well
c***                   
c***references         
c
c***routines called    
c***end prologue       vhmo
      subroutine vhmo(v,pt,fac,n,prn)
      implicit integer (a-z)
      real*8 v, pt, fac
      logical prn
      character*80 title
      dimension v(*), pt(*)
      common/io/inp, iout
      do 10 i=1,n
         v(i) = v(i) + fac*pt(i)*pt(i)
 10   continue   
      if(prn) then
         title='potential'
         call prntfm(title,v,n,1,n,1,iout)  
      endif
      return
      end       



