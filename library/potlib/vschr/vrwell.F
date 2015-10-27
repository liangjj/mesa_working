*deck vrwell.f
c***begin prologue     vrwell
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
c***end prologue       vrwell
      subroutine vrwell(v,pt,nwell,awell,n,prn)
      implicit integer (a-z)
      real*8 v, pt, spt, awell, one
      logical prn
      character*80 title
      data one/1.d0/
      dimension v(*), pt(*)
      common/io/inp, iout
      do 10 i=1,n
         v(i) = v(i) - one/sqrt( awell + pt(i)**nwell )
 10   continue   
      if(prn) then
         title='potential'
         call prntfm(title,v,n,1,n,1,iout)  
      endif
      return
      end       



