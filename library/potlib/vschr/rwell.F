*deck rwell.f
c***begin prologue     rwell
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            rounded potential well
c***                   
c***references         
c
c***routines called    
c***end prologue       rwell
      subroutine rwell(v,pt,nwell,awell,n,prn)
      implicit integer (a-z)
      real*8 v, pt, awell, one
      logical prn
      character*80 title
      data one/1.d0/
      dimension v(*), pt(*)
      common/io/inp, iout
      do 10 i=1,n
         v(i) = v(i) - one/sqrt( awell + pt(i)**nwell )
 10   continue   
      if(prn) then
         title='potential for this region'
         call prntrm(title,v,n,1,n,1,iout)
      endif
      return
      end       



