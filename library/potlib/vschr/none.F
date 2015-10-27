*deck none.f
c***begin prologue     none
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            no potential
c***                   
c***references         
c
c***routines called    
c***end prologue       none
      subroutine none(v,pt,n,prn)
      implicit integer (a-z)
      real*8 v, pt
      logical prn
      character*80 title
      dimension v(*), pt(*)
      common/io/inp, iout
      if(prn) then
         title='potential'
         call prntrm(title,v,n,1,n,1,iout)
      endif
      return
      end       



