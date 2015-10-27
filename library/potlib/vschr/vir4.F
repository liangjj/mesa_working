*deck vir4.f
c***begin prologue     vir4
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
c***end prologue       vir4
      subroutine vir4(v,pt,n,prn)
      implicit integer (a-z)
      real*8 v, pt, one
      logical prn
      character*80 title
      data one/1.d0/
      dimension v(*), pt(*)
      common/io/inp, iout
      do 10 i=1,n
	 v(i) = v(i) + one/( one + pt(i) )**4
 10   continue   
      if(prn) then
         title='potential'
         call prntfm(title,v,n,1,n,1,iout)  
      endif
      return
      end       



