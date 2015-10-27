*deck vres.f
c***begin prologue     vres
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            resonant exponential potential
c***                   
c***references         
c
c***routines called    
c***end prologue       vres
      subroutine vres(v,pt,a,s,shift,n,prn)
      implicit integer (a-z)
      real*8 v, pt, a, s, shift
      logical prn
      character*80 title
      dimension v(*), pt(*), a(2), s(2)
      common/io/inp, iout
      do 10 i=1,n
         v(i) = v(i) + 
     1          a(1)*exp(-s(1)*pt(i)) + 
     2          a(2)*exp(-s(2)*( pt(i) - shift )*( pt(i) - shift ) ) 
 10   continue   
      if(prn) then
         title='potential'
         call prntfm(title,v,n,1,n,1,iout)  
      endif
      return
      end       



