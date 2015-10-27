*deck vanhmo.f
c***begin prologue     vanhmo
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            anharmonic oscillator
c***                   
c***references         
c***routines called    
c***end prologue       vhmo
      subroutine vanhmo(v,pt,n,prn)
      implicit integer (a-z)
      real*8 v, pt, ptsq
      logical prn
      character*80 title
      dimension v(*), pt(*)
      common/io/inp, iout
      do 10 i=1,n
         ptsq=pt(i)*pt(i)
         v(i) = v(i) + .5d0*( ptsq + ptsq*ptsq )
 10   continue   
      if(prn) then
         title='potential'
         call prntfm(title,v,n,1,n,1,iout)  
      endif
      return
      end       



