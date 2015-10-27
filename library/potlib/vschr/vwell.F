*deck vwell.f
c***begin prologue     vwell
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
c***end prologue       vwell
      subroutine vwell(v,pt,len,depth,n,prn)
      implicit integer (a-z)
      real*8 v, pt, depth, len
      logical prn
      character*80 title
      dimension v(*), pt(*)
      common/io/inp, iout
      do 10 i=1,n
         if(pt(i).le.len) then
            v(i) = v(i) + depth
         endif            
 10   continue   
      if(prn) then
         title='potential'
         call prntrm(title,v,n,1,n,1,iout)
      endif
      return
      end       



