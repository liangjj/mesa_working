*deck addang.f
c***begin prologue     addang
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            centrifugal potential
c***                   
c***references         
c
c***routines called    
c***end prologue       addang
      subroutine addang(v,pt,l,coord,n,prn)
      implicit integer (a-z)
      real*8 v, pt
      logical prn
      character*80 title
      character*(*) coord
      character*3 itoc
      dimension v(n), pt(n)
      common/io/inp, iout
      lnum=l*(l+1)
      if(coord.eq.'rho') then
         lnum=l*l
      endif
      do 10 i=1,n
         v(i) = v(i) + lnum/(2.d0*pt(i)*pt(i))
 10   continue   
      if(prn) then
         title='potential l = '//itoc(l)
         call prntfm(title,v,n,1,n,1,iout)
      endif
      return
      end       



