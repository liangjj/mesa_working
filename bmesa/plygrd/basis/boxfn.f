*deck boxfn
      subroutine boxfn(sn,cn,x,left,right,n)
c***begin prologue     boxfn
c***date written       940504   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose             
c***                   .
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       boxfn
c
      implicit integer (a-z)
      real*8 sn, cn, x, left, right
      real*8 pi, prefac, norm, arg 
      character*80 title
      dimension sn(n,n), cn(n,n), x(n)
      common /io/ inp, iout
      data pi/3.1415926535897932384d0/
      prefac = pi/(right-left)
      norm=sqrt(2.d0/(right-left))
      do 20 i=1,n
         do 30 j=1,n
            arg = i*prefac*( x(j) - left )
            sn(j,i) = sin(arg )
            cn(j,i) = cos(arg)
 30      continue
 20   continue
      call vscale(sn,sn,norm,n*n)
      call vscale(cn,cn,norm,n*n)
      return
      end















