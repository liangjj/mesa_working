*deck @(#)testnt.f	1.1 9/8/91
      subroutine testnt(ints,mmax)
      implicit integer (a-z)
      real*8 pi, step, val, ints
      character*80 title
      dimension ints(0:mmax,2)
      common/io/inp,iout
      data pi / 3.141592653589793d0 /
      do 40 m=0,mmax
         nmin=max(m+1,3)
         n=2*(nmin/2)+1
         n=9
         step=2.d0*pi/(n-1)
         val=step
         do 50 i=2,n-1,2
            ints(m,1)=ints(m,1)+4.d0*sin(m*val)
            ints(m,2)=ints(m,2)+4.d0*cos(m*val)
            val=val+step+step
   50    continue
         val=step+step
         do 60 i=3,n-2,2
            ints(m,1)=ints(m,1)+2.d0*sin(m*val)
            ints(m,2)=ints(m,2)+2.d0*cos(m*val)
            val=val+step+step
   60    continue
         ints(m,2)=ints(m,2)+2.d0
         ints(m,1)=step*ints(m,1)/3.d0
         ints(m,2)=step*ints(m,2)/3.d0
   40 continue
      title='trigonometric sine and cosine integrals'
      call prntfm(title,ints,mmax+1,2,mmax+1,2,iout)
      return
      end


