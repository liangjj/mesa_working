*deck @(#)cheby.f	1.1 9/8/91
      subroutine cheby(n,mmax)
      implicit integer (a-z)
      character*80 title
      real*8 pi, integ, val, step, fac1, fac2, fac3
      dimension integ(0:mmax,2)
      data pi /3.141592653d0/
      common /io/ inp,iout
      fac1=1.d0/3.d0
      fac2=4.d0/3.d0
      fac3=2.d0/3.d0
      nq=max(n,3)
      nq=nq/2
      nq=2*nq+1
      do 10 m=0,mmax
         integ(m,1)=0.d0
         integ(m,2)=0.d0
         do 20 i=1,nq
            val=(2*i-1)*pi/(2*nq)
            integ(m,1)=integ(m,1)+ sin(m*val)
            integ(m,2)=integ(m,2)+ cos(m*val)
   20    continue
         integ(m,1)=pi*integ(m,1)/nq
         integ(m,2)=pi*integ(m,2)/nq  
   10 continue
      title='trigonometric sine and cosine integrals-cheby'
      call prntfm(title,integ,mmax+1,2,mmax+1,2,iout)      
      step=pi/(nq-1)
      write(iout,*) ' step = ',step
      do 30 m=0,mmax
         integ(m,2)=fac1
         integ(m,1)=0.d0
         val=step
         do 40 i=2,nq-1,2
            integ(m,1)=integ(m,1)+ fac2*sin(m*val)
            integ(m,2)=integ(m,2)+ fac2*cos(m*val)
            val=val+step+step
   40    continue
         val=step+step
         do 50 i=3,nq-2,2
            integ(m,1)=integ(m,1)+fac3*sin(m*val)
            integ(m,2)=integ(m,2)+fac3*cos(m*val)
            val=val+step+step
   50    continue
         integ(m,1)=integ(m,1)+fac1*sin(m*pi)
         integ(m,2)=integ(m,2)+fac1*cos(m*pi)
   30 continue
      call vscale(integ(0,1),integ(0,1),step,mmax+1)
      call vscale(integ(0,2),integ(0,2),step,mmax+1)
      title='trigonometric sine and cosine integrals-simpson'
      call prntfm(title,integ,mmax+1,2,mmax+1,2,iout)      
      return
      end
