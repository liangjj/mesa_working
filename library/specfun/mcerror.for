	program mcerror
c
c       ============================================================
c       purpose: this program computes the error function erf(z) 
c                for a complex argument using subroutine cerror
c       input :  x   --- real part of z
c                y   --- imaginary part of z  ( y ó 3.0 )
c       output:  err --- real part of erf(z)
c                eri --- imaginary part of erf(z)
c       example:
c                   x       y       re[erf(z)]      im[erf(z)]
c                 ---------------------------------------------
c                  1.0     2.0      -.53664357     -5.04914370
c                  2.0     2.0      1.15131087       .12729163
c                  3.0     2.0       .99896328      -.00001155
c                  4.0     2.0      1.00000057      -.00000051
c                  5.0     2.0      1.00000000       .00000000
c       ============================================================
c
	implicit complex *16 (c,z)  
	double precision x,y
	write(*,*)'x,y=?'
	read(*,*)x,y
	write(*,*)'   x      y      re[erf(z)]      im[erf(z)]'
	write(*,*)' ---------------------------------------------'
	z=cmplx(x,y)
	call cerror(z,cer)
	write(*,10) z,cer
10      format(1x,f5.1,2x,f5.1,1x,2e16.8)
	end


	subroutine cerror(z,cer)
c
c       ====================================================
c       purpose: compute error function erf(z) for a complex
c                argument (z=x+iy)
c       input :  z   --- complex argument
c       output:  cer --- erf(z)
c       ====================================================
c
	implicit complex *16 (c,z)
	double precision a0,pi
	a0=cdabs(z)
	c0=cdexp(-z*z)
	pi=3.141592653589793d0
	z1=z
	if (real(z).lt.0.0) then
	   z1=-z
	endif
	if (a0.le.5.8d0) then    
	   cs=z1
	   cr=z1
	   do 10 k=1,120
	      cr=cr*z1*z1/(k+0.5d0)
	      cs=cs+cr
	      if (cdabs(cr/cs).lt.1.0d-15) go to 15
10         continue
15         cer=2.0d0*c0*cs/dsqrt(pi)
	else                              
	   cl=1.0d0/z1              
	   cr=cl
	   do 20 k=1,13
	      cr=-cr*(k-0.5d0)/(z1*z1)
	      cl=cl+cr
	      if (cdabs(cr/cl).lt.1.0d-15) go to 25
20         continue
25         cer=1.0d0-c0*cl/dsqrt(pi)
	endif
	if (real(z).lt.0.0) then
	   cer=-cer
	endif
	return
	end


