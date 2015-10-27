*deck srot
      subroutine srot(n,sx,incx,sy,incy,sc,ss)
      real*8  sx,sy,sc,ss,zero,one,w,z
      dimension sx(1),sy(1)
      data zero,one/0.d0,1.d0/
      if(n .le. 0 .or. (ss .eq. zero .and. sc .eq. one)) go to 40
      if(.not. (incx .eq. incy .and. incx .gt. 0)) go to 20
           nsteps=incx*n
           do 10 i=1,nsteps,incx
                w=sx(i)
                z=sy(i)
                sx(i)=sc*w+ss*z
                sy(i)=-ss*w+sc*z
   10           continue
           go to 40
   20 continue
           kx=1
           ky=1
           if(incx .lt. 0) kx=1-(n-1)*incx
           if(incy .lt. 0) ky=1-(n-1)*incy
           do 30 i=1,n
                w=sx(kx)
                z=sy(ky)
                sx(kx)=sc*w+ss*z
                sy(ky)=-ss*w+sc*z
                kx=kx+incx
                ky=ky+incy
   30           continue
   40 continue
      return
      end
