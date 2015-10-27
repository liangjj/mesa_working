*deck @(#)prspln.f	1.1 9/8/91
      subroutine prspln(x,cj,cy,nr,lmax)
      implicit integer (a-z)
      common/ io/ inp, iout
      complex *16  cj
      real *8 x, cy
      dimension x(nr), cj(nr,0:lmax), cy(nr,0:lmax)
      do 10 l=0,lmax
         write (iout,100) l
         write (iout,200) (cj(i,l), i=1,nr)
         write (iout,300) l
         write (iout,200) (cy(i,l), i=1,nr)
   10 continue
  100 format(/,5x,'cj for l',1x,i3)
  300 format(/,5x,'cy for l',1x,i3)
  200 format(/,5x,4d15.8)
      return
      end
