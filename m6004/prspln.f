      subroutine prspln(x,hs,hsder,cj,cy,nr,lmax)
      implicit integer (a-z)
      common/ io/ inp, iout
      complex *16 hs, hsder, cj, cy
      real *8 x
      dimension x(nr), hs(nr,0:lmax), hsder(nr,0:lmax)
      dimension cj(nr,0:lmax), cy(nr,0:lmax)
      do 10 l=0,lmax
         write (iout,100) l
         write (iout,200) (hs(i,l), i=1,nr)
         write (iout,110) l
         write (iout,200) (hsder(i,l), i=1,nr)
         write (iout,120) l
         write (iout,200) (cj(i,l), i=1,nr)
         write (iout,130) l
         write (iout,200) (cy(i,l), i=1,nr)
   10 continue
  100 format(/,5x,'hs for l',1x,i3)
  110 format(/,5x,'hsder for l',1x,i3)
  120 format(/,5x,'cj for l',1x,i3)
  130 format(/,5x,'cy for l',1x,i3)
  200 format(/,5x,4d15.8)
      return
      end
