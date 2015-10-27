*deck @(#)prbes.f	1.1 9/8/91
      subroutine prbes(x,hs,hsder,nr,lmax)
      implicit integer (a-z)
      common/ io/ inp, iout
      complex *16 hs
      real *8 x, hsder
      dimension x(nr), hs(nr,0:lmax), hsder(nr,0:lmax)
      do 10 l=0,lmax
         write (iout,100) l
         write (iout,200) (hs(i,l), i=1,nr)
         write (iout,300) l
         write (iout,200) (hsder(i,l), i=1,nr)
   10 continue
  100 format(/,5x,'hs for l',1x,i3)
  300 format(/,5x,'hsder for l',1x,i3)
  200 format(/,5x,4d15.8)
      return
      end
