*deck @(#)matout.f	5.1  11/6/94
      subroutine matout(xm,mr,mc,nr,nc,iunit)
      implicit real*8(a-h,o-z)
      dimension xm(mr,mc)
c
      istep=5
      do 1 i=1,nc,istep
         ie=min(istep,nc-i+1)
         ie=i+ie-1
         write(iunit,10)(j,j=i,ie)
         do 2 j=1,nr
            write(iunit,11)j,(xm(j,k),k=i,ie)
  2      continue
  1   continue
c
  10  format(5x,5(4x,i5,3x))
  11  format(1x,i3,1x,5(f12.7))
c
c
      return
      end
