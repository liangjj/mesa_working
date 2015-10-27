*deck %W%  %G%
      subroutine prntfm(title,a,n,m,ia,ja,iout)
      implicit integer (a-z)
      character *80 title
      real *8 a
      dimension a(ia,ja)
      write(iout,1) title
    1 format(a80)
      ibeg=0
      do 10 i=1,m,6
      iend=min(ibeg+6,m)
      ibeg=ibeg+1
      write (iout,20) (ii,ii=ibeg,iend)
   20 format (1x,'col',6(5x,i6,4x))
      do 30 j=1,n
   30 write (iout,40) (a(j,k),k=ibeg,iend)
      ibeg=iend
   10 continue
   40 format(4x,6e15.8)
      return
      end
