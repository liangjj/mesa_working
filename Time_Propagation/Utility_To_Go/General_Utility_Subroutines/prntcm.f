*deck prntcm
      subroutine prntcm(title,a,n,m,ia,ja,iout)
      implicit integer (a-z)
      complex *16 a
      character *80 title
      dimension a(ia,ja)
      write(iout,1) title
    1 format(/,a80,/)
      ibeg=0
      do 10 i=1,m,2
      iend=min(ibeg+2,m)
      ibeg=ibeg+1
      write (iout,20) (ii,ii=ibeg,iend)
   20 format(1x,'col',7x,2(5x,i6,19x))
      do 30 j=1,n
         write (iout,40) (a(j,k),k=ibeg,iend)
   30 continue
      ibeg=iend
   10 continue
   40 format(4x,2(: '(',f15.8,',',f15.8,')' ) )
      return
      end
