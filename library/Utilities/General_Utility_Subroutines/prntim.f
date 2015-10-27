*deck prntim
      subroutine prntim(title,imat,n,m,ia,ja,iout)
      implicit integer (a-z)
      character *80 title
      integer imat
      dimension imat(ia,ja)
      write(iout,1) title
    1 format(a80)
      ibeg=0
      do 10 i=1,m,5
      iend=min(ibeg+5,m)
      ibeg=ibeg+1
      write (iout,20) (ii,ii=ibeg,iend)
   20 format (1x,'col',5(5x,i6,4x))
      do 30 j=1,n
   30 write (iout,40) (imat(j,k),k=ibeg,iend)
      ibeg=iend
   10 continue
   40 format(4x,5i10)
      return
      end
